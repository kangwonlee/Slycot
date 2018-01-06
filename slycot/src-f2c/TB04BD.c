#line 1 "TB04BD.f"
/* TB04BD.f -- translated by f2c (version 20100827).
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

#line 1 "TB04BD.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;

/* Subroutine */ int tb04bd_(char *jobd, char *order, char *equil, integer *n,
	 integer *m, integer *p, integer *md, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, integer *ign, integer *ldign, integer *
	igd, integer *ldigd, doublereal *gn, doublereal *gd, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen jobd_len, ftnlen order_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, igd_dim1, igd_offset, ign_dim1, ign_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal x, z__[1];
    static integer ia, ib, ic, jj, im, ip, iz, iac, icc;
    static doublereal dij;
    static integer ias, iip, irp, ipm1, ierr, itau;
    static doublereal epsn;
    static integer itau1;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mc01pd_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *), tb01id_(char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb04bx_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, integer *
	    ), tb01zd_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), mc01py_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *);
    static doublereal anorm;
    static logical dijnz, withd;
    static integer ncont;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer jwork, jwork1;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    static logical fndeig, ascend;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    extern /* Subroutine */ int dhseqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
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

/*     To compute the transfer function matrix G of a state-space */
/*     representation (A,B,C,D) of a linear time-invariant multivariable */
/*     system, using the pole-zeros method. Each element of the transfer */
/*     function matrix is returned in a cancelled, minimal form, with */
/*     numerator and denominator polynomials stored either in increasing */
/*     or decreasing order of the powers of the indeterminate. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state-space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed to be a zero matrix. */

/*     ORDER   CHARACTER*1 */
/*             Specifies the order in which the polynomial coefficients */
/*             are stored, as follows: */
/*             = 'I':  Increasing order of powers of the indeterminate; */
/*             = 'D':  Decreasing order of powers of the indeterminate. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the triplet (A,B,C) as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system (A,B,C,D).  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of the system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of the system outputs.  P >= 0. */

/*     MD      (input) INTEGER */
/*             The maximum degree of the polynomials in G, plus 1. An */
/*             upper bound for MD is N+1.  MD >= 1. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, if EQUIL = 'S', the leading N-by-N part of this */
/*             array contains the balanced matrix inv(S)*A*S, as returned */
/*             by SLICOT Library routine TB01ID. */
/*             If EQUIL = 'N', this array is unchanged on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B. */
/*             On exit, the contents of B are destroyed: all elements but */
/*             those in the first row are set to zero. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C. */
/*             On exit, if EQUIL = 'S', the leading P-by-N part of this */
/*             array contains the balanced matrix C*S, as returned by */
/*             SLICOT Library routine TB01ID. */
/*             If EQUIL = 'N', this array is unchanged on exit. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             If JOBD = 'D', the leading P-by-M part of this array must */
/*             contain the matrix D. */
/*             If JOBD = 'Z', the array D is not referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P), if JOBD = 'D'; */
/*             LDD >= 1,        if JOBD = 'Z'. */

/*     IGN     (output) INTEGER array, dimension (LDIGN,M) */
/*             The leading P-by-M part of this array contains the degrees */
/*             of the numerator polynomials in the transfer function */
/*             matrix G. Specifically, the (i,j) element of IGN contains */
/*             the degree of the numerator polynomial of the transfer */
/*             function G(i,j) from the j-th input to the i-th output. */

/*     LDIGN   INTEGER */
/*             The leading dimension of array IGN.  LDIGN >= max(1,P). */

/*     IGD     (output) INTEGER array, dimension (LDIGD,M) */
/*             The leading P-by-M part of this array contains the degrees */
/*             of the denominator polynomials in the transfer function */
/*             matrix G. Specifically, the (i,j) element of IGD contains */
/*             the degree of the denominator polynomial of the transfer */
/*             function G(i,j). */

/*     LDIGD   INTEGER */
/*             The leading dimension of array IGD.  LDIGD >= max(1,P). */

/*     GN      (output) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             This array contains the coefficients of the numerator */
/*             polynomials, Num(i,j), of the transfer function matrix G. */
/*             The polynomials are stored in a column-wise order, i.e., */
/*             Num(1,1), Num(2,1), ..., Num(P,1), Num(1,2), Num(2,2), */
/*             ..., Num(P,2), ..., Num(1,M), Num(2,M), ..., Num(P,M); */
/*             MD memory locations are reserved for each polynomial, */
/*             hence, the (i,j) polynomial is stored starting from the */
/*             location ((j-1)*P+i-1)*MD+1. The coefficients appear in */
/*             increasing or decreasing order of the powers of the */
/*             indeterminate, according to ORDER. */

/*     GD      (output) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             This array contains the coefficients of the denominator */
/*             polynomials, Den(i,j), of the transfer function matrix G. */
/*             The polynomials are stored in the same way as the */
/*             numerator polynomials. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the */
/*             controllability of a single-input system (A,b) or (A',c'), */
/*             where b and c' are columns in B and C' (C transposed). If */
/*             the user sets TOL > 0, then the given value of TOL is used */
/*             as an absolute tolerance; elements with absolute value */
/*             less than TOL are considered neglijible. If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by TOLDEF = N*EPS*MAX( NORM(A), NORM(bc) ) is used */
/*             instead, where EPS is the machine precision (see LAPACK */
/*             Library routine DLAMCH), and bc denotes the currently used */
/*             column in B or C' (see METHOD). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N*(N+P) + */
/*                              MAX( N + MAX( N,P ), N*(2*N+5))) */
/*             If N >= P, N >= 1, the formula above can be written as */
/*             LDWORK >= N*(3*N + P + 5). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the QR algorithm failed to converge when trying to */
/*                   compute the zeros of a transfer function; */
/*             = 2:  the QR algorithm failed to converge when trying to */
/*                   compute the poles of a transfer function. */
/*                   The errors INFO = 1 or 2 are unlikely to appear. */

/*     METHOD */

/*     The routine implements the pole-zero method proposed in [1]. */
/*     This method is based on an algorithm for computing the transfer */
/*     function of a single-input single-output (SISO) system. */
/*     Let (A,b,c,d) be a SISO system. Its transfer function is computed */
/*     as follows: */

/*     1) Find a controllable realization (Ac,bc,cc) of (A,b,c). */
/*     2) Find an observable realization (Ao,bo,co) of (Ac,bc,cc). */
/*     3) Compute the r eigenvalues of Ao (the poles of (Ao,bo,co)). */
/*     4) Compute the zeros of (Ao,bo,co,d). */
/*     5) Compute the gain of (Ao,bo,co,d). */

/*     This algorithm can be implemented using only orthogonal */
/*     transformations [1]. However, for better efficiency, the */
/*     implementation in TB04BD uses one elementary transformation */
/*     in Step 4 and r elementary transformations in Step 5 (to reduce */
/*     an upper Hessenberg matrix to upper triangular form). These */
/*     special elementary transformations are numerically stable */
/*     in practice. */

/*     In the multi-input multi-output (MIMO) case, the algorithm */
/*     computes each element (i,j) of the transfer function matrix G, */
/*     for i = 1 : P, and for j = 1 : M. For efficiency reasons, Step 1 */
/*     is performed once for each value of j (each column of B). The */
/*     matrices Ac and Ao result in Hessenberg form. */

/*     REFERENCES */

/*     [1] Varga, A. and Sima, V. */
/*         Numerically Stable Algorithm for Transfer Function Matrix */
/*         Evaluation. */
/*         Int. J. Control, vol. 33, nr. 6, pp. 1123-1133, 1981. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable in practice and requires about */
/*     20*N**3 floating point operations at most, but usually much less. */

/*     FURTHER COMMENTS */

/*     For maximum efficiency of index calculations, GN and GD are */
/*     implemented as one-dimensional arrays. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2002. */
/*     Partly based on the BIMASC Library routine TSMT1 by A. Varga. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Eigenvalue, state-space representation, transfer function, zeros. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar parameters. */

#line 285 "TB04BD.f"
    /* Parameter adjustments */
#line 285 "TB04BD.f"
    a_dim1 = *lda;
#line 285 "TB04BD.f"
    a_offset = 1 + a_dim1;
#line 285 "TB04BD.f"
    a -= a_offset;
#line 285 "TB04BD.f"
    b_dim1 = *ldb;
#line 285 "TB04BD.f"
    b_offset = 1 + b_dim1;
#line 285 "TB04BD.f"
    b -= b_offset;
#line 285 "TB04BD.f"
    c_dim1 = *ldc;
#line 285 "TB04BD.f"
    c_offset = 1 + c_dim1;
#line 285 "TB04BD.f"
    c__ -= c_offset;
#line 285 "TB04BD.f"
    d_dim1 = *ldd;
#line 285 "TB04BD.f"
    d_offset = 1 + d_dim1;
#line 285 "TB04BD.f"
    d__ -= d_offset;
#line 285 "TB04BD.f"
    ign_dim1 = *ldign;
#line 285 "TB04BD.f"
    ign_offset = 1 + ign_dim1;
#line 285 "TB04BD.f"
    ign -= ign_offset;
#line 285 "TB04BD.f"
    igd_dim1 = *ldigd;
#line 285 "TB04BD.f"
    igd_offset = 1 + igd_dim1;
#line 285 "TB04BD.f"
    igd -= igd_offset;
#line 285 "TB04BD.f"
    --gn;
#line 285 "TB04BD.f"
    --gd;
#line 285 "TB04BD.f"
    --iwork;
#line 285 "TB04BD.f"
    --dwork;
#line 285 "TB04BD.f"

#line 285 "TB04BD.f"
    /* Function Body */
#line 285 "TB04BD.f"
    *info = 0;
#line 286 "TB04BD.f"
    withd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
#line 287 "TB04BD.f"
    ascend = lsame_(order, "I", (ftnlen)1, (ftnlen)1);
#line 288 "TB04BD.f"
    if (! withd && ! lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1)) {
#line 289 "TB04BD.f"
	*info = -1;
#line 290 "TB04BD.f"
    } else if (! ascend && ! lsame_(order, "D", (ftnlen)1, (ftnlen)1)) {
#line 291 "TB04BD.f"
	*info = -2;
#line 292 "TB04BD.f"
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 294 "TB04BD.f"
	*info = -3;
#line 295 "TB04BD.f"
    } else if (*n < 0) {
#line 296 "TB04BD.f"
	*info = -4;
#line 297 "TB04BD.f"
    } else if (*m < 0) {
#line 298 "TB04BD.f"
	*info = -5;
#line 299 "TB04BD.f"
    } else if (*p < 0) {
#line 300 "TB04BD.f"
	*info = -6;
#line 301 "TB04BD.f"
    } else if (*md < 1) {
#line 302 "TB04BD.f"
	*info = -7;
#line 303 "TB04BD.f"
    } else if (*lda < max(1,*n)) {
#line 304 "TB04BD.f"
	*info = -9;
#line 305 "TB04BD.f"
    } else if (*ldb < max(1,*n)) {
#line 306 "TB04BD.f"
	*info = -11;
#line 307 "TB04BD.f"
    } else if (*ldc < max(1,*p)) {
#line 308 "TB04BD.f"
	*info = -13;
#line 309 "TB04BD.f"
    } else if (*ldd < 1 || withd && *ldd < *p) {
#line 310 "TB04BD.f"
	*info = -15;
#line 311 "TB04BD.f"
    } else if (*ldign < max(1,*p)) {
#line 312 "TB04BD.f"
	*info = -17;
#line 313 "TB04BD.f"
    } else if (*ldigd < max(1,*p)) {
#line 314 "TB04BD.f"
	*info = -19;
#line 315 "TB04BD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 315 "TB04BD.f"
	i__3 = *n + max(*n,*p), i__4 = *n * ((*n << 1) + 5);
#line 315 "TB04BD.f"
	i__1 = 1, i__2 = *n * (*n + *p) + max(i__3,i__4);
#line 315 "TB04BD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 318 "TB04BD.f"
	    *info = -25;
#line 319 "TB04BD.f"
	}
#line 319 "TB04BD.f"
    }

#line 321 "TB04BD.f"
    if (*info != 0) {

/*        Error return. */

#line 325 "TB04BD.f"
	i__1 = -(*info);
#line 325 "TB04BD.f"
	xerbla_("TB04BD", &i__1, (ftnlen)6);
#line 326 "TB04BD.f"
	return 0;
#line 327 "TB04BD.f"
    }

/*     Initialize GN and GD to zero. */

#line 331 "TB04BD.f"
    z__[0] = 0.;
#line 332 "TB04BD.f"
    i__1 = *p * *m * *md;
#line 332 "TB04BD.f"
    dcopy_(&i__1, z__, &c__0, &gn[1], &c__1);
#line 333 "TB04BD.f"
    i__1 = *p * *m * *md;
#line 333 "TB04BD.f"
    dcopy_(&i__1, z__, &c__0, &gd[1], &c__1);

/*     Quick return if possible. */

/* Computing MIN */
#line 337 "TB04BD.f"
    i__1 = min(*n,*p);
#line 337 "TB04BD.f"
    if (min(i__1,*m) == 0) {
#line 338 "TB04BD.f"
	if (min(*p,*m) > 0) {
#line 339 "TB04BD.f"
	    k = 1;

#line 341 "TB04BD.f"
	    i__1 = *m;
#line 341 "TB04BD.f"
	    for (j = 1; j <= i__1; ++j) {

#line 343 "TB04BD.f"
		i__2 = *p;
#line 343 "TB04BD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 344 "TB04BD.f"
		    ign[i__ + j * ign_dim1] = 0;
#line 345 "TB04BD.f"
		    igd[i__ + j * igd_dim1] = 0;
#line 346 "TB04BD.f"
		    if (withd) {
#line 346 "TB04BD.f"
			gn[k] = d__[i__ + j * d_dim1];
#line 346 "TB04BD.f"
		    }
#line 348 "TB04BD.f"
		    gd[k] = 1.;
#line 349 "TB04BD.f"
		    k += *md;
#line 350 "TB04BD.f"
/* L10: */
#line 350 "TB04BD.f"
		}

#line 352 "TB04BD.f"
/* L20: */
#line 352 "TB04BD.f"
	    }

#line 354 "TB04BD.f"
	}
#line 355 "TB04BD.f"
	dwork[1] = 1.;
#line 356 "TB04BD.f"
	return 0;
#line 357 "TB04BD.f"
    }

/*     Prepare the computation of the default tolerance. */

#line 361 "TB04BD.f"
    toldef = *tol;
#line 362 "TB04BD.f"
    if (toldef <= 0.) {
#line 363 "TB04BD.f"
	epsn = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
#line 364 "TB04BD.f"
	anorm = dlange_("Frobenius", n, n, &a[a_offset], lda, &dwork[1], (
		ftnlen)9);
#line 365 "TB04BD.f"
    }

/*     Initializations. */

#line 369 "TB04BD.f"
    ia = 1;
#line 370 "TB04BD.f"
    ic = ia + *n * *n;
#line 371 "TB04BD.f"
    itau = ic + *p * *n;
#line 372 "TB04BD.f"
    jwork = itau + *n;
#line 373 "TB04BD.f"
    iac = itau;

#line 375 "TB04BD.f"
    k = 1;
#line 376 "TB04BD.f"
    dij = 0.;

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance.) */

#line 382 "TB04BD.f"
    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(S)*A*S,  B <- inv(S)*B and C <- C*S, where S is a */
/*        diagonal scaling matrix. */
/*        Workspace: need   N. */

#line 389 "TB04BD.f"
	maxred = 100.;
#line 390 "TB04BD.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], &ierr, (ftnlen)3);
#line 392 "TB04BD.f"
    }

/*     Compute the transfer function matrix of the system (A,B,C,D). */

#line 396 "TB04BD.f"
    i__1 = *m;
#line 396 "TB04BD.f"
    for (j = 1; j <= i__1; ++j) {

/*        Save A and C. */
/*        Workspace: need   W1 = N*(N+P). */

#line 401 "TB04BD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)4);
#line 402 "TB04BD.f"
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ic], p, (ftnlen)4);

/*        Remove the uncontrollable part of the system (A,B(J),C). */
/*        Workspace: need   W1+N+MAX(N,P); */
/*                   prefer larger. */

#line 408 "TB04BD.f"
	i__2 = *ldwork - jwork + 1;
#line 408 "TB04BD.f"
	tb01zd_("No Z", n, p, &dwork[ia], n, &b[j * b_dim1 + 1], &dwork[ic], 
		p, &ncont, z__, &c__1, &dwork[itau], tol, &dwork[jwork], &
		i__2, &ierr, (ftnlen)4);
#line 411 "TB04BD.f"
	if (j == 1) {
#line 411 "TB04BD.f"
	    wrkopt = (integer) dwork[jwork] + jwork - 1;
#line 411 "TB04BD.f"
	}

#line 414 "TB04BD.f"
	ib = iac + ncont * ncont;
#line 415 "TB04BD.f"
	icc = ib + ncont;
#line 416 "TB04BD.f"
	itau1 = icc + ncont;
#line 417 "TB04BD.f"
	irp = itau1;
#line 418 "TB04BD.f"
	iip = irp + ncont;
#line 419 "TB04BD.f"
	ias = iip + ncont;
#line 420 "TB04BD.f"
	jwork1 = ias + ncont * ncont;

#line 422 "TB04BD.f"
	i__2 = *p;
#line 422 "TB04BD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 423 "TB04BD.f"
	    if (withd) {
#line 423 "TB04BD.f"
		dij = d__[i__ + j * d_dim1];
#line 423 "TB04BD.f"
	    }
#line 425 "TB04BD.f"
	    if (ncont > 0) {

/*              Form the matrices of the state-space representation of */
/*              the dual system for the controllable part. */
/*              Workspace: need   W2 = W1+N*(N+2). */

#line 431 "TB04BD.f"
		ma02ad_("Full", &ncont, &ncont, &dwork[ia], n, &dwork[iac], &
			ncont, (ftnlen)4);
#line 433 "TB04BD.f"
		dcopy_(&ncont, &b[j * b_dim1 + 1], &c__1, &dwork[ib], &c__1);
#line 434 "TB04BD.f"
		dcopy_(&ncont, &dwork[ic + i__ - 1], p, &dwork[icc], &c__1);

/*              Remove the unobservable part of the system (A,B(J),C(I)). */
/*              Workspace: need   W2+2*N; */
/*                         prefer larger. */

#line 440 "TB04BD.f"
		i__3 = *ldwork - iip + 1;
#line 440 "TB04BD.f"
		tb01zd_("No Z", &ncont, &c__1, &dwork[iac], &ncont, &dwork[
			icc], &dwork[ib], &c__1, &ip, z__, &c__1, &dwork[
			itau1], tol, &dwork[iip], &i__3, &ierr, (ftnlen)4);
#line 444 "TB04BD.f"
		if (i__ == 1) {
/* Computing MAX */
#line 444 "TB04BD.f"
		    i__3 = wrkopt, i__4 = (integer) dwork[iip] + iip - 1;
#line 444 "TB04BD.f"
		    wrkopt = max(i__3,i__4);
#line 444 "TB04BD.f"
		}

#line 447 "TB04BD.f"
		if (ip > 0) {

/*                 Save the state matrix of the minimal part. */
/*                 Workspace: need   W3 = W2+N*(N+2). */

#line 452 "TB04BD.f"
		    dlacpy_("Full", &ip, &ip, &dwork[iac], &ncont, &dwork[ias]
			    , &ip, (ftnlen)4);

/*                 Compute the poles of the transfer function. */
/*                 Workspace: need   W3+N; */
/*                            prefer larger. */

#line 459 "TB04BD.f"
		    i__3 = *ldwork - jwork1 + 1;
#line 459 "TB04BD.f"
		    dhseqr_("Eigenvalues", "No vectors", &ip, &c__1, &ip, &
			    dwork[iac], &ncont, &dwork[irp], &dwork[iip], z__,
			     &c__1, &dwork[jwork1], &i__3, &ierr, (ftnlen)11, 
			    (ftnlen)10);
#line 463 "TB04BD.f"
		    if (ierr != 0) {
#line 464 "TB04BD.f"
			*info = 2;
#line 465 "TB04BD.f"
			return 0;
#line 466 "TB04BD.f"
		    }
/* Computing MAX */
#line 467 "TB04BD.f"
		    i__3 = wrkopt, i__4 = (integer) dwork[jwork1] + jwork1 - 
			    1;
#line 467 "TB04BD.f"
		    wrkopt = max(i__3,i__4);

/*                 Compute the zeros of the transfer function. */

#line 472 "TB04BD.f"
		    ipm1 = ip - 1;
#line 473 "TB04BD.f"
		    dijnz = withd && dij != 0.;
#line 474 "TB04BD.f"
		    fndeig = dijnz || ipm1 > 0;
#line 475 "TB04BD.f"
		    if (! fndeig) {
#line 476 "TB04BD.f"
			iz = 0;
#line 477 "TB04BD.f"
		    } else if (dijnz) {

/*                    Add the contribution due to D(i,j). */
/*                    Note that the matrix whose eigenvalues have to */
/*                    be computed remains in an upper Hessenberg form. */

#line 483 "TB04BD.f"
			iz = ip;
#line 484 "TB04BD.f"
			dlacpy_("Full", &iz, &iz, &dwork[ias], &ip, &dwork[
				iac], &ncont, (ftnlen)4);
#line 486 "TB04BD.f"
			d__1 = -dwork[icc] / dij;
#line 486 "TB04BD.f"
			daxpy_(&iz, &d__1, &dwork[ib], &c__1, &dwork[iac], &
				ncont);
#line 488 "TB04BD.f"
		    } else {
#line 489 "TB04BD.f"
			if (*tol <= 0.) {
/* Computing MAX */
#line 489 "TB04BD.f"
			    d__1 = anorm, d__2 = dlange_("Frobenius", &ip, &
				    c__1, &dwork[ib], &c__1, &dwork[1], (
				    ftnlen)9);
#line 489 "TB04BD.f"
			    toldef = epsn * max(d__1,d__2);
#line 489 "TB04BD.f"
			}

#line 495 "TB04BD.f"
			i__3 = ipm1;
#line 495 "TB04BD.f"
			for (im = 1; im <= i__3; ++im) {
#line 496 "TB04BD.f"
			    if ((d__1 = dwork[ib + im - 1], abs(d__1)) > 
				    toldef) {
#line 496 "TB04BD.f"
				goto L40;
#line 496 "TB04BD.f"
			    }
#line 497 "TB04BD.f"
/* L30: */
#line 497 "TB04BD.f"
			}

#line 499 "TB04BD.f"
			iz = 0;
#line 500 "TB04BD.f"
			goto L50;

#line 502 "TB04BD.f"
L40:

/*                    Restore (part of) the saved state matrix. */

#line 506 "TB04BD.f"
			iz = ip - im;
#line 507 "TB04BD.f"
			dlacpy_("Full", &iz, &iz, &dwork[ias + im * (ip + 1)],
				 &ip, &dwork[iac], &ncont, (ftnlen)4);

/*                    Apply the output injection. */

#line 512 "TB04BD.f"
			d__1 = -dwork[ias + im * (ip + 1) - ip] / dwork[ib + 
				im - 1];
#line 512 "TB04BD.f"
			daxpy_(&iz, &d__1, &dwork[ib + im], &c__1, &dwork[iac]
				, &ncont);
#line 515 "TB04BD.f"
		    }

#line 517 "TB04BD.f"
		    if (fndeig) {

/*                    Find the zeros. */
/*                    Workspace: need   W3+N; */
/*                               prefer larger. */

#line 523 "TB04BD.f"
			i__3 = *ldwork - jwork1 + 1;
#line 523 "TB04BD.f"
			dhseqr_("Eigenvalues", "No vectors", &iz, &c__1, &iz, 
				&dwork[iac], &ncont, &gn[k], &gd[k], z__, &
				c__1, &dwork[jwork1], &i__3, &ierr, (ftnlen)
				11, (ftnlen)10);
#line 527 "TB04BD.f"
			if (ierr != 0) {
#line 528 "TB04BD.f"
			    *info = 1;
#line 529 "TB04BD.f"
			    return 0;
#line 530 "TB04BD.f"
			}
#line 531 "TB04BD.f"
		    }

/*                 Compute the gain. */

#line 535 "TB04BD.f"
L50:
#line 536 "TB04BD.f"
		    if (dijnz) {
#line 537 "TB04BD.f"
			x = dij;
#line 538 "TB04BD.f"
		    } else {
#line 539 "TB04BD.f"
			tb04bx_(&ip, &iz, &dwork[ias], &ip, &dwork[icc], &
				dwork[ib], &dij, &dwork[irp], &dwork[iip], &
				gn[k], &gd[k], &x, &iwork[1]);
#line 542 "TB04BD.f"
		    }

/*                 Form the numerator coefficients in increasing or */
/*                 decreasing powers of the indeterminate. */
/*                 IAS is used here as pointer to the workspace. */

#line 548 "TB04BD.f"
		    if (ascend) {
#line 549 "TB04BD.f"
			mc01pd_(&iz, &gn[k], &gd[k], &dwork[ib], &dwork[ias], 
				&ierr);
#line 551 "TB04BD.f"
		    } else {
#line 552 "TB04BD.f"
			mc01py_(&iz, &gn[k], &gd[k], &dwork[ib], &dwork[ias], 
				&ierr);
#line 554 "TB04BD.f"
		    }
#line 555 "TB04BD.f"
		    jj = k;

#line 557 "TB04BD.f"
		    i__3 = ib + iz;
#line 557 "TB04BD.f"
		    for (l = ib; l <= i__3; ++l) {
#line 558 "TB04BD.f"
			gn[jj] = dwork[l] * x;
#line 559 "TB04BD.f"
			++jj;
#line 560 "TB04BD.f"
/* L60: */
#line 560 "TB04BD.f"
		    }

/*                 Form the denominator coefficients. */

#line 564 "TB04BD.f"
		    if (ascend) {
#line 565 "TB04BD.f"
			mc01pd_(&ip, &dwork[irp], &dwork[iip], &gd[k], &dwork[
				ias], &ierr);
#line 567 "TB04BD.f"
		    } else {
#line 568 "TB04BD.f"
			mc01py_(&ip, &dwork[irp], &dwork[iip], &gd[k], &dwork[
				ias], &ierr);
#line 570 "TB04BD.f"
		    }
#line 571 "TB04BD.f"
		    ign[i__ + j * ign_dim1] = iz;
#line 572 "TB04BD.f"
		    igd[i__ + j * igd_dim1] = ip;
#line 573 "TB04BD.f"
		} else {

/*                 Null element. */

#line 577 "TB04BD.f"
		    ign[i__ + j * ign_dim1] = 0;
#line 578 "TB04BD.f"
		    igd[i__ + j * igd_dim1] = 0;
#line 579 "TB04BD.f"
		    gn[k] = dij;
#line 580 "TB04BD.f"
		    gd[k] = 1.;
#line 581 "TB04BD.f"
		}

#line 583 "TB04BD.f"
	    } else {

/*              Null element. */

#line 587 "TB04BD.f"
		ign[i__ + j * ign_dim1] = 0;
#line 588 "TB04BD.f"
		igd[i__ + j * igd_dim1] = 0;
#line 589 "TB04BD.f"
		gn[k] = dij;
#line 590 "TB04BD.f"
		gd[k] = 1.;
#line 591 "TB04BD.f"
	    }

#line 593 "TB04BD.f"
	    k += *md;
#line 594 "TB04BD.f"
/* L70: */
#line 594 "TB04BD.f"
	}

#line 596 "TB04BD.f"
/* L80: */
#line 596 "TB04BD.f"
    }

#line 598 "TB04BD.f"
    return 0;
/* *** Last line of TB04BD *** */
} /* tb04bd_ */

