#line 1 "AB08NZ.f"
/* AB08NZ.f -- translated by f2c (version 20100827).
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

#line 1 "AB08NZ.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static doublecomplex c_b2 = {1.,0.};
static integer c__0 = 0;
static integer c_n1 = -1;
static integer c__1 = 1;

/* Subroutine */ int ab08nz_(char *equil, integer *n, integer *m, integer *p, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *c__, integer *ldc, doublecomplex *d__, integer *ldd, 
	integer *nu, integer *rank, integer *dinfz, integer *nkror, integer *
	nkrol, integer *infz, integer *kronr, integer *kronl, doublecomplex *
	af, integer *ldaf, doublecomplex *bf, integer *ldbf, doublereal *tol, 
	integer *iwork, doublereal *dwork, doublecomplex *zwork, integer *
	lzwork, integer *info, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, bf_dim1, 
	    bf_offset, c_dim1, c_offset, d_dim1, d_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6, i__7;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, i1, nb, ii, mm, nn, pp, ro, mu, nu1, mnu, numu, 
	    numu1, sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01iz_(char *, integer *, integer *, integer 
	    *, doublereal *, doublecomplex *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublereal *, integer *, 
	    ftnlen);
    static integer ninfz;
    static doublereal toler;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), ab8nxz_(integer *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublecomplex *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublecomplex *,
	     integer *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal maxred;
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static logical lequil;
    static doublereal thresh;
    extern /* Subroutine */ int zlacpy_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, ftnlen), 
	    zlaset_(char *, integer *, integer *, doublecomplex *, 
	    doublecomplex *, doublecomplex *, integer *, ftnlen);
    static doublereal svlmax;
    static logical lquery;
    static integer wrkopt;
    extern /* Subroutine */ int zunmrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublecomplex *, integer *, doublecomplex *,
	     doublecomplex *, integer *, doublecomplex *, integer *, integer *
	    , ftnlen, ftnlen), ztzrzf_(integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, doublecomplex *, integer *, integer *)
	    ;


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

/*     To construct for a linear multivariable system described by a */
/*     state-space model (A,B,C,D) a regular pencil (A - lambda*B ) which */
/*                                                    f          f */
/*     has the invariant zeros of the system as generalized eigenvalues. */
/*     The routine also computes the orders of the infinite zeros and the */
/*     right and left Kronecker indices of the system (A,B,C,D). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to balance the compound */
/*             matrix (see METHOD) as follows: */
/*             = 'S':  Perform balancing (scaling); */
/*             = 'N':  Do not perform balancing. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of state variables, i.e., the order of the */
/*             matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     A       (input) COMPLEX*16 array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state dynamics matrix A of the system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) COMPLEX*16 array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input/state matrix B of the system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) COMPLEX*16 array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             state/output matrix C of the system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) COMPLEX*16 array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             direct transmission matrix D of the system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     NU      (output) INTEGER */
/*             The number of (finite) invariant zeros. */

/*     RANK    (output) INTEGER */
/*             The normal rank of the transfer function matrix. */

/*     DINFZ   (output) INTEGER */
/*             The maximum degree of infinite elementary divisors. */

/*     NKROR   (output) INTEGER */
/*             The number of right Kronecker indices. */

/*     NKROL   (output) INTEGER */
/*             The number of left Kronecker indices. */

/*     INFZ    (output) INTEGER array, dimension (N) */
/*             The leading DINFZ elements of INFZ contain information */
/*             on the infinite elementary divisors as follows: */
/*             the system has INFZ(i) infinite elementary divisors */
/*             of degree i, where i = 1,2,...,DINFZ. */

/*     KRONR   (output) INTEGER array, dimension (MAX(N,M)+1) */
/*             The leading NKROR elements of this array contain the */
/*             right Kronecker (column) indices. */

/*     KRONL   (output) INTEGER array, dimension (MAX(N,P)+1) */
/*             The leading NKROL elements of this array contain the */
/*             left Kronecker (row) indices. */

/*     AF      (output) COMPLEX*16 array, dimension (LDAF,N+MIN(P,M)) */
/*             The leading NU-by-NU part of this array contains the */
/*             coefficient matrix A  of the reduced pencil. The remainder */
/*                                 f */
/*             of the leading (N+M)-by-(N+MIN(P,M)) part is used as */
/*             internal workspace. */

/*     LDAF    INTEGER */
/*             The leading dimension of array AF.  LDAF >= MAX(1,N+M). */

/*     BF      (output) COMPLEX*16 array, dimension (LDBF,N+M) */
/*             The leading NU-by-NU part of this array contains the */
/*             coefficient matrix B  of the reduced pencil. The */
/*                                 f */
/*             remainder of the leading (N+P)-by-(N+M) part is used as */
/*             internal workspace. */

/*     LDBF    INTEGER */
/*             The leading dimension of array BF.  LDBF >= MAX(1,N+P). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance used in rank decisions to determine the */
/*             effective rank, which is defined as the order of the */
/*             largest leading (or trailing) triangular submatrix in the */
/*             QR (or RQ) factorization with column (or row) pivoting */
/*             whose estimated condition number is less than 1/TOL. */
/*             If the user sets TOL to be less than SQRT((N+P)*(N+M))*EPS */
/*             then the tolerance is taken as SQRT((N+P)*(N+M))*EPS, */
/*             where EPS is the machine precision (see LAPACK Library */
/*             Routine DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (MAX(M,P)) */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(N,2*MAX(P,M))) */

/*     ZWORK   DOUBLE PRECISION array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the optimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N), */
/*                               MIN(P,N) + MAX(3*P-1,N+P,N+M), */
/*                               MIN(M,N) + MAX(3*M-1,N+M) ). */
/*             An upper bound is MAX(s,N) + MAX(3*s-1,N+s), with */
/*             s = MAX(M,P). */
/*             For optimum performance LZWORK should be larger. */

/*             If LZWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             ZWORK array, returns this value as the first entry of */
/*             the ZWORK array, and no error message related to LZWORK */
/*             is issued by XERBLA. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine extracts from the system matrix of a state-space */
/*     system (A,B,C,D) a regular pencil A - lambda*B  which has the */
/*                                        f          f */
/*     invariant zeros of the system as generalized eigenvalues as */
/*     follows: */

/*        (a) construct the (N+P)-by-(N+M) compound matrix (B  A); */
/*                                                         (D  C) */

/*        (b) reduce the above system to one with the same invariant */
/*            zeros and with D of full row rank; */

/*        (c) pertranspose the system; */

/*        (d) reduce the system to one with the same invariant zeros and */
/*            with D square invertible; */

/*        (e) perform a unitary transformation on the columns of */
/*            (A - lambda*I  B) in order to reduce it to */
/*            (      C       D) */

/*            (A  - lambda*B   X) */
/*            ( f           f   ), with Y and B  square invertible; */
/*            (     0          Y)              f */

/*        (f) compute the right and left Kronecker indices of the system */
/*            (A,B,C,D), which together with the orders of the infinite */
/*            zeros (determined by steps (a) - (e)) constitute the */
/*            complete set of structural invariants under strict */
/*            equivalence transformations of a linear system. */

/*     REFERENCES */

/*     [1] Svaricek, F. */
/*         Computation of the Structural Invariants of Linear */
/*         Multivariable Systems with an Extended Version of */
/*         the Program ZEROS. */
/*         System & Control Letters, 6, pp. 261-266, 1985. */

/*     [2] Emami-Naeini, A. and Van Dooren, P. */
/*         Computation of Zeros of Linear Multivariable Systems. */
/*         Automatica, 18, pp. 415-430, 1982. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable (see [2] and [1]). */

/*     FURTHER COMMENTS */

/*     In order to compute the invariant zeros of the system explicitly, */
/*     a call to this routine may be followed by a call to the LAPACK */
/*     Library routine ZGGEV with A = A , B = B  and N = NU. */
/*                                     f       f */
/*     If RANK = 0, the routine ZGEEV can be used (since B = I). */
/*                                                        f */
/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2009, */
/*     Apr. 2009. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, multivariable */
/*     system, unitary transformation, structural invariant. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 285 "AB08NZ.f"
    /* Parameter adjustments */
#line 285 "AB08NZ.f"
    a_dim1 = *lda;
#line 285 "AB08NZ.f"
    a_offset = 1 + a_dim1;
#line 285 "AB08NZ.f"
    a -= a_offset;
#line 285 "AB08NZ.f"
    b_dim1 = *ldb;
#line 285 "AB08NZ.f"
    b_offset = 1 + b_dim1;
#line 285 "AB08NZ.f"
    b -= b_offset;
#line 285 "AB08NZ.f"
    c_dim1 = *ldc;
#line 285 "AB08NZ.f"
    c_offset = 1 + c_dim1;
#line 285 "AB08NZ.f"
    c__ -= c_offset;
#line 285 "AB08NZ.f"
    d_dim1 = *ldd;
#line 285 "AB08NZ.f"
    d_offset = 1 + d_dim1;
#line 285 "AB08NZ.f"
    d__ -= d_offset;
#line 285 "AB08NZ.f"
    --infz;
#line 285 "AB08NZ.f"
    --kronr;
#line 285 "AB08NZ.f"
    --kronl;
#line 285 "AB08NZ.f"
    af_dim1 = *ldaf;
#line 285 "AB08NZ.f"
    af_offset = 1 + af_dim1;
#line 285 "AB08NZ.f"
    af -= af_offset;
#line 285 "AB08NZ.f"
    bf_dim1 = *ldbf;
#line 285 "AB08NZ.f"
    bf_offset = 1 + bf_dim1;
#line 285 "AB08NZ.f"
    bf -= bf_offset;
#line 285 "AB08NZ.f"
    --iwork;
#line 285 "AB08NZ.f"
    --dwork;
#line 285 "AB08NZ.f"
    --zwork;
#line 285 "AB08NZ.f"

#line 285 "AB08NZ.f"
    /* Function Body */
#line 285 "AB08NZ.f"
    *info = 0;
#line 286 "AB08NZ.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 287 "AB08NZ.f"
    lquery = *lzwork == -1;

/*     Test the input scalar arguments. */

#line 291 "AB08NZ.f"
    if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
#line 292 "AB08NZ.f"
	*info = -1;
#line 293 "AB08NZ.f"
    } else if (*n < 0) {
#line 294 "AB08NZ.f"
	*info = -2;
#line 295 "AB08NZ.f"
    } else if (*m < 0) {
#line 296 "AB08NZ.f"
	*info = -3;
#line 297 "AB08NZ.f"
    } else if (*p < 0) {
#line 298 "AB08NZ.f"
	*info = -4;
#line 299 "AB08NZ.f"
    } else if (*lda < max(1,*n)) {
#line 300 "AB08NZ.f"
	*info = -6;
#line 301 "AB08NZ.f"
    } else if (*ldb < max(1,*n)) {
#line 302 "AB08NZ.f"
	*info = -8;
#line 303 "AB08NZ.f"
    } else if (*ldc < max(1,*p)) {
#line 304 "AB08NZ.f"
	*info = -10;
#line 305 "AB08NZ.f"
    } else if (*ldd < max(1,*p)) {
#line 306 "AB08NZ.f"
	*info = -12;
#line 307 "AB08NZ.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 307 "AB08NZ.f"
	i__1 = 1, i__2 = *n + *m;
#line 307 "AB08NZ.f"
	if (*ldaf < max(i__1,i__2)) {
#line 308 "AB08NZ.f"
	    *info = -22;
#line 309 "AB08NZ.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 309 "AB08NZ.f"
	    i__1 = 1, i__2 = *n + *p;
#line 309 "AB08NZ.f"
	    if (*ldbf < max(i__1,i__2)) {
#line 310 "AB08NZ.f"
		*info = -24;
#line 311 "AB08NZ.f"
	    } else {
#line 312 "AB08NZ.f"
		ii = min(*p,*m);
/* Computing MAX */
/* Computing MAX */
#line 313 "AB08NZ.f"
		i__3 = *m * 3 - 1;
/* Computing MAX */
#line 313 "AB08NZ.f"
		i__4 = *p * 3 - 1, i__5 = *n + *p, i__4 = max(i__4,i__5), 
			i__5 = *n + *m;
/* Computing MAX */
#line 313 "AB08NZ.f"
		i__6 = *m * 3 - 1, i__7 = *n + *m;
#line 313 "AB08NZ.f"
		i__1 = ii + max(i__3,*n), i__2 = min(*p,*n) + max(i__4,i__5), 
			i__1 = max(i__1,i__2), i__2 = min(*m,*n) + max(i__6,
			i__7), i__1 = max(i__1,i__2);
#line 313 "AB08NZ.f"
		i__ = max(i__1,1);
#line 316 "AB08NZ.f"
		if (lquery) {
#line 317 "AB08NZ.f"
		    svlmax = 0.;
#line 318 "AB08NZ.f"
		    ninfz = 0;
#line 319 "AB08NZ.f"
		    ab8nxz_(n, m, p, p, &c__0, &svlmax, &bf[bf_offset], ldbf, 
			    &ninfz, &infz[1], &kronl[1], &mu, nu, nkrol, tol, 
			    &iwork[1], &dwork[1], &zwork[1], &c_n1, info);
/* Computing MAX */
#line 322 "AB08NZ.f"
		    i__1 = i__, i__2 = (integer) zwork[1].r;
#line 322 "AB08NZ.f"
		    wrkopt = max(i__1,i__2);
#line 323 "AB08NZ.f"
		    i__1 = *m - ii;
#line 323 "AB08NZ.f"
		    ab8nxz_(n, &ii, m, &i__1, &ii, &svlmax, &af[af_offset], 
			    ldaf, &ninfz, &infz[1], &kronl[1], &mu, nu, nkrol,
			     tol, &iwork[1], &dwork[1], &zwork[1], &c_n1, 
			    info);
/* Computing MAX */
#line 326 "AB08NZ.f"
		    i__1 = wrkopt, i__2 = (integer) zwork[1].r;
#line 326 "AB08NZ.f"
		    wrkopt = max(i__1,i__2);
#line 327 "AB08NZ.f"
		    i__1 = *n + ii;
#line 327 "AB08NZ.f"
		    nb = ilaenv_(&c__1, "ZGERQF", " ", &ii, &i__1, &c_n1, &
			    c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 328 "AB08NZ.f"
		    i__1 = wrkopt, i__2 = ii + ii * nb;
#line 328 "AB08NZ.f"
		    wrkopt = max(i__1,i__2);
/* Computing MIN */
#line 329 "AB08NZ.f"
		    i__3 = *n + ii;
#line 329 "AB08NZ.f"
		    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", "RC", n, &i__3,
			     &ii, &c_n1, (ftnlen)6, (ftnlen)2);
#line 329 "AB08NZ.f"
		    nb = min(i__1,i__2);
/* Computing MAX */
#line 330 "AB08NZ.f"
		    i__1 = wrkopt, i__2 = ii + *n * nb;
#line 330 "AB08NZ.f"
		    wrkopt = max(i__1,i__2);
#line 331 "AB08NZ.f"
		} else if (*lzwork < i__) {
#line 332 "AB08NZ.f"
		    *info = -29;
#line 333 "AB08NZ.f"
		}
#line 334 "AB08NZ.f"
	    }
#line 334 "AB08NZ.f"
	}
#line 334 "AB08NZ.f"
    }

#line 336 "AB08NZ.f"
    if (*info != 0) {

/*        Error return. */

#line 340 "AB08NZ.f"
	i__1 = -(*info);
#line 340 "AB08NZ.f"
	xerbla_("AB08NZ", &i__1, (ftnlen)6);
#line 341 "AB08NZ.f"
	return 0;
#line 342 "AB08NZ.f"
    } else if (lquery) {
#line 343 "AB08NZ.f"
	zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
#line 344 "AB08NZ.f"
	return 0;
#line 345 "AB08NZ.f"
    }

#line 347 "AB08NZ.f"
    *dinfz = 0;
#line 348 "AB08NZ.f"
    *nkrol = 0;
#line 349 "AB08NZ.f"
    *nkror = 0;

/*     Quick return if possible. */

#line 353 "AB08NZ.f"
    if (*n == 0) {
#line 354 "AB08NZ.f"
	if (min(*m,*p) == 0) {
#line 355 "AB08NZ.f"
	    *nu = 0;
#line 356 "AB08NZ.f"
	    *rank = 0;
#line 357 "AB08NZ.f"
	    zwork[1].r = 1., zwork[1].i = 0.;
#line 358 "AB08NZ.f"
	    return 0;
#line 359 "AB08NZ.f"
	}
#line 360 "AB08NZ.f"
    }

#line 362 "AB08NZ.f"
    mm = *m;
#line 363 "AB08NZ.f"
    nn = *n;
#line 364 "AB08NZ.f"
    pp = *p;

#line 366 "AB08NZ.f"
    i__1 = *n;
#line 366 "AB08NZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 367 "AB08NZ.f"
	infz[i__] = 0;
#line 368 "AB08NZ.f"
/* L20: */
#line 368 "AB08NZ.f"
    }

#line 370 "AB08NZ.f"
    if (*m > 0) {
#line 371 "AB08NZ.f"
	i__1 = *n + 1;
#line 371 "AB08NZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 372 "AB08NZ.f"
	    kronr[i__] = 0;
#line 373 "AB08NZ.f"
/* L40: */
#line 373 "AB08NZ.f"
	}
#line 374 "AB08NZ.f"
    }

#line 376 "AB08NZ.f"
    if (*p > 0) {
#line 377 "AB08NZ.f"
	i__1 = *n + 1;
#line 377 "AB08NZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 378 "AB08NZ.f"
	    kronl[i__] = 0;
#line 379 "AB08NZ.f"
/* L60: */
#line 379 "AB08NZ.f"
	}
#line 380 "AB08NZ.f"
    }

/*     (Note: Comments in the code beginning "CWorkspace:" and */
/*     "RWorkspace:" describe the minimal amount of complex and real */
/*     workspace, respectively, needed at that point in the code, as */
/*     well as the preferred amount for good performance.) */

#line 387 "AB08NZ.f"
    wrkopt = 1;

/*     Construct the compound matrix  ( B  A ), dimension (N+P)-by-(M+N). */
/*                                    ( D  C ) */

#line 392 "AB08NZ.f"
    zlacpy_("Full", &nn, &mm, &b[b_offset], ldb, &bf[bf_offset], ldbf, (
	    ftnlen)4);
#line 393 "AB08NZ.f"
    if (pp > 0) {
#line 393 "AB08NZ.f"
	zlacpy_("Full", &pp, &mm, &d__[d_offset], ldd, &bf[nn + 1 + bf_dim1], 
		ldbf, (ftnlen)4);
#line 393 "AB08NZ.f"
    }
#line 395 "AB08NZ.f"
    if (nn > 0) {
#line 396 "AB08NZ.f"
	zlacpy_("Full", &nn, &nn, &a[a_offset], lda, &bf[(mm + 1) * bf_dim1 + 
		1], ldbf, (ftnlen)4);
#line 397 "AB08NZ.f"
	if (pp > 0) {
#line 397 "AB08NZ.f"
	    zlacpy_("Full", &pp, &nn, &c__[c_offset], ldc, &bf[nn + 1 + (mm + 
		    1) * bf_dim1], ldbf, (ftnlen)4);
#line 397 "AB08NZ.f"
	}
#line 399 "AB08NZ.f"
    }

/*     If required, balance the compound matrix (default MAXRED). */
/*     RWorkspace: need   N. */

#line 404 "AB08NZ.f"
    if (lequil && nn > 0 && pp > 0) {
#line 405 "AB08NZ.f"
	maxred = 0.;
#line 406 "AB08NZ.f"
	tb01iz_("A", &nn, &mm, &pp, &maxred, &bf[(mm + 1) * bf_dim1 + 1], 
		ldbf, &bf[bf_offset], ldbf, &bf[nn + 1 + (mm + 1) * bf_dim1], 
		ldbf, &dwork[1], info, (ftnlen)1);
#line 408 "AB08NZ.f"
    }

/*     If required, set tolerance. */

#line 412 "AB08NZ.f"
    thresh = sqrt((doublereal) ((*n + *p) * (*n + *m))) * dlamch_("Precision",
	     (ftnlen)9);
#line 413 "AB08NZ.f"
    toler = *tol;
#line 414 "AB08NZ.f"
    if (toler < thresh) {
#line 414 "AB08NZ.f"
	toler = thresh;
#line 414 "AB08NZ.f"
    }
#line 415 "AB08NZ.f"
    i__1 = nn + pp;
#line 415 "AB08NZ.f"
    i__2 = nn + mm;
#line 415 "AB08NZ.f"
    svlmax = zlange_("Frobenius", &i__1, &i__2, &bf[bf_offset], ldbf, &dwork[
	    1], (ftnlen)9);

/*     Reduce this system to one with the same invariant zeros and with */
/*     D upper triangular of full row rank MU (the normal rank of the */
/*     original system). */
/*     RWorkspace: need   2*MAX(M,P); */
/*     CWorkspace: need   MAX( 1, MIN(P,M) + MAX(3*M-1,N), */
/*                                MIN(P,N) + MAX(3*P-1,N+P,N+M) ); */
/*                 prefer larger. */

#line 425 "AB08NZ.f"
    ro = pp;
#line 426 "AB08NZ.f"
    sigma = 0;
#line 427 "AB08NZ.f"
    ninfz = 0;
#line 428 "AB08NZ.f"
    ab8nxz_(&nn, &mm, &pp, &ro, &sigma, &svlmax, &bf[bf_offset], ldbf, &ninfz,
	     &infz[1], &kronl[1], &mu, nu, nkrol, &toler, &iwork[1], &dwork[1]
	    , &zwork[1], lzwork, info);
/* Computing MAX */
#line 431 "AB08NZ.f"
    i__1 = wrkopt, i__2 = (integer) zwork[1].r;
#line 431 "AB08NZ.f"
    wrkopt = max(i__1,i__2);
#line 432 "AB08NZ.f"
    *rank = mu;

/*     Pertranspose the system. */

#line 436 "AB08NZ.f"
    numu = *nu + mu;
#line 437 "AB08NZ.f"
    if (numu != 0) {
#line 438 "AB08NZ.f"
	mnu = mm + *nu;
#line 439 "AB08NZ.f"
	numu1 = numu + 1;

#line 441 "AB08NZ.f"
	i__1 = numu;
#line 441 "AB08NZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 442 "AB08NZ.f"
	    zcopy_(&mnu, &bf[i__ + bf_dim1], ldbf, &af[(numu1 - i__) * 
		    af_dim1 + 1], &c_n1);
#line 443 "AB08NZ.f"
/* L80: */
#line 443 "AB08NZ.f"
	}

#line 445 "AB08NZ.f"
	if (mu != mm) {

/*           Here MU < MM and MM > 0 (since MM = 0 implies MU = 0 = MM). */

#line 449 "AB08NZ.f"
	    pp = mm;
#line 450 "AB08NZ.f"
	    nn = *nu;
#line 451 "AB08NZ.f"
	    mm = mu;

/*           Reduce the system to one with the same invariant zeros and */
/*           with D square invertible. */
/*           RWorkspace: need   2*M. */
/*           CWorkspace: need   MAX( 1, MU + MAX(3*MU-1,N), */
/*                                   MIN(M,N) + MAX(3*M-1,N+M) ); */
/*                       prefer larger. Note that MU <= MIN(M,P). */

#line 460 "AB08NZ.f"
	    ro = pp - mm;
#line 461 "AB08NZ.f"
	    sigma = mm;
#line 462 "AB08NZ.f"
	    ab8nxz_(&nn, &mm, &pp, &ro, &sigma, &svlmax, &af[af_offset], ldaf,
		     &ninfz, &infz[1], &kronr[1], &mu, nu, nkror, &toler, &
		    iwork[1], &dwork[1], &zwork[1], lzwork, info);
/* Computing MAX */
#line 465 "AB08NZ.f"
	    i__1 = wrkopt, i__2 = (integer) zwork[1].r;
#line 465 "AB08NZ.f"
	    wrkopt = max(i__1,i__2);
#line 466 "AB08NZ.f"
	}

#line 468 "AB08NZ.f"
	if (*nu != 0) {

/*           Perform a unitary transformation on the columns of */
/*                     ( B  A-lambda*I ) */
/*                     ( D       C     ) */
/*           in order to reduce it to */
/*                     ( X  AF-lambda*BF ) */
/*                     ( Y       0       ) */
/*           with Y and BF square invertible. */

#line 478 "AB08NZ.f"
	    zlaset_("Full", nu, &mu, &c_b1, &c_b1, &bf[bf_offset], ldbf, (
		    ftnlen)4);
#line 479 "AB08NZ.f"
	    zlaset_("Full", nu, nu, &c_b1, &c_b2, &bf[(mu + 1) * bf_dim1 + 1],
		     ldbf, (ftnlen)4);

#line 481 "AB08NZ.f"
	    if (*rank != 0) {
#line 482 "AB08NZ.f"
		nu1 = *nu + 1;
#line 483 "AB08NZ.f"
		i1 = *nu + mu;

/*              CWorkspace: need   2*MIN(M,P); */
/*                          prefer MIN(M,P) + MIN(M,P)*NB. */

#line 488 "AB08NZ.f"
		i__1 = *lzwork - mu;
#line 488 "AB08NZ.f"
		ztzrzf_(&mu, &i1, &af[nu1 + af_dim1], ldaf, &zwork[1], &zwork[
			mu + 1], &i__1, info);
/* Computing MAX */
#line 490 "AB08NZ.f"
		i__3 = mu + 1;
#line 490 "AB08NZ.f"
		i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + mu;
#line 490 "AB08NZ.f"
		wrkopt = max(i__1,i__2);

/*              CWorkspace: need   MIN(M,P) + N; */
/*                          prefer MIN(M,P) + N*NB. */

#line 495 "AB08NZ.f"
		i__1 = *lzwork - mu;
#line 495 "AB08NZ.f"
		zunmrz_("Right", "Conjugate transpose", nu, &i1, &mu, nu, &af[
			nu1 + af_dim1], ldaf, &zwork[1], &af[af_offset], ldaf,
			 &zwork[mu + 1], &i__1, info, (ftnlen)5, (ftnlen)19);
/* Computing MAX */
#line 498 "AB08NZ.f"
		i__3 = mu + 1;
#line 498 "AB08NZ.f"
		i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + mu;
#line 498 "AB08NZ.f"
		wrkopt = max(i__1,i__2);

#line 500 "AB08NZ.f"
		i__1 = *lzwork - mu;
#line 500 "AB08NZ.f"
		zunmrz_("Right", "Conjugate transpose", nu, &i1, &mu, nu, &af[
			nu1 + af_dim1], ldaf, &zwork[1], &bf[bf_offset], ldbf,
			 &zwork[mu + 1], &i__1, info, (ftnlen)5, (ftnlen)19);

#line 504 "AB08NZ.f"
	    }

/*           Move AF and BF in the first columns. This assumes that */
/*           ZLACPY moves column by column. */

#line 509 "AB08NZ.f"
	    zlacpy_("Full", nu, nu, &af[(mu + 1) * af_dim1 + 1], ldaf, &af[
		    af_offset], ldaf, (ftnlen)4);
#line 510 "AB08NZ.f"
	    if (*rank != 0) {
#line 510 "AB08NZ.f"
		zlacpy_("Full", nu, nu, &bf[(mu + 1) * bf_dim1 + 1], ldbf, &
			bf[bf_offset], ldbf, (ftnlen)4);
#line 510 "AB08NZ.f"
	    }

#line 513 "AB08NZ.f"
	}
#line 514 "AB08NZ.f"
    }

/*     Set right Kronecker indices (column indices). */

#line 518 "AB08NZ.f"
    if (*nkror > 0) {
#line 519 "AB08NZ.f"
	j = 1;

#line 521 "AB08NZ.f"
	i__1 = *n + 1;
#line 521 "AB08NZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 523 "AB08NZ.f"
	    i__2 = j + kronr[i__] - 1;
#line 523 "AB08NZ.f"
	    for (ii = j; ii <= i__2; ++ii) {
#line 524 "AB08NZ.f"
		iwork[ii] = i__ - 1;
#line 525 "AB08NZ.f"
/* L100: */
#line 525 "AB08NZ.f"
	    }

#line 527 "AB08NZ.f"
	    j += kronr[i__];
#line 528 "AB08NZ.f"
	    kronr[i__] = 0;
#line 529 "AB08NZ.f"
/* L120: */
#line 529 "AB08NZ.f"
	}

#line 531 "AB08NZ.f"
	*nkror = j - 1;

#line 533 "AB08NZ.f"
	i__1 = *nkror;
#line 533 "AB08NZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 534 "AB08NZ.f"
	    kronr[i__] = iwork[i__];
#line 535 "AB08NZ.f"
/* L140: */
#line 535 "AB08NZ.f"
	}

#line 537 "AB08NZ.f"
    }

/*     Set left Kronecker indices (row indices). */

#line 541 "AB08NZ.f"
    if (*nkrol > 0) {
#line 542 "AB08NZ.f"
	j = 1;

#line 544 "AB08NZ.f"
	i__1 = *n + 1;
#line 544 "AB08NZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 546 "AB08NZ.f"
	    i__2 = j + kronl[i__] - 1;
#line 546 "AB08NZ.f"
	    for (ii = j; ii <= i__2; ++ii) {
#line 547 "AB08NZ.f"
		iwork[ii] = i__ - 1;
#line 548 "AB08NZ.f"
/* L160: */
#line 548 "AB08NZ.f"
	    }

#line 550 "AB08NZ.f"
	    j += kronl[i__];
#line 551 "AB08NZ.f"
	    kronl[i__] = 0;
#line 552 "AB08NZ.f"
/* L180: */
#line 552 "AB08NZ.f"
	}

#line 554 "AB08NZ.f"
	*nkrol = j - 1;

#line 556 "AB08NZ.f"
	i__1 = *nkrol;
#line 556 "AB08NZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 557 "AB08NZ.f"
	    kronl[i__] = iwork[i__];
#line 558 "AB08NZ.f"
/* L200: */
#line 558 "AB08NZ.f"
	}

#line 560 "AB08NZ.f"
    }

#line 562 "AB08NZ.f"
    if (*n > 0) {
#line 563 "AB08NZ.f"
	*dinfz = *n;

#line 565 "AB08NZ.f"
L220:
#line 566 "AB08NZ.f"
	if (infz[*dinfz] == 0) {
#line 567 "AB08NZ.f"
	    --(*dinfz);
#line 568 "AB08NZ.f"
	    if (*dinfz > 0) {
#line 568 "AB08NZ.f"
		goto L220;
#line 568 "AB08NZ.f"
	    }
#line 570 "AB08NZ.f"
	}
#line 571 "AB08NZ.f"
    }

#line 573 "AB08NZ.f"
    zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
#line 574 "AB08NZ.f"
    return 0;
/* *** Last line of AB08NZ *** */
} /* ab08nz_ */

