#line 1 "AB08ND.f"
/* AB08ND.f -- translated by f2c (version 20100827).
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

#line 1 "AB08ND.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c_n1 = -1;
static integer c__1 = 1;
static doublereal c_b30 = 0.;
static doublereal c_b34 = 1.;

/* Subroutine */ int ab08nd_(char *equil, integer *n, integer *m, integer *p, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, integer *nu, 
	integer *rank, integer *dinfz, integer *nkror, integer *nkrol, 
	integer *infz, integer *kronr, integer *kronl, doublereal *af, 
	integer *ldaf, doublereal *bf, integer *ldbf, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, af_dim1, af_offset, b_dim1, b_offset, bf_dim1, 
	    bf_offset, c_dim1, c_offset, d_dim1, d_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6, i__7;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, i1, nb, ii, mm, nn, pp, ro, mu, nu1, mnu, numu, 
	    numu1;
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ab08nx_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), 
	    dcopy_(integer *, doublereal *, integer *, doublereal *, integer *
	    );
    static integer ninfz;
    static doublereal toler;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical lequil;
    static doublereal thresh, svlmax;
    extern /* Subroutine */ int dormrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static logical lquery;
    extern /* Subroutine */ int dtzrzf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
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

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state dynamics matrix A of the system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input/state matrix B of the system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             state/output matrix C of the system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
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

/*     AF      (output) DOUBLE PRECISION array, dimension */
/*             (LDAF,N+MIN(P,M)) */
/*             The leading NU-by-NU part of this array contains the */
/*             coefficient matrix A  of the reduced pencil. The remainder */
/*                                 f */
/*             of the leading (N+M)-by-(N+MIN(P,M)) part is used as */
/*             internal workspace. */

/*     LDAF    INTEGER */
/*             The leading dimension of array AF.  LDAF >= MAX(1,N+M). */

/*     BF      (output) DOUBLE PRECISION array, dimension (LDBF,N+M) */
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

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX( 1, MIN(P,M) + MAX(3*M-1,N), */
/*                               MIN(P,N) + MAX(3*P-1,N+P,N+M), */
/*                               MIN(M,N) + MAX(3*M-1,N+M) ). */
/*             An upper bound is MAX(s,N) + MAX(3*s-1,N+s), with */
/*             s = MAX(M,P). */
/*             For optimum performance LDWORK should be larger. */

/*             If LDWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             DWORK array, returns this value as the first entry of */
/*             the DWORK array, and no error message related to LDWORK */
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
/*            (      0         Y)              f */

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
/*     Library routine DGGEV with A = A , B = B  and N = NU. */
/*                                     f       f */
/*     If RANK = 0, the routine DGEEV can be used (since B = I). */
/*                                                        f */
/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */
/*     Supersedes Release 2.0 routine AB08BD by F. Svaricek. */

/*     REVISIONS */

/*     Oct. 1997, Feb. 1998, Dec. 2003, March 2004, Jan. 2009, Mar. 2009, */
/*     Apr. 2009. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, multivariable */
/*     system, orthogonal transformation, structural invariant. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 279 "AB08ND.f"
    /* Parameter adjustments */
#line 279 "AB08ND.f"
    a_dim1 = *lda;
#line 279 "AB08ND.f"
    a_offset = 1 + a_dim1;
#line 279 "AB08ND.f"
    a -= a_offset;
#line 279 "AB08ND.f"
    b_dim1 = *ldb;
#line 279 "AB08ND.f"
    b_offset = 1 + b_dim1;
#line 279 "AB08ND.f"
    b -= b_offset;
#line 279 "AB08ND.f"
    c_dim1 = *ldc;
#line 279 "AB08ND.f"
    c_offset = 1 + c_dim1;
#line 279 "AB08ND.f"
    c__ -= c_offset;
#line 279 "AB08ND.f"
    d_dim1 = *ldd;
#line 279 "AB08ND.f"
    d_offset = 1 + d_dim1;
#line 279 "AB08ND.f"
    d__ -= d_offset;
#line 279 "AB08ND.f"
    --infz;
#line 279 "AB08ND.f"
    --kronr;
#line 279 "AB08ND.f"
    --kronl;
#line 279 "AB08ND.f"
    af_dim1 = *ldaf;
#line 279 "AB08ND.f"
    af_offset = 1 + af_dim1;
#line 279 "AB08ND.f"
    af -= af_offset;
#line 279 "AB08ND.f"
    bf_dim1 = *ldbf;
#line 279 "AB08ND.f"
    bf_offset = 1 + bf_dim1;
#line 279 "AB08ND.f"
    bf -= bf_offset;
#line 279 "AB08ND.f"
    --iwork;
#line 279 "AB08ND.f"
    --dwork;
#line 279 "AB08ND.f"

#line 279 "AB08ND.f"
    /* Function Body */
#line 279 "AB08ND.f"
    *info = 0;
#line 280 "AB08ND.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 281 "AB08ND.f"
    lquery = *ldwork == -1;

/*     Test the input scalar arguments. */

#line 285 "AB08ND.f"
    if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
#line 286 "AB08ND.f"
	*info = -1;
#line 287 "AB08ND.f"
    } else if (*n < 0) {
#line 288 "AB08ND.f"
	*info = -2;
#line 289 "AB08ND.f"
    } else if (*m < 0) {
#line 290 "AB08ND.f"
	*info = -3;
#line 291 "AB08ND.f"
    } else if (*p < 0) {
#line 292 "AB08ND.f"
	*info = -4;
#line 293 "AB08ND.f"
    } else if (*lda < max(1,*n)) {
#line 294 "AB08ND.f"
	*info = -6;
#line 295 "AB08ND.f"
    } else if (*ldb < max(1,*n)) {
#line 296 "AB08ND.f"
	*info = -8;
#line 297 "AB08ND.f"
    } else if (*ldc < max(1,*p)) {
#line 298 "AB08ND.f"
	*info = -10;
#line 299 "AB08ND.f"
    } else if (*ldd < max(1,*p)) {
#line 300 "AB08ND.f"
	*info = -12;
#line 301 "AB08ND.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 301 "AB08ND.f"
	i__1 = 1, i__2 = *n + *m;
#line 301 "AB08ND.f"
	if (*ldaf < max(i__1,i__2)) {
#line 302 "AB08ND.f"
	    *info = -22;
#line 303 "AB08ND.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 303 "AB08ND.f"
	    i__1 = 1, i__2 = *n + *p;
#line 303 "AB08ND.f"
	    if (*ldbf < max(i__1,i__2)) {
#line 304 "AB08ND.f"
		*info = -24;
#line 305 "AB08ND.f"
	    } else {
#line 306 "AB08ND.f"
		ii = min(*p,*m);
/* Computing MAX */
/* Computing MAX */
#line 307 "AB08ND.f"
		i__3 = *m * 3 - 1;
/* Computing MAX */
#line 307 "AB08ND.f"
		i__4 = *p * 3 - 1, i__5 = *n + *p, i__4 = max(i__4,i__5), 
			i__5 = *n + *m;
/* Computing MAX */
#line 307 "AB08ND.f"
		i__6 = *m * 3 - 1, i__7 = *n + *m;
#line 307 "AB08ND.f"
		i__1 = ii + max(i__3,*n), i__2 = min(*p,*n) + max(i__4,i__5), 
			i__1 = max(i__1,i__2), i__2 = min(*m,*n) + max(i__6,
			i__7), i__1 = max(i__1,i__2);
#line 307 "AB08ND.f"
		i__ = max(i__1,1);
#line 310 "AB08ND.f"
		if (lquery) {
#line 311 "AB08ND.f"
		    svlmax = 0.;
#line 312 "AB08ND.f"
		    ninfz = 0;
#line 313 "AB08ND.f"
		    ab08nx_(n, m, p, p, &c__0, &svlmax, &bf[bf_offset], ldbf, 
			    &ninfz, &infz[1], &kronl[1], &mu, nu, nkrol, tol, 
			    &iwork[1], &dwork[1], &c_n1, info);
/* Computing MAX */
#line 316 "AB08ND.f"
		    i__1 = i__, i__2 = (integer) dwork[1];
#line 316 "AB08ND.f"
		    wrkopt = max(i__1,i__2);
#line 317 "AB08ND.f"
		    i__1 = *m - ii;
#line 317 "AB08ND.f"
		    ab08nx_(n, &ii, m, &i__1, &ii, &svlmax, &af[af_offset], 
			    ldaf, &ninfz, &infz[1], &kronl[1], &mu, nu, nkrol,
			     tol, &iwork[1], &dwork[1], &c_n1, info);
/* Computing MAX */
#line 320 "AB08ND.f"
		    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 320 "AB08ND.f"
		    wrkopt = max(i__1,i__2);
#line 321 "AB08ND.f"
		    i__1 = *n + ii;
#line 321 "AB08ND.f"
		    nb = ilaenv_(&c__1, "DGERQF", " ", &ii, &i__1, &c_n1, &
			    c_n1, (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 322 "AB08ND.f"
		    i__1 = wrkopt, i__2 = ii + ii * nb;
#line 322 "AB08ND.f"
		    wrkopt = max(i__1,i__2);
/* Computing MIN */
#line 323 "AB08ND.f"
		    i__3 = *n + ii;
#line 323 "AB08ND.f"
		    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RT", n, &i__3,
			     &ii, &c_n1, (ftnlen)6, (ftnlen)2);
#line 323 "AB08ND.f"
		    nb = min(i__1,i__2);
/* Computing MAX */
#line 324 "AB08ND.f"
		    i__1 = wrkopt, i__2 = ii + *n * nb;
#line 324 "AB08ND.f"
		    wrkopt = max(i__1,i__2);
#line 325 "AB08ND.f"
		} else if (*ldwork < i__) {
#line 326 "AB08ND.f"
		    *info = -28;
#line 327 "AB08ND.f"
		}
#line 328 "AB08ND.f"
	    }
#line 328 "AB08ND.f"
	}
#line 328 "AB08ND.f"
    }

#line 330 "AB08ND.f"
    if (*info != 0) {

/*        Error return. */

#line 334 "AB08ND.f"
	i__1 = -(*info);
#line 334 "AB08ND.f"
	xerbla_("AB08ND", &i__1, (ftnlen)6);
#line 335 "AB08ND.f"
	return 0;
#line 336 "AB08ND.f"
    } else if (lquery) {
#line 337 "AB08ND.f"
	dwork[1] = (doublereal) wrkopt;
#line 338 "AB08ND.f"
	return 0;
#line 339 "AB08ND.f"
    }

#line 341 "AB08ND.f"
    *dinfz = 0;
#line 342 "AB08ND.f"
    *nkrol = 0;
#line 343 "AB08ND.f"
    *nkror = 0;

/*     Quick return if possible. */

#line 347 "AB08ND.f"
    if (*n == 0) {
#line 348 "AB08ND.f"
	if (min(*m,*p) == 0) {
#line 349 "AB08ND.f"
	    *nu = 0;
#line 350 "AB08ND.f"
	    *rank = 0;
#line 351 "AB08ND.f"
	    dwork[1] = 1.;
#line 352 "AB08ND.f"
	    return 0;
#line 353 "AB08ND.f"
	}
#line 354 "AB08ND.f"
    }

#line 356 "AB08ND.f"
    mm = *m;
#line 357 "AB08ND.f"
    nn = *n;
#line 358 "AB08ND.f"
    pp = *p;

#line 360 "AB08ND.f"
    i__1 = *n;
#line 360 "AB08ND.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 361 "AB08ND.f"
	infz[i__] = 0;
#line 362 "AB08ND.f"
/* L20: */
#line 362 "AB08ND.f"
    }

#line 364 "AB08ND.f"
    if (*m > 0) {
#line 365 "AB08ND.f"
	i__1 = *n + 1;
#line 365 "AB08ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 366 "AB08ND.f"
	    kronr[i__] = 0;
#line 367 "AB08ND.f"
/* L40: */
#line 367 "AB08ND.f"
	}
#line 368 "AB08ND.f"
    }

#line 370 "AB08ND.f"
    if (*p > 0) {
#line 371 "AB08ND.f"
	i__1 = *n + 1;
#line 371 "AB08ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 372 "AB08ND.f"
	    kronl[i__] = 0;
#line 373 "AB08ND.f"
/* L60: */
#line 373 "AB08ND.f"
	}
#line 374 "AB08ND.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance.) */

#line 380 "AB08ND.f"
    wrkopt = 1;

/*     Construct the compound matrix  ( B  A ), dimension (N+P)-by-(M+N). */
/*                                    ( D  C ) */

#line 385 "AB08ND.f"
    dlacpy_("Full", &nn, &mm, &b[b_offset], ldb, &bf[bf_offset], ldbf, (
	    ftnlen)4);
#line 386 "AB08ND.f"
    if (pp > 0) {
#line 386 "AB08ND.f"
	dlacpy_("Full", &pp, &mm, &d__[d_offset], ldd, &bf[nn + 1 + bf_dim1], 
		ldbf, (ftnlen)4);
#line 386 "AB08ND.f"
    }
#line 388 "AB08ND.f"
    if (nn > 0) {
#line 389 "AB08ND.f"
	dlacpy_("Full", &nn, &nn, &a[a_offset], lda, &bf[(mm + 1) * bf_dim1 + 
		1], ldbf, (ftnlen)4);
#line 390 "AB08ND.f"
	if (pp > 0) {
#line 390 "AB08ND.f"
	    dlacpy_("Full", &pp, &nn, &c__[c_offset], ldc, &bf[nn + 1 + (mm + 
		    1) * bf_dim1], ldbf, (ftnlen)4);
#line 390 "AB08ND.f"
	}
#line 392 "AB08ND.f"
    }

/*     If required, balance the compound matrix (default MAXRED). */
/*     Workspace: need   N. */

#line 397 "AB08ND.f"
    if (lequil && nn > 0 && pp > 0) {
#line 398 "AB08ND.f"
	maxred = 0.;
#line 399 "AB08ND.f"
	tb01id_("A", &nn, &mm, &pp, &maxred, &bf[(mm + 1) * bf_dim1 + 1], 
		ldbf, &bf[bf_offset], ldbf, &bf[nn + 1 + (mm + 1) * bf_dim1], 
		ldbf, &dwork[1], info, (ftnlen)1);
#line 401 "AB08ND.f"
	wrkopt = *n;
#line 402 "AB08ND.f"
    }

/*     If required, set tolerance. */

#line 406 "AB08ND.f"
    thresh = sqrt((doublereal) ((*n + *p) * (*n + *m))) * dlamch_("Precision",
	     (ftnlen)9);
#line 407 "AB08ND.f"
    toler = *tol;
#line 408 "AB08ND.f"
    if (toler < thresh) {
#line 408 "AB08ND.f"
	toler = thresh;
#line 408 "AB08ND.f"
    }
#line 409 "AB08ND.f"
    i__1 = nn + pp;
#line 409 "AB08ND.f"
    i__2 = nn + mm;
#line 409 "AB08ND.f"
    svlmax = dlange_("Frobenius", &i__1, &i__2, &bf[bf_offset], ldbf, &dwork[
	    1], (ftnlen)9);

/*     Reduce this system to one with the same invariant zeros and with */
/*     D upper triangular of full row rank MU (the normal rank of the */
/*     original system). */
/*     Workspace: need   MAX( 1, MIN(P,M) + MAX(3*M-1,N), */
/*                               MIN(P,N) + MAX(3*P-1,N+P,N+M) ); */
/*                prefer larger. */

#line 418 "AB08ND.f"
    ro = pp;
#line 419 "AB08ND.f"
    sigma = 0;
#line 420 "AB08ND.f"
    ninfz = 0;
#line 421 "AB08ND.f"
    ab08nx_(&nn, &mm, &pp, &ro, &sigma, &svlmax, &bf[bf_offset], ldbf, &ninfz,
	     &infz[1], &kronl[1], &mu, nu, nkrol, &toler, &iwork[1], &dwork[1]
	    , ldwork, info);
/* Computing MAX */
#line 424 "AB08ND.f"
    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 424 "AB08ND.f"
    wrkopt = max(i__1,i__2);
#line 425 "AB08ND.f"
    *rank = mu;

/*     Pertranspose the system. */

#line 429 "AB08ND.f"
    numu = *nu + mu;
#line 430 "AB08ND.f"
    if (numu != 0) {
#line 431 "AB08ND.f"
	mnu = mm + *nu;
#line 432 "AB08ND.f"
	numu1 = numu + 1;

#line 434 "AB08ND.f"
	i__1 = numu;
#line 434 "AB08ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 435 "AB08ND.f"
	    dcopy_(&mnu, &bf[i__ + bf_dim1], ldbf, &af[(numu1 - i__) * 
		    af_dim1 + 1], &c_n1);
#line 436 "AB08ND.f"
/* L80: */
#line 436 "AB08ND.f"
	}

#line 438 "AB08ND.f"
	if (mu != mm) {

/*           Here MU < MM and MM > 0 (since MM = 0 implies MU = 0 = MM). */

#line 442 "AB08ND.f"
	    pp = mm;
#line 443 "AB08ND.f"
	    nn = *nu;
#line 444 "AB08ND.f"
	    mm = mu;

/*           Reduce the system to one with the same invariant zeros and */
/*           with D square invertible. */
/*           Workspace: need  MAX( 1, MU + MAX(3*MU-1,N), */
/*                                 MIN(M,N) + MAX(3*M-1,N+M) ); */
/*                prefer larger. Note that MU <= MIN(P,M). */

#line 452 "AB08ND.f"
	    ro = pp - mm;
#line 453 "AB08ND.f"
	    sigma = mm;
#line 454 "AB08ND.f"
	    ab08nx_(&nn, &mm, &pp, &ro, &sigma, &svlmax, &af[af_offset], ldaf,
		     &ninfz, &infz[1], &kronr[1], &mu, nu, nkror, &toler, &
		    iwork[1], &dwork[1], ldwork, info);
/* Computing MAX */
#line 457 "AB08ND.f"
	    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 457 "AB08ND.f"
	    wrkopt = max(i__1,i__2);
#line 458 "AB08ND.f"
	}

#line 460 "AB08ND.f"
	if (*nu != 0) {

/*           Perform a unitary transformation on the columns of */
/*                     ( B  A-lambda*I ) */
/*                     ( D       C     ) */
/*           in order to reduce it to */
/*                     ( X  AF-lambda*BF ) */
/*                     ( Y       0       ) */
/*           with Y and BF square invertible. */

#line 470 "AB08ND.f"
	    dlaset_("Full", nu, &mu, &c_b30, &c_b30, &bf[bf_offset], ldbf, (
		    ftnlen)4);
#line 471 "AB08ND.f"
	    dlaset_("Full", nu, nu, &c_b30, &c_b34, &bf[(mu + 1) * bf_dim1 + 
		    1], ldbf, (ftnlen)4);

#line 473 "AB08ND.f"
	    if (*rank != 0) {
#line 474 "AB08ND.f"
		nu1 = *nu + 1;
#line 475 "AB08ND.f"
		i1 = *nu + mu;

/*              Workspace: need   2*MIN(M,P); */
/*                         prefer MIN(M,P) + MIN(M,P)*NB. */

#line 480 "AB08ND.f"
		i__1 = *ldwork - mu;
#line 480 "AB08ND.f"
		dtzrzf_(&mu, &i1, &af[nu1 + af_dim1], ldaf, &dwork[1], &dwork[
			mu + 1], &i__1, info);
/* Computing MAX */
#line 482 "AB08ND.f"
		i__1 = wrkopt, i__2 = (integer) dwork[mu + 1] + mu;
#line 482 "AB08ND.f"
		wrkopt = max(i__1,i__2);

/*              Workspace: need   MIN(M,P) + N; */
/*                         prefer MIN(M,P) + N*NB. */

#line 487 "AB08ND.f"
		i__1 = *ldwork - mu;
#line 487 "AB08ND.f"
		dormrz_("Right", "Transpose", nu, &i1, &mu, nu, &af[nu1 + 
			af_dim1], ldaf, &dwork[1], &af[af_offset], ldaf, &
			dwork[mu + 1], &i__1, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 490 "AB08ND.f"
		i__1 = wrkopt, i__2 = (integer) dwork[mu + 1] + mu;
#line 490 "AB08ND.f"
		wrkopt = max(i__1,i__2);

#line 492 "AB08ND.f"
		i__1 = *ldwork - mu;
#line 492 "AB08ND.f"
		dormrz_("Right", "Transpose", nu, &i1, &mu, nu, &af[nu1 + 
			af_dim1], ldaf, &dwork[1], &bf[bf_offset], ldbf, &
			dwork[mu + 1], &i__1, info, (ftnlen)5, (ftnlen)9);

#line 496 "AB08ND.f"
	    }

/*           Move AF and BF in the first columns. This assumes that */
/*           DLACPY moves column by column. */

#line 501 "AB08ND.f"
	    dlacpy_("Full", nu, nu, &af[(mu + 1) * af_dim1 + 1], ldaf, &af[
		    af_offset], ldaf, (ftnlen)4);
#line 502 "AB08ND.f"
	    if (*rank != 0) {
#line 502 "AB08ND.f"
		dlacpy_("Full", nu, nu, &bf[(mu + 1) * bf_dim1 + 1], ldbf, &
			bf[bf_offset], ldbf, (ftnlen)4);
#line 502 "AB08ND.f"
	    }

#line 505 "AB08ND.f"
	}
#line 506 "AB08ND.f"
    }

/*     Set right Kronecker indices (column indices). */

#line 510 "AB08ND.f"
    if (*nkror > 0) {
#line 511 "AB08ND.f"
	j = 1;

#line 513 "AB08ND.f"
	i__1 = *n + 1;
#line 513 "AB08ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 515 "AB08ND.f"
	    i__2 = j + kronr[i__] - 1;
#line 515 "AB08ND.f"
	    for (ii = j; ii <= i__2; ++ii) {
#line 516 "AB08ND.f"
		iwork[ii] = i__ - 1;
#line 517 "AB08ND.f"
/* L100: */
#line 517 "AB08ND.f"
	    }

#line 519 "AB08ND.f"
	    j += kronr[i__];
#line 520 "AB08ND.f"
	    kronr[i__] = 0;
#line 521 "AB08ND.f"
/* L120: */
#line 521 "AB08ND.f"
	}

#line 523 "AB08ND.f"
	*nkror = j - 1;

#line 525 "AB08ND.f"
	i__1 = *nkror;
#line 525 "AB08ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 526 "AB08ND.f"
	    kronr[i__] = iwork[i__];
#line 527 "AB08ND.f"
/* L140: */
#line 527 "AB08ND.f"
	}

#line 529 "AB08ND.f"
    }

/*     Set left Kronecker indices (row indices). */

#line 533 "AB08ND.f"
    if (*nkrol > 0) {
#line 534 "AB08ND.f"
	j = 1;

#line 536 "AB08ND.f"
	i__1 = *n + 1;
#line 536 "AB08ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

#line 538 "AB08ND.f"
	    i__2 = j + kronl[i__] - 1;
#line 538 "AB08ND.f"
	    for (ii = j; ii <= i__2; ++ii) {
#line 539 "AB08ND.f"
		iwork[ii] = i__ - 1;
#line 540 "AB08ND.f"
/* L160: */
#line 540 "AB08ND.f"
	    }

#line 542 "AB08ND.f"
	    j += kronl[i__];
#line 543 "AB08ND.f"
	    kronl[i__] = 0;
#line 544 "AB08ND.f"
/* L180: */
#line 544 "AB08ND.f"
	}

#line 546 "AB08ND.f"
	*nkrol = j - 1;

#line 548 "AB08ND.f"
	i__1 = *nkrol;
#line 548 "AB08ND.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 549 "AB08ND.f"
	    kronl[i__] = iwork[i__];
#line 550 "AB08ND.f"
/* L200: */
#line 550 "AB08ND.f"
	}

#line 552 "AB08ND.f"
    }

#line 554 "AB08ND.f"
    if (*n > 0) {
#line 555 "AB08ND.f"
	*dinfz = *n;

#line 557 "AB08ND.f"
L220:
#line 558 "AB08ND.f"
	if (infz[*dinfz] == 0) {
#line 559 "AB08ND.f"
	    --(*dinfz);
#line 560 "AB08ND.f"
	    if (*dinfz > 0) {
#line 560 "AB08ND.f"
		goto L220;
#line 560 "AB08ND.f"
	    }
#line 562 "AB08ND.f"
	}
#line 563 "AB08ND.f"
    }

#line 565 "AB08ND.f"
    dwork[1] = (doublereal) wrkopt;
#line 566 "AB08ND.f"
    return 0;
/* *** Last line of AB08ND *** */
} /* ab08nd_ */

