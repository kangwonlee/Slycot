#line 1 "AG08BD.f"
/* AG08BD.f -- translated by f2c (version 20100827).
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

#line 1 "AG08BD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static doublereal c_b25 = 0.;
static integer c__0 = 0;

/* Subroutine */ int ag08bd_(char *equil, integer *l, integer *n, integer *m, 
	integer *p, doublereal *a, integer *lda, doublereal *e, integer *lde, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, integer *nfz, integer *nrank, integer *
	niz, integer *dinfz, integer *nkror, integer *ninfe, integer *nkrol, 
	integer *infz, integer *kronr, integer *infe, integer *kronl, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, i0, i1, n2, nb, ii, mm, nn, pp, mu, nu, ipd;
    static doublereal dum[1];
    static integer ldw, itau, numu, kabcd;
    extern /* Subroutine */ int ma02bd_(char *, integer *, integer *, 
	    doublereal *, integer *, ftnlen), ma02cd_(integer *, integer *, 
	    integer *, doublereal *, integer *), tg01ad_(char *, integer *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, ftnlen), tg01fd_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen), ag08by_(logical *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer labcd2;
    static doublereal toler;
    static integer jwork, ldabcd;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer nsinfe;
    static logical lequil;
    static doublereal svlmax;
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

/*     To extract from the system pencil */

/*                       ( A-lambda*E B ) */
/*           S(lambda) = (              ) */
/*                       (      C     D ) */

/*     a regular pencil Af-lambda*Ef which has the finite Smith zeros of */
/*     S(lambda) as generalized eigenvalues. The routine also computes */
/*     the orders of the infinite Smith zeros and determines the singular */
/*     and infinite Kronecker structure of system pencil, i.e., the right */
/*     and left Kronecker indices, and the multiplicities of infinite */
/*     eigenvalues. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to balance the system */
/*             matrix as follows: */
/*             = 'S':  Perform balancing (scaling); */
/*             = 'N':  Do not perform balancing. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The number of rows of matrices A, B, and E.  L >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of matrices A, E, and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of matrix C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A of the system. */
/*             On exit, the leading NFZ-by-NFZ part of this array */
/*             contains the matrix Af of the reduced pencil. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E of the system. */
/*             On exit, the leading NFZ-by-NFZ part of this array */
/*             contains the matrix Ef of the reduced pencil. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B of the system. */
/*             On exit, this matrix does not contain useful information. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0; */
/*             LDB >= 1        if M = 0. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C of the system. */
/*             On exit, this matrix does not contain useful information. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             direct transmission matrix D of the system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     NFZ     (output) INTEGER */
/*             The number of finite zeros. */

/*     NRANK   (output) INTEGER */
/*             The normal rank of the system pencil. */

/*     NIZ     (output) INTEGER */
/*             The number of infinite zeros. */

/*     DINFZ   (output) INTEGER */
/*             The maximal multiplicity of infinite Smith zeros. */

/*     NKROR   (output) INTEGER */
/*             The number of right Kronecker indices. */

/*     NINFE   (output) INTEGER */
/*             The number of elementary infinite blocks. */

/*     NKROL   (output) INTEGER */
/*             The number of left Kronecker indices. */

/*     INFZ    (output) INTEGER array, dimension (N+1) */
/*             The leading DINFZ elements of INFZ contain information */
/*             on the infinite elementary divisors as follows: */
/*             the system has INFZ(i) infinite elementary divisors of */
/*             degree i in the Smith form, where i = 1,2,...,DINFZ. */

/*     KRONR   (output) INTEGER array, dimension (N+M+1) */
/*             The leading NKROR elements of this array contain the */
/*             right Kronecker (column) indices. */

/*     INFE    (output) INTEGER array, dimension (1+MIN(L+P,N+M)) */
/*             The leading NINFE elements of INFE contain the */
/*             multiplicities of infinite eigenvalues. */

/*     KRONL   (output) INTEGER array, dimension (L+P+1) */
/*             The leading NKROL elements of this array contain the */
/*             left Kronecker (row) indices. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance used in rank decisions to determine the */
/*             effective rank, which is defined as the order of the */
/*             largest leading (or trailing) triangular submatrix in the */
/*             QR (or RQ) factorization with column (or row) pivoting */
/*             whose estimated condition number is less than 1/TOL. */
/*             If the user sets TOL <= 0, then default tolerances are */
/*             used instead, as follows: TOLDEF = L*N*EPS in TG01FD */
/*             (to determine the rank of E) and TOLDEF = (L+P)*(N+M)*EPS */
/*             in the rest, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH).  TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension N+max(1,M) */
/*             On output, IWORK(1) contains the normal rank of the */
/*             transfer function matrix. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( 4*(L+N), LDW ), if EQUIL = 'S', */
/*             LDWORK >= LDW,                 if EQUIL = 'N', where */
/*             LDW = max(L+P,M+N)*(M+N) + max(1,5*max(L+P,M+N)). */
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

/*     The routine extracts from the system matrix of a descriptor */
/*     system (A-lambda*E,B,C,D) a regular pencil Af-lambda*Ef which */
/*     has the finite zeros of the system as generalized eigenvalues. */
/*     The procedure has the following main computational steps: */

/*        (a) construct the (L+P)-by-(N+M) system pencil */

/*             S(lambda) = ( B  A )-lambda*( 0  E ); */
/*                         ( D  C )        ( 0  0 ) */

/*        (b) reduce S(lambda) to S1(lambda) with the same finite */
/*            zeros and right Kronecker structure but with E */
/*            upper triangular and nonsingular; */

/*        (c) reduce S1(lambda) to S2(lambda) with the same finite */
/*            zeros and right Kronecker structure but with D of */
/*            full row rank; */

/*        (d) reduce S2(lambda) to S3(lambda) with the same finite zeros */
/*            and with D square invertible; */

/*        (e) perform a unitary transformation on the columns of */

/*            S3(lambda) = (A-lambda*E   B) in order to reduce it to */
/*                         (     C       D) */

/*            (Af-lambda*Ef   X), with Y and Ef square invertible; */
/*            (     0         Y) */

/*        (f) compute the right and left Kronecker indices of the system */
/*            matrix, which together with the multiplicities of the */
/*            finite and infinite eigenvalues constitute the */
/*            complete set of structural invariants under strict */
/*            equivalence transformations of a linear system. */

/*     REFERENCES */

/*     [1] P. Misra, P. Van Dooren and A. Varga. */
/*         Computation of structural invariants of generalized */
/*         state-space systems. */
/*         Automatica, 30, pp. 1921-1936, 1994. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable (see [1]). */

/*     FURTHER COMMENTS */

/*     In order to compute the finite Smith zeros of the system */
/*     explicitly, a call to this routine may be followed by a */
/*     call to the LAPACK Library routines DGEGV or DGGEV. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen, */
/*     May 1999. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Sep. 1999, */
/*     Jan. 2009, Mar. 2009, Apr. 2009. */
/*     A. Varga, DLR Oberpfaffenhofen, Nov. 1999, Feb. 2002, Mar. 2002. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, multivariable */
/*     system, orthogonal transformation, structural invariant. */

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

#line 293 "AG08BD.f"
    /* Parameter adjustments */
#line 293 "AG08BD.f"
    a_dim1 = *lda;
#line 293 "AG08BD.f"
    a_offset = 1 + a_dim1;
#line 293 "AG08BD.f"
    a -= a_offset;
#line 293 "AG08BD.f"
    e_dim1 = *lde;
#line 293 "AG08BD.f"
    e_offset = 1 + e_dim1;
#line 293 "AG08BD.f"
    e -= e_offset;
#line 293 "AG08BD.f"
    b_dim1 = *ldb;
#line 293 "AG08BD.f"
    b_offset = 1 + b_dim1;
#line 293 "AG08BD.f"
    b -= b_offset;
#line 293 "AG08BD.f"
    c_dim1 = *ldc;
#line 293 "AG08BD.f"
    c_offset = 1 + c_dim1;
#line 293 "AG08BD.f"
    c__ -= c_offset;
#line 293 "AG08BD.f"
    d_dim1 = *ldd;
#line 293 "AG08BD.f"
    d_offset = 1 + d_dim1;
#line 293 "AG08BD.f"
    d__ -= d_offset;
#line 293 "AG08BD.f"
    --infz;
#line 293 "AG08BD.f"
    --kronr;
#line 293 "AG08BD.f"
    --infe;
#line 293 "AG08BD.f"
    --kronl;
#line 293 "AG08BD.f"
    --iwork;
#line 293 "AG08BD.f"
    --dwork;
#line 293 "AG08BD.f"

#line 293 "AG08BD.f"
    /* Function Body */
#line 293 "AG08BD.f"
    *info = 0;
/* Computing MAX */
#line 294 "AG08BD.f"
    i__1 = *l + *p, i__2 = *n + *m;
#line 294 "AG08BD.f"
    ldabcd = max(i__1,i__2);
#line 295 "AG08BD.f"
    labcd2 = ldabcd * (*n + *m);
#line 296 "AG08BD.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 297 "AG08BD.f"
    lquery = *ldwork == -1;

/*     Test the input scalar arguments. */

#line 301 "AG08BD.f"
    if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
#line 302 "AG08BD.f"
	*info = -1;
#line 303 "AG08BD.f"
    } else if (*l < 0) {
#line 304 "AG08BD.f"
	*info = -2;
#line 305 "AG08BD.f"
    } else if (*n < 0) {
#line 306 "AG08BD.f"
	*info = -3;
#line 307 "AG08BD.f"
    } else if (*m < 0) {
#line 308 "AG08BD.f"
	*info = -4;
#line 309 "AG08BD.f"
    } else if (*p < 0) {
#line 310 "AG08BD.f"
	*info = -5;
#line 311 "AG08BD.f"
    } else if (*lda < max(1,*l)) {
#line 312 "AG08BD.f"
	*info = -7;
#line 313 "AG08BD.f"
    } else if (*lde < max(1,*l)) {
#line 314 "AG08BD.f"
	*info = -9;
#line 315 "AG08BD.f"
    } else if (*ldb < 1 || *m > 0 && *ldb < *l) {
#line 316 "AG08BD.f"
	*info = -11;
#line 317 "AG08BD.f"
    } else if (*ldc < max(1,*p)) {
#line 318 "AG08BD.f"
	*info = -13;
#line 319 "AG08BD.f"
    } else if (*ldd < max(1,*p)) {
#line 320 "AG08BD.f"
	*info = -15;
#line 321 "AG08BD.f"
    } else if (*tol >= 1.) {
#line 322 "AG08BD.f"
	*info = -27;
#line 323 "AG08BD.f"
    } else {
/* Computing MIN */
#line 324 "AG08BD.f"
	i__1 = *l + *p, i__2 = *m + *n;
#line 324 "AG08BD.f"
	i0 = min(i__1,i__2);
#line 325 "AG08BD.f"
	i1 = min(*l,*n);
#line 326 "AG08BD.f"
	ii = min(*m,*p);
/* Computing MAX */
#line 327 "AG08BD.f"
	i__1 = 1, i__2 = ldabcd * 5;
#line 327 "AG08BD.f"
	ldw = labcd2 + max(i__1,i__2);
#line 328 "AG08BD.f"
	if (lequil) {
/* Computing MAX */
#line 328 "AG08BD.f"
	    i__1 = *l + *n << 2;
#line 328 "AG08BD.f"
	    ldw = max(i__1,ldw);
#line 328 "AG08BD.f"
	}
#line 330 "AG08BD.f"
	if (lquery) {
#line 331 "AG08BD.f"
	    tg01fd_("N", "N", "N", l, n, m, p, &a[a_offset], lda, &e[e_offset]
		    , lde, &b[b_offset], ldb, &c__[c_offset], ldc, dum, &c__1,
		     dum, &c__1, &nn, &n2, tol, &iwork[1], &dwork[1], &c_n1, 
		    info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 334 "AG08BD.f"
	    i__1 = ldw, i__2 = (integer) dwork[1];
#line 334 "AG08BD.f"
	    wrkopt = max(i__1,i__2);
#line 335 "AG08BD.f"
	    svlmax = 0.;
#line 336 "AG08BD.f"
	    i__1 = *m + *n;
#line 336 "AG08BD.f"
	    i__2 = *p + *l;
#line 336 "AG08BD.f"
	    i__3 = ldabcd + i1;
#line 336 "AG08BD.f"
	    ag08by_(&c_true, &i1, &i__1, &i__2, &svlmax, &dwork[1], &i__3, &e[
		    e_offset], lde, &nu, &mu, niz, dinfz, nkrol, &infz[1], &
		    kronl[1], tol, &iwork[1], &dwork[1], &c_n1, info);
/* Computing MAX */
#line 339 "AG08BD.f"
	    i__1 = wrkopt, i__2 = labcd2 + (integer) dwork[1];
#line 339 "AG08BD.f"
	    wrkopt = max(i__1,i__2);
#line 340 "AG08BD.f"
	    i__1 = *m + *n;
#line 340 "AG08BD.f"
	    i__2 = ldabcd + i1;
#line 340 "AG08BD.f"
	    ag08by_(&c_false, &i1, &ii, &i__1, &svlmax, &dwork[1], &i__2, &e[
		    e_offset], lde, &nu, &mu, niz, dinfz, nkrol, &infz[1], &
		    kronl[1], tol, &iwork[1], &dwork[1], &c_n1, info);
/* Computing MAX */
#line 343 "AG08BD.f"
	    i__1 = wrkopt, i__2 = labcd2 + (integer) dwork[1];
#line 343 "AG08BD.f"
	    wrkopt = max(i__1,i__2);
#line 344 "AG08BD.f"
	    i__1 = i1 + ii;
#line 344 "AG08BD.f"
	    nb = ilaenv_(&c__1, "ZGERQF", " ", &ii, &i__1, &c_n1, &c_n1, (
		    ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 345 "AG08BD.f"
	    i__1 = wrkopt, i__2 = labcd2 + ii + ii * nb;
#line 345 "AG08BD.f"
	    wrkopt = max(i__1,i__2);
/* Computing MIN */
#line 346 "AG08BD.f"
	    i__3 = i1 + ii;
#line 346 "AG08BD.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMRQ", "RT", &i1, &i__3, &ii, 
		    &c_n1, (ftnlen)6, (ftnlen)2);
#line 346 "AG08BD.f"
	    nb = min(i__1,i__2);
/* Computing MAX */
#line 348 "AG08BD.f"
	    i__1 = wrkopt, i__2 = labcd2 + ii + i1 * nb;
#line 348 "AG08BD.f"
	    wrkopt = max(i__1,i__2);
#line 349 "AG08BD.f"
	} else if (*ldwork < ldw) {
#line 350 "AG08BD.f"
	    *info = -30;
#line 351 "AG08BD.f"
	}
#line 352 "AG08BD.f"
    }

#line 354 "AG08BD.f"
    if (*info != 0) {

/*        Error return. */

#line 358 "AG08BD.f"
	i__1 = -(*info);
#line 358 "AG08BD.f"
	xerbla_("AG08BD", &i__1, (ftnlen)6);
#line 359 "AG08BD.f"
	return 0;
#line 360 "AG08BD.f"
    } else if (lquery) {
#line 361 "AG08BD.f"
	dwork[1] = (doublereal) wrkopt;
#line 362 "AG08BD.f"
	return 0;
#line 363 "AG08BD.f"
    }

#line 365 "AG08BD.f"
    *niz = 0;
#line 366 "AG08BD.f"
    *nkrol = 0;
#line 367 "AG08BD.f"
    *nkror = 0;

/*     Quick return if possible. */

/* Computing MAX */
#line 371 "AG08BD.f"
    i__1 = max(*l,*n), i__1 = max(i__1,*m);
#line 371 "AG08BD.f"
    if (max(i__1,*p) == 0) {
#line 372 "AG08BD.f"
	*nfz = 0;
#line 373 "AG08BD.f"
	*dinfz = 0;
#line 374 "AG08BD.f"
	*ninfe = 0;
#line 375 "AG08BD.f"
	*nrank = 0;
#line 376 "AG08BD.f"
	iwork[1] = 0;
#line 377 "AG08BD.f"
	dwork[1] = 1.;
#line 378 "AG08BD.f"
	return 0;
#line 379 "AG08BD.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance.) */

#line 385 "AG08BD.f"
    wrkopt = 1;
#line 386 "AG08BD.f"
    kabcd = 1;
#line 387 "AG08BD.f"
    jwork = kabcd + labcd2;

/*     If required, balance the system pencil. */
/*     Workspace: need   4*(L+N). */

#line 392 "AG08BD.f"
    if (lequil) {
#line 393 "AG08BD.f"
	tg01ad_("A", l, n, m, p, &c_b25, &a[a_offset], lda, &e[e_offset], lde,
		 &b[b_offset], ldb, &c__[c_offset], ldc, &dwork[1], &dwork[*l 
		+ 1], &dwork[*l + *n + 1], info, (ftnlen)1);
#line 395 "AG08BD.f"
	wrkopt = *l + *n << 2;
#line 396 "AG08BD.f"
    }
#line 397 "AG08BD.f"
    svlmax = dlange_("Frobenius", l, n, &e[e_offset], lde, &dwork[1], (ftnlen)
	    9);

/*     Reduce the system matrix to QR form, */

/*          ( A11-lambda*E11 A12 B1 ) */
/*          (     A21        A22 B2 ) , */
/*          (     C1         C2  D  ) */

/*     with E11 invertible and upper triangular. */
/*     Real workspace: need   max( 1, N+P, min(L,N)+max(3*N-1,M,L) ); */
/*                     prefer larger. */
/*     Integer workspace: N. */

#line 410 "AG08BD.f"
    tg01fd_("N", "N", "N", l, n, m, p, &a[a_offset], lda, &e[e_offset], lde, &
	    b[b_offset], ldb, &c__[c_offset], ldc, dum, &c__1, dum, &c__1, &
	    nn, &n2, tol, &iwork[1], &dwork[1], ldwork, info, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 413 "AG08BD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 413 "AG08BD.f"
    wrkopt = max(i__1,i__2);

/*     Construct the system pencil */

/*                          MM         NN */
/*                      ( B1 A12 A11-lambda*E11 ) NN */
/*        S1(lambda) =  ( B2 A22      A21       ) L-NN */
/*                      ( D  C2       C1        ) P */

/*     of dimension (L+P)-by-(M+N). */
/*     Workspace: need  LABCD2 = max( L+P, N+M )*( N+M ). */

#line 425 "AG08BD.f"
    n2 = *n - nn;
#line 426 "AG08BD.f"
    mm = *m + n2;
#line 427 "AG08BD.f"
    pp = *p + (*l - nn);
#line 428 "AG08BD.f"
    dlacpy_("Full", l, m, &b[b_offset], ldb, &dwork[kabcd], &ldabcd, (ftnlen)
	    4);
#line 429 "AG08BD.f"
    dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[kabcd + *l], &ldabcd, (
	    ftnlen)4);
#line 430 "AG08BD.f"
    dlacpy_("Full", l, &n2, &a[(nn + 1) * a_dim1 + 1], lda, &dwork[kabcd + 
	    ldabcd * *m], &ldabcd, (ftnlen)4);
#line 432 "AG08BD.f"
    dlacpy_("Full", p, &n2, &c__[(nn + 1) * c_dim1 + 1], ldc, &dwork[kabcd + 
	    ldabcd * *m + *l], &ldabcd, (ftnlen)4);
#line 434 "AG08BD.f"
    dlacpy_("Full", l, &nn, &a[a_offset], lda, &dwork[kabcd + ldabcd * mm], &
	    ldabcd, (ftnlen)4);
#line 436 "AG08BD.f"
    dlacpy_("Full", p, &nn, &c__[c_offset], ldc, &dwork[kabcd + ldabcd * mm + 
	    *l], &ldabcd, (ftnlen)4);

/*     If required, set tolerance. */

#line 441 "AG08BD.f"
    toler = *tol;
#line 442 "AG08BD.f"
    if (toler <= 0.) {
#line 443 "AG08BD.f"
	toler = (doublereal) ((*l + *p) * (*m + *n)) * dlamch_("Precision", (
		ftnlen)9);
#line 444 "AG08BD.f"
    }
/* Computing MAX */
#line 445 "AG08BD.f"
    i__1 = nn + pp;
#line 445 "AG08BD.f"
    i__2 = nn + mm;
#line 445 "AG08BD.f"
    d__1 = svlmax, d__2 = dlange_("Frobenius", &i__1, &i__2, &dwork[kabcd], &
	    ldabcd, &dwork[jwork], (ftnlen)9);
#line 445 "AG08BD.f"
    svlmax = max(d__1,d__2);

/*     Extract the reduced pencil S2(lambda) */

/*             ( Bc  Ac-lambda*Ec ) */
/*             ( Dc      Cc       ) */

/*     having the same finite Smith zeros as the system pencil */
/*     S(lambda) but with Dc, a MU-by-MM full row rank */
/*     left upper trapezoidal matrix, and Ec, an NU-by-NU */
/*     upper triangular nonsingular matrix. */

/*     Real workspace: need   max( min(P+L,M+N)+max(min(L,N),3*(M+N)-1), */
/*                                  5*(P+L), 1 ) + LABCD2; */
/*                     prefer larger. */
/*     Integer workspace: MM, MM <= M+N; PP <= P+L. */

#line 464 "AG08BD.f"
    i__1 = *ldwork - jwork + 1;
#line 464 "AG08BD.f"
    ag08by_(&c_true, &nn, &mm, &pp, &svlmax, &dwork[kabcd], &ldabcd, &e[
	    e_offset], lde, &nu, &mu, niz, dinfz, nkrol, &infz[1], &kronl[1], 
	    &toler, &iwork[1], &dwork[jwork], &i__1, info);

/* Computing MAX */
#line 468 "AG08BD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 468 "AG08BD.f"
    wrkopt = max(i__1,i__2);

/*     Set the number of simple (nondynamic) infinite eigenvalues */
/*     and the normal rank of the system pencil. */

#line 473 "AG08BD.f"
    nsinfe = mu;
#line 474 "AG08BD.f"
    *nrank = nn + mu;

/*     Pertranspose the system. */

/* Computing MAX */
#line 478 "AG08BD.f"
    i__2 = 0, i__3 = nu - 1;
#line 478 "AG08BD.f"
    i__1 = max(i__2,i__3);
/* Computing MAX */
#line 478 "AG08BD.f"
    i__5 = 0, i__6 = nu - 1;
#line 478 "AG08BD.f"
    i__4 = max(i__5,i__6);
#line 478 "AG08BD.f"
    tb01xd_("D", &nu, &mm, &mm, &i__1, &i__4, &dwork[kabcd + ldabcd * mm], &
	    ldabcd, &dwork[kabcd], &ldabcd, &dwork[kabcd + ldabcd * mm + nu], 
	    &ldabcd, &dwork[kabcd + nu], &ldabcd, info, (ftnlen)1);
#line 483 "AG08BD.f"
    i__1 = nu + mm;
#line 483 "AG08BD.f"
    ma02bd_("Right", &i__1, &mm, &dwork[kabcd], &ldabcd, (ftnlen)5);
#line 484 "AG08BD.f"
    i__1 = nu + mm;
#line 484 "AG08BD.f"
    ma02bd_("Left", &mm, &i__1, &dwork[kabcd + nu], &ldabcd, (ftnlen)4);
/* Computing MAX */
#line 485 "AG08BD.f"
    i__2 = 0, i__3 = nu - 1;
#line 485 "AG08BD.f"
    i__1 = max(i__2,i__3);
#line 485 "AG08BD.f"
    ma02cd_(&nu, &c__0, &i__1, &e[e_offset], lde);

#line 487 "AG08BD.f"
    if (mu != mm) {
#line 488 "AG08BD.f"
	nn = nu;
#line 489 "AG08BD.f"
	pp = mm;
#line 490 "AG08BD.f"
	mm = mu;
#line 491 "AG08BD.f"
	kabcd += (pp - mm) * ldabcd;

/*        Extract the reduced pencil S3(lambda), */

/*             ( Br  Ar-lambda*Er ) , */
/*             ( Dr      Cr       ) */

/*        having the same finite Smith zeros as the pencil S(lambda), */
/*        but with Dr, an MU-by-MU invertible upper triangular matrix, */
/*        and Er, an NU-by-NU upper triangular nonsingular matrix. */

/*        Workspace: need   max( 1, 5*(M+N) ) + LABCD2. */
/*                   prefer larger. */
/*        No integer workspace necessary. */

#line 506 "AG08BD.f"
	i__1 = *ldwork - jwork + 1;
#line 506 "AG08BD.f"
	ag08by_(&c_false, &nn, &mm, &pp, &svlmax, &dwork[kabcd], &ldabcd, &e[
		e_offset], lde, &nu, &mu, &i0, &i1, nkror, &iwork[1], &kronr[
		1], &toler, &iwork[1], &dwork[jwork], &i__1, info);

/* Computing MAX */
#line 510 "AG08BD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 510 "AG08BD.f"
	wrkopt = max(i__1,i__2);
#line 511 "AG08BD.f"
    }

#line 513 "AG08BD.f"
    if (nu != 0) {

/*        Perform a unitary transformation on the columns of */
/*                     ( Br Ar-lambda*Er ) */
/*                     ( Dr     Cr       ) */
/*        in order to reduce it to */
/*                     ( *  Af-lambda*Ef ) */
/*                     ( Y       0       ) */
/*        with Y and Ef square invertible. */

/*        Compute Af by reducing  ( Br Ar ) to  ( *  Af ) . */
/*                                ( Dr Cr )     ( Y   0 ) */

#line 526 "AG08BD.f"
	numu = nu + mu;
#line 527 "AG08BD.f"
	ipd = kabcd + nu;
#line 528 "AG08BD.f"
	itau = jwork;
#line 529 "AG08BD.f"
	jwork = itau + mu;

/*        Workspace: need   LABCD2 + 2*min(M,P); */
/*                   prefer LABCD2 + min(M,P) + min(M,P)*NB. */

#line 534 "AG08BD.f"
	i__1 = *ldwork - jwork + 1;
#line 534 "AG08BD.f"
	dtzrzf_(&mu, &numu, &dwork[ipd], &ldabcd, &dwork[itau], &dwork[jwork],
		 &i__1, info);
/* Computing MAX */
#line 536 "AG08BD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 536 "AG08BD.f"
	wrkopt = max(i__1,i__2);

/*        Workspace: need   LABCD2 + min(M,P) + min(L,N); */
/*                   prefer LABCD2 + min(M,P) + min(L,N)*NB. */

#line 541 "AG08BD.f"
	i__1 = *ldwork - jwork + 1;
#line 541 "AG08BD.f"
	dormrz_("Right", "Transpose", &nu, &numu, &mu, &nu, &dwork[ipd], &
		ldabcd, &dwork[itau], &dwork[kabcd], &ldabcd, &dwork[jwork], &
		i__1, info, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 544 "AG08BD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 544 "AG08BD.f"
	wrkopt = max(i__1,i__2);

/*        Save Af. */

#line 548 "AG08BD.f"
	dlacpy_("Full", &nu, &nu, &dwork[kabcd + ldabcd * mu], &ldabcd, &a[
		a_offset], lda, (ftnlen)4);

/*        Compute Ef by applying the saved transformations from previous */
/*        reduction to ( 0  Er ) . */

#line 554 "AG08BD.f"
	dlaset_("Full", &nu, &mu, &c_b25, &c_b25, &dwork[kabcd], &ldabcd, (
		ftnlen)4);
#line 555 "AG08BD.f"
	dlacpy_("Full", &nu, &nu, &e[e_offset], lde, &dwork[kabcd + ldabcd * 
		mu], &ldabcd, (ftnlen)4);

#line 558 "AG08BD.f"
	i__1 = *ldwork - jwork + 1;
#line 558 "AG08BD.f"
	dormrz_("Right", "Transpose", &nu, &numu, &mu, &nu, &dwork[ipd], &
		ldabcd, &dwork[itau], &dwork[kabcd], &ldabcd, &dwork[jwork], &
		i__1, info, (ftnlen)5, (ftnlen)9);

/*        Save Ef. */

#line 564 "AG08BD.f"
	dlacpy_("Full", &nu, &nu, &dwork[kabcd + ldabcd * mu], &ldabcd, &e[
		e_offset], lde, (ftnlen)4);
#line 566 "AG08BD.f"
    }

#line 568 "AG08BD.f"
    *nfz = nu;

/*     Set right Kronecker indices (column indices). */

#line 572 "AG08BD.f"
    i__1 = *nkror;
#line 572 "AG08BD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 573 "AG08BD.f"
	iwork[i__] = kronr[i__];
#line 574 "AG08BD.f"
/* L10: */
#line 574 "AG08BD.f"
    }

#line 576 "AG08BD.f"
    j = 0;
#line 577 "AG08BD.f"
    i__1 = *nkror;
#line 577 "AG08BD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 578 "AG08BD.f"
	i__2 = j + iwork[i__];
#line 578 "AG08BD.f"
	for (ii = j + 1; ii <= i__2; ++ii) {
#line 579 "AG08BD.f"
	    kronr[ii] = i__ - 1;
#line 580 "AG08BD.f"
/* L20: */
#line 580 "AG08BD.f"
	}
#line 581 "AG08BD.f"
	j += iwork[i__];
#line 582 "AG08BD.f"
/* L30: */
#line 582 "AG08BD.f"
    }

#line 584 "AG08BD.f"
    *nkror = j;

/*     Set left Kronecker indices (row indices). */

#line 588 "AG08BD.f"
    i__1 = *nkrol;
#line 588 "AG08BD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 589 "AG08BD.f"
	iwork[i__] = kronl[i__];
#line 590 "AG08BD.f"
/* L40: */
#line 590 "AG08BD.f"
    }

#line 592 "AG08BD.f"
    j = 0;
#line 593 "AG08BD.f"
    i__1 = *nkrol;
#line 593 "AG08BD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 594 "AG08BD.f"
	i__2 = j + iwork[i__];
#line 594 "AG08BD.f"
	for (ii = j + 1; ii <= i__2; ++ii) {
#line 595 "AG08BD.f"
	    kronl[ii] = i__ - 1;
#line 596 "AG08BD.f"
/* L50: */
#line 596 "AG08BD.f"
	}
#line 597 "AG08BD.f"
	j += iwork[i__];
#line 598 "AG08BD.f"
/* L60: */
#line 598 "AG08BD.f"
    }

#line 600 "AG08BD.f"
    *nkrol = j;

/*     Determine the number of simple infinite blocks */
/*     as the difference between the number of infinite blocks */
/*     of order greater than one and the order of Dr. */

#line 606 "AG08BD.f"
    *ninfe = 0;
#line 607 "AG08BD.f"
    i__1 = *dinfz;
#line 607 "AG08BD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 608 "AG08BD.f"
	*ninfe += infz[i__];
#line 609 "AG08BD.f"
/* L70: */
#line 609 "AG08BD.f"
    }
#line 610 "AG08BD.f"
    *ninfe = nsinfe - *ninfe;
#line 611 "AG08BD.f"
    i__1 = *ninfe;
#line 611 "AG08BD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 612 "AG08BD.f"
	infe[i__] = 1;
#line 613 "AG08BD.f"
/* L80: */
#line 613 "AG08BD.f"
    }

/*     Set the structure of infinite eigenvalues. */

#line 617 "AG08BD.f"
    i__1 = *dinfz;
#line 617 "AG08BD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 618 "AG08BD.f"
	i__2 = *ninfe + infz[i__];
#line 618 "AG08BD.f"
	for (ii = *ninfe + 1; ii <= i__2; ++ii) {
#line 619 "AG08BD.f"
	    infe[ii] = i__ + 1;
#line 620 "AG08BD.f"
/* L90: */
#line 620 "AG08BD.f"
	}
#line 621 "AG08BD.f"
	*ninfe += infz[i__];
#line 622 "AG08BD.f"
/* L100: */
#line 622 "AG08BD.f"
    }

#line 624 "AG08BD.f"
    iwork[1] = nsinfe;
#line 625 "AG08BD.f"
    dwork[1] = (doublereal) wrkopt;
#line 626 "AG08BD.f"
    return 0;
/* *** Last line of AG08BD *** */
} /* ag08bd_ */

