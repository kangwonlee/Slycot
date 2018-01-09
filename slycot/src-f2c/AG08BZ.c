#line 1 "AG08BZ.f"
/* AG08BZ.f -- translated by f2c (version 20100827).
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

#line 1 "AG08BZ.f"
/* Table of constant values */

static doublecomplex c_b1 = {0.,0.};
static integer c__1 = 1;
static integer c_n1 = -1;
static logical c_true = TRUE_;
static logical c_false = FALSE_;
static doublereal c_b26 = 0.;
static integer c__0 = 0;

/* Subroutine */ int ag08bz_(char *equil, integer *l, integer *n, integer *m, 
	integer *p, doublecomplex *a, integer *lda, doublecomplex *e, integer 
	*lde, doublecomplex *b, integer *ldb, doublecomplex *c__, integer *
	ldc, doublecomplex *d__, integer *ldd, integer *nfz, integer *nrank, 
	integer *niz, integer *dinfz, integer *nkror, integer *ninfe, integer 
	*nkrol, integer *infz, integer *kronr, integer *infe, integer *kronl, 
	doublereal *tol, integer *iwork, doublereal *dwork, doublecomplex *
	zwork, integer *lzwork, integer *info, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j, i0, i1, n2, nb, ii, mm, nn, pp, mu, nu, ipd;
    static doublecomplex dum[1];
    static integer lzw, itau, numu, kabcd;
    extern /* Subroutine */ int ma02bz_(char *, integer *, integer *, 
	    doublecomplex *, integer *, ftnlen), ma02cz_(integer *, integer *,
	     integer *, doublecomplex *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tg01az_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), tg01fz_(char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublecomplex *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer labcd2;
    static doublereal toler;
    extern /* Subroutine */ int tb01xz_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen);
    static integer jwork;
    extern /* Subroutine */ int ag8byz_(logical *, integer *, integer *, 
	    integer *, doublereal *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, integer *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublecomplex *, integer *, integer *);
    static integer ldabcd;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern doublereal zlange_(char *, integer *, integer *, doublecomplex *, 
	    integer *, doublereal *, ftnlen);
    static integer nsinfe;
    static logical lequil;
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

/*     A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the state dynamics matrix A of the system. */
/*             On exit, the leading NFZ-by-NFZ part of this array */
/*             contains the matrix Af of the reduced pencil. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,L). */

/*     E       (input/output) COMPLEX*16 array, dimension (LDE,N) */
/*             On entry, the leading L-by-N part of this array must */
/*             contain the descriptor matrix E of the system. */
/*             On exit, the leading NFZ-by-NFZ part of this array */
/*             contains the matrix Ef of the reduced pencil. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,L). */

/*     B       (input/output) COMPLEX*16 array, dimension (LDB,M) */
/*             On entry, the leading L-by-M part of this array must */
/*             contain the input/state matrix B of the system. */
/*             On exit, this matrix does not contain useful information. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,L) if M > 0; */
/*             LDB >= 1        if M = 0. */

/*     C       (input/output) COMPLEX*16 array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C of the system. */
/*             On exit, this matrix does not contain useful information. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) COMPLEX*16 array, dimension (LDD,M) */
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
/*             used instead, as follows: TOLDEF = L*N*EPS in TG01FZ */
/*             (to determine the rank of E) and TOLDEF = (L+P)*(N+M)*EPS */
/*             in the rest, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH).  TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension N+max(1,M) */
/*             On output, IWORK(1) contains the normal rank of the */
/*             transfer function matrix. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             LDWORK >= max(4*(L+N), 2*max(L+P,M+N))), if EQUIL = 'S', */
/*             LDWORK >= 2*max(L+P,M+N)),               if EQUIL = 'N'. */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the optimal value */
/*             of LZWORK. */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK >= max( max(L+P,M+N)*(M+N) + */
/*                            max(min(L+P,M+N) + max(min(L,N),3*(M+N)-1), */
/*                                3*(L+P), 1)) */
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
/*     call to the LAPACK Library routines ZGEGV or ZGGEV. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, DLR Oberpfaffenhofen, */
/*     May 1999. */
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
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 301 "AG08BZ.f"
    /* Parameter adjustments */
#line 301 "AG08BZ.f"
    a_dim1 = *lda;
#line 301 "AG08BZ.f"
    a_offset = 1 + a_dim1;
#line 301 "AG08BZ.f"
    a -= a_offset;
#line 301 "AG08BZ.f"
    e_dim1 = *lde;
#line 301 "AG08BZ.f"
    e_offset = 1 + e_dim1;
#line 301 "AG08BZ.f"
    e -= e_offset;
#line 301 "AG08BZ.f"
    b_dim1 = *ldb;
#line 301 "AG08BZ.f"
    b_offset = 1 + b_dim1;
#line 301 "AG08BZ.f"
    b -= b_offset;
#line 301 "AG08BZ.f"
    c_dim1 = *ldc;
#line 301 "AG08BZ.f"
    c_offset = 1 + c_dim1;
#line 301 "AG08BZ.f"
    c__ -= c_offset;
#line 301 "AG08BZ.f"
    d_dim1 = *ldd;
#line 301 "AG08BZ.f"
    d_offset = 1 + d_dim1;
#line 301 "AG08BZ.f"
    d__ -= d_offset;
#line 301 "AG08BZ.f"
    --infz;
#line 301 "AG08BZ.f"
    --kronr;
#line 301 "AG08BZ.f"
    --infe;
#line 301 "AG08BZ.f"
    --kronl;
#line 301 "AG08BZ.f"
    --iwork;
#line 301 "AG08BZ.f"
    --dwork;
#line 301 "AG08BZ.f"
    --zwork;
#line 301 "AG08BZ.f"

#line 301 "AG08BZ.f"
    /* Function Body */
#line 301 "AG08BZ.f"
    *info = 0;
/* Computing MAX */
#line 302 "AG08BZ.f"
    i__1 = *l + *p, i__2 = *n + *m;
#line 302 "AG08BZ.f"
    ldabcd = max(i__1,i__2);
#line 303 "AG08BZ.f"
    labcd2 = ldabcd * (*n + *m);
#line 304 "AG08BZ.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 305 "AG08BZ.f"
    lquery = *lzwork == -1;

/*     Test the input scalar arguments. */

#line 309 "AG08BZ.f"
    if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
#line 310 "AG08BZ.f"
	*info = -1;
#line 311 "AG08BZ.f"
    } else if (*l < 0) {
#line 312 "AG08BZ.f"
	*info = -2;
#line 313 "AG08BZ.f"
    } else if (*n < 0) {
#line 314 "AG08BZ.f"
	*info = -3;
#line 315 "AG08BZ.f"
    } else if (*m < 0) {
#line 316 "AG08BZ.f"
	*info = -4;
#line 317 "AG08BZ.f"
    } else if (*p < 0) {
#line 318 "AG08BZ.f"
	*info = -5;
#line 319 "AG08BZ.f"
    } else if (*lda < max(1,*l)) {
#line 320 "AG08BZ.f"
	*info = -7;
#line 321 "AG08BZ.f"
    } else if (*lde < max(1,*l)) {
#line 322 "AG08BZ.f"
	*info = -9;
#line 323 "AG08BZ.f"
    } else if (*ldb < 1 || *m > 0 && *ldb < *l) {
#line 324 "AG08BZ.f"
	*info = -11;
#line 325 "AG08BZ.f"
    } else if (*ldc < max(1,*p)) {
#line 326 "AG08BZ.f"
	*info = -13;
#line 327 "AG08BZ.f"
    } else if (*ldd < max(1,*p)) {
#line 328 "AG08BZ.f"
	*info = -15;
#line 329 "AG08BZ.f"
    } else if (*tol >= 1.) {
#line 330 "AG08BZ.f"
	*info = -27;
#line 331 "AG08BZ.f"
    } else {
/* Computing MIN */
#line 332 "AG08BZ.f"
	i__1 = *l + *p, i__2 = *m + *n;
#line 332 "AG08BZ.f"
	i0 = min(i__1,i__2);
#line 333 "AG08BZ.f"
	i1 = min(*l,*n);
#line 334 "AG08BZ.f"
	ii = min(*m,*p);
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 335 "AG08BZ.f"
	i__5 = i1, i__6 = (*m + *n) * 3 - 1;
#line 335 "AG08BZ.f"
	i__3 = i0 + max(i__5,i__6), i__4 = (*l + *p) * 3;
#line 335 "AG08BZ.f"
	i__1 = 1, i__2 = labcd2 + max(i__3,i__4);
#line 335 "AG08BZ.f"
	lzw = max(i__1,i__2);
#line 337 "AG08BZ.f"
	if (lquery) {
#line 338 "AG08BZ.f"
	    tg01fz_("N", "N", "N", l, n, m, p, &a[a_offset], lda, &e[e_offset]
		    , lde, &b[b_offset], ldb, &c__[c_offset], ldc, dum, &c__1,
		     dum, &c__1, &nn, &n2, tol, &iwork[1], &dwork[1], &zwork[
		    1], &c_n1, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 341 "AG08BZ.f"
	    i__1 = lzw, i__2 = (integer) zwork[1].r;
#line 341 "AG08BZ.f"
	    wrkopt = max(i__1,i__2);
#line 342 "AG08BZ.f"
	    svlmax = 0.;
#line 343 "AG08BZ.f"
	    i__1 = *m + *n;
#line 343 "AG08BZ.f"
	    i__2 = *p + *l;
#line 343 "AG08BZ.f"
	    i__3 = ldabcd + i1;
#line 343 "AG08BZ.f"
	    ag8byz_(&c_true, &i1, &i__1, &i__2, &svlmax, &zwork[1], &i__3, &e[
		    e_offset], lde, &nu, &mu, niz, dinfz, nkrol, &infz[1], &
		    kronl[1], tol, &iwork[1], &dwork[1], &zwork[1], &c_n1, 
		    info);
/* Computing MAX */
#line 346 "AG08BZ.f"
	    i__1 = wrkopt, i__2 = labcd2 + (integer) zwork[1].r;
#line 346 "AG08BZ.f"
	    wrkopt = max(i__1,i__2);
#line 347 "AG08BZ.f"
	    i__1 = *m + *n;
#line 347 "AG08BZ.f"
	    i__2 = ldabcd + i1;
#line 347 "AG08BZ.f"
	    ag8byz_(&c_false, &i1, &ii, &i__1, &svlmax, &zwork[1], &i__2, &e[
		    e_offset], lde, &nu, &mu, niz, dinfz, nkrol, &infz[1], &
		    kronl[1], tol, &iwork[1], &dwork[1], &zwork[1], &c_n1, 
		    info);
/* Computing MAX */
#line 350 "AG08BZ.f"
	    i__1 = wrkopt, i__2 = labcd2 + (integer) zwork[1].r;
#line 350 "AG08BZ.f"
	    wrkopt = max(i__1,i__2);
#line 351 "AG08BZ.f"
	    i__1 = i1 + ii;
#line 351 "AG08BZ.f"
	    nb = ilaenv_(&c__1, "ZGERQF", " ", &ii, &i__1, &c_n1, &c_n1, (
		    ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 352 "AG08BZ.f"
	    i__1 = wrkopt, i__2 = labcd2 + ii + ii * nb;
#line 352 "AG08BZ.f"
	    wrkopt = max(i__1,i__2);
/* Computing MIN */
#line 353 "AG08BZ.f"
	    i__3 = i1 + ii;
#line 353 "AG08BZ.f"
	    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMRQ", "RC", &i1, &i__3, &ii, 
		    &c_n1, (ftnlen)6, (ftnlen)2);
#line 353 "AG08BZ.f"
	    nb = min(i__1,i__2);
/* Computing MAX */
#line 355 "AG08BZ.f"
	    i__1 = wrkopt, i__2 = labcd2 + ii + i1 * nb;
#line 355 "AG08BZ.f"
	    wrkopt = max(i__1,i__2);
#line 356 "AG08BZ.f"
	} else if (*lzwork < lzw) {
#line 357 "AG08BZ.f"
	    *info = -31;
#line 358 "AG08BZ.f"
	}
#line 359 "AG08BZ.f"
    }

#line 361 "AG08BZ.f"
    if (*info != 0) {

/*        Error return. */

#line 365 "AG08BZ.f"
	i__1 = -(*info);
#line 365 "AG08BZ.f"
	xerbla_("AG08BZ", &i__1, (ftnlen)6);
#line 366 "AG08BZ.f"
	return 0;
#line 367 "AG08BZ.f"
    } else if (lquery) {
#line 368 "AG08BZ.f"
	zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
#line 369 "AG08BZ.f"
	return 0;
#line 370 "AG08BZ.f"
    }

#line 372 "AG08BZ.f"
    *niz = 0;
#line 373 "AG08BZ.f"
    *nkrol = 0;
#line 374 "AG08BZ.f"
    *nkror = 0;

/*     Quick return if possible. */

/* Computing MAX */
#line 378 "AG08BZ.f"
    i__1 = max(*l,*n), i__1 = max(i__1,*m);
#line 378 "AG08BZ.f"
    if (max(i__1,*p) == 0) {
#line 379 "AG08BZ.f"
	*nfz = 0;
#line 380 "AG08BZ.f"
	*dinfz = 0;
#line 381 "AG08BZ.f"
	*ninfe = 0;
#line 382 "AG08BZ.f"
	*nrank = 0;
#line 383 "AG08BZ.f"
	iwork[1] = 0;
#line 384 "AG08BZ.f"
	zwork[1].r = 1., zwork[1].i = 0.;
#line 385 "AG08BZ.f"
	return 0;
#line 386 "AG08BZ.f"
    }

/*     (Note: Comments in the code beginning "CWorkspace:", "RWorkspace:" */
/*     and "IWorkspace:" describe the minimal amount of complex, real and */
/*     integer workspace, respectively, needed at that point in the code, */
/*     as well as the preferred amount for good performance.) */

#line 393 "AG08BZ.f"
    wrkopt = 1;
#line 394 "AG08BZ.f"
    kabcd = 1;
#line 395 "AG08BZ.f"
    jwork = kabcd + labcd2;

/*     If required, balance the system pencil. */
/*     RWorkspace: need   4*(L+N). */

#line 400 "AG08BZ.f"
    if (lequil) {
#line 401 "AG08BZ.f"
	tg01az_("A", l, n, m, p, &c_b26, &a[a_offset], lda, &e[e_offset], lde,
		 &b[b_offset], ldb, &c__[c_offset], ldc, &dwork[1], &dwork[*l 
		+ 1], &dwork[*l + *n + 1], info, (ftnlen)1);
#line 403 "AG08BZ.f"
    }
#line 404 "AG08BZ.f"
    svlmax = zlange_("Frobenius", l, n, &e[e_offset], lde, &dwork[1], (ftnlen)
	    9);

/*     Reduce the system matrix to QR form, */

/*          ( A11-lambda*E11 A12 B1 ) */
/*          (     A21        A22 B2 ) , */
/*          (     C1         C2  D  ) */

/*     with E11 invertible and upper triangular. */
/*     IWorkspace: need   N. */
/*     RWorkspace: need   2*N. */
/*     CWorkspace: need   max( 1, N+P, min(L,N)+max(3*N-1,M,L) ); */
/*                 prefer larger. */

#line 418 "AG08BZ.f"
    tg01fz_("N", "N", "N", l, n, m, p, &a[a_offset], lda, &e[e_offset], lde, &
	    b[b_offset], ldb, &c__[c_offset], ldc, dum, &c__1, dum, &c__1, &
	    nn, &n2, tol, &iwork[1], &dwork[1], &zwork[1], lzwork, info, (
	    ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 421 "AG08BZ.f"
    i__1 = wrkopt, i__2 = (integer) zwork[1].r;
#line 421 "AG08BZ.f"
    wrkopt = max(i__1,i__2);

/*     Construct the system pencil */

/*                          MM         NN */
/*                      ( B1 A12 A11-lambda*E11 ) NN */
/*        S1(lambda) =  ( B2 A22      A21       ) L-NN */
/*                      ( D  C2       C1        ) P */

/*     of dimension (L+P)-by-(M+N). */
/*     CWorkspace: need  LABCD2 = max( L+P, N+M )*( N+M ). */

#line 433 "AG08BZ.f"
    n2 = *n - nn;
#line 434 "AG08BZ.f"
    mm = *m + n2;
#line 435 "AG08BZ.f"
    pp = *p + (*l - nn);
#line 436 "AG08BZ.f"
    zlacpy_("Full", l, m, &b[b_offset], ldb, &zwork[kabcd], &ldabcd, (ftnlen)
	    4);
#line 437 "AG08BZ.f"
    zlacpy_("Full", p, m, &d__[d_offset], ldd, &zwork[kabcd + *l], &ldabcd, (
	    ftnlen)4);
#line 438 "AG08BZ.f"
    zlacpy_("Full", l, &n2, &a[(nn + 1) * a_dim1 + 1], lda, &zwork[kabcd + 
	    ldabcd * *m], &ldabcd, (ftnlen)4);
#line 440 "AG08BZ.f"
    zlacpy_("Full", p, &n2, &c__[(nn + 1) * c_dim1 + 1], ldc, &zwork[kabcd + 
	    ldabcd * *m + *l], &ldabcd, (ftnlen)4);
#line 442 "AG08BZ.f"
    zlacpy_("Full", l, &nn, &a[a_offset], lda, &zwork[kabcd + ldabcd * mm], &
	    ldabcd, (ftnlen)4);
#line 444 "AG08BZ.f"
    zlacpy_("Full", p, &nn, &c__[c_offset], ldc, &zwork[kabcd + ldabcd * mm + 
	    *l], &ldabcd, (ftnlen)4);

/*     If required, set tolerance. */

#line 449 "AG08BZ.f"
    toler = *tol;
#line 450 "AG08BZ.f"
    if (toler <= 0.) {
#line 451 "AG08BZ.f"
	toler = (doublereal) ((*l + *p) * (*m + *n)) * dlamch_("Precision", (
		ftnlen)9);
#line 452 "AG08BZ.f"
    }
/* Computing MAX */
#line 453 "AG08BZ.f"
    i__1 = nn + pp;
#line 453 "AG08BZ.f"
    i__2 = nn + mm;
#line 453 "AG08BZ.f"
    d__1 = svlmax, d__2 = zlange_("Frobenius", &i__1, &i__2, &zwork[kabcd], &
	    ldabcd, &dwork[1], (ftnlen)9);
#line 453 "AG08BZ.f"
    svlmax = max(d__1,d__2);

/*     Extract the reduced pencil S2(lambda) */

/*             ( Bc  Ac-lambda*Ec ) */
/*             ( Dc      Cc       ) */

/*     having the same finite Smith zeros as the system pencil */
/*     S(lambda) but with Dc, a MU-by-MM full row rank */
/*     left upper trapezoidal matrix, and Ec, an NU-by-NU */
/*     upper triangular nonsingular matrix. */

/*     IWorkspace: need   MM,            MM <= M+N; */
/*     RWorkspace: need   2*max(MM,PP);  PP <= P+L; */
/*     CWorkspace: need   max( min(P+L,M+N)+max(min(L,N),3*(M+N)-1), */
/*                              3*(P+L), 1 ) + LABCD2; */
/*                 prefer larger. */

#line 473 "AG08BZ.f"
    i__1 = *lzwork - jwork + 1;
#line 473 "AG08BZ.f"
    ag8byz_(&c_true, &nn, &mm, &pp, &svlmax, &zwork[kabcd], &ldabcd, &e[
	    e_offset], lde, &nu, &mu, niz, dinfz, nkrol, &infz[1], &kronl[1], 
	    &toler, &iwork[1], &dwork[1], &zwork[jwork], &i__1, info);

/* Computing MAX */
#line 478 "AG08BZ.f"
    i__3 = jwork;
#line 478 "AG08BZ.f"
    i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + jwork - 1;
#line 478 "AG08BZ.f"
    wrkopt = max(i__1,i__2);

/*     Set the number of simple (nondynamic) infinite eigenvalues */
/*     and the normal rank of the system pencil. */

#line 483 "AG08BZ.f"
    nsinfe = mu;
#line 484 "AG08BZ.f"
    *nrank = nn + mu;

/*     Pertranspose the system. */

/* Computing MAX */
#line 488 "AG08BZ.f"
    i__2 = 0, i__3 = nu - 1;
#line 488 "AG08BZ.f"
    i__1 = max(i__2,i__3);
/* Computing MAX */
#line 488 "AG08BZ.f"
    i__5 = 0, i__6 = nu - 1;
#line 488 "AG08BZ.f"
    i__4 = max(i__5,i__6);
#line 488 "AG08BZ.f"
    tb01xz_("D", &nu, &mm, &mm, &i__1, &i__4, &zwork[kabcd + ldabcd * mm], &
	    ldabcd, &zwork[kabcd], &ldabcd, &zwork[kabcd + ldabcd * mm + nu], 
	    &ldabcd, &zwork[kabcd + nu], &ldabcd, info, (ftnlen)1);
#line 493 "AG08BZ.f"
    i__1 = nu + mm;
#line 493 "AG08BZ.f"
    ma02bz_("Right", &i__1, &mm, &zwork[kabcd], &ldabcd, (ftnlen)5);
#line 494 "AG08BZ.f"
    i__1 = nu + mm;
#line 494 "AG08BZ.f"
    ma02bz_("Left", &mm, &i__1, &zwork[kabcd + nu], &ldabcd, (ftnlen)4);
/* Computing MAX */
#line 495 "AG08BZ.f"
    i__2 = 0, i__3 = nu - 1;
#line 495 "AG08BZ.f"
    i__1 = max(i__2,i__3);
#line 495 "AG08BZ.f"
    ma02cz_(&nu, &c__0, &i__1, &e[e_offset], lde);

#line 497 "AG08BZ.f"
    if (mu != mm) {
#line 498 "AG08BZ.f"
	nn = nu;
#line 499 "AG08BZ.f"
	pp = mm;
#line 500 "AG08BZ.f"
	mm = mu;
#line 501 "AG08BZ.f"
	kabcd += (pp - mm) * ldabcd;

/*        Extract the reduced pencil S3(lambda), */

/*             ( Br  Ar-lambda*Er ) , */
/*             ( Dr      Cr       ) */

/*        having the same finite Smith zeros as the pencil S(lambda), */
/*        but with Dr, an MU-by-MU invertible upper triangular matrix, */
/*        and Er, an NU-by-NU upper triangular nonsingular matrix. */

/*        IWorkspace: need   0; */
/*        RWorkspace: need   2*(M+N); */
/*        CWorkspace: need   max( 1, 3*(M+N) ) + LABCD2. */
/*                    prefer larger. */

#line 517 "AG08BZ.f"
	i__1 = *lzwork - jwork + 1;
#line 517 "AG08BZ.f"
	ag8byz_(&c_false, &nn, &mm, &pp, &svlmax, &zwork[kabcd], &ldabcd, &e[
		e_offset], lde, &nu, &mu, &i0, &i1, nkror, &iwork[1], &kronr[
		1], &toler, &iwork[1], &dwork[1], &zwork[jwork], &i__1, info);

/* Computing MAX */
#line 522 "AG08BZ.f"
	i__3 = jwork;
#line 522 "AG08BZ.f"
	i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + jwork - 1;
#line 522 "AG08BZ.f"
	wrkopt = max(i__1,i__2);
#line 523 "AG08BZ.f"
    }

#line 525 "AG08BZ.f"
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

#line 538 "AG08BZ.f"
	numu = nu + mu;
#line 539 "AG08BZ.f"
	ipd = kabcd + nu;
#line 540 "AG08BZ.f"
	itau = jwork;
#line 541 "AG08BZ.f"
	jwork = itau + mu;

/*        CWorkspace: need   LABCD2 + 2*min(M,P); */
/*                    prefer LABCD2 + min(M,P) + min(M,P)*NB. */

#line 546 "AG08BZ.f"
	i__1 = *lzwork - jwork + 1;
#line 546 "AG08BZ.f"
	ztzrzf_(&mu, &numu, &zwork[ipd], &ldabcd, &zwork[itau], &zwork[jwork],
		 &i__1, info);
/* Computing MAX */
#line 548 "AG08BZ.f"
	i__3 = jwork;
#line 548 "AG08BZ.f"
	i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + jwork - 1;
#line 548 "AG08BZ.f"
	wrkopt = max(i__1,i__2);

/*        CWorkspace: need   LABCD2 + min(M,P) + min(L,N); */
/*                    prefer LABCD2 + min(M,P) + min(L,N)*NB. */

#line 553 "AG08BZ.f"
	i__1 = *lzwork - jwork + 1;
#line 553 "AG08BZ.f"
	zunmrz_("Right", "Conjugate transpose", &nu, &numu, &mu, &nu, &zwork[
		ipd], &ldabcd, &zwork[itau], &zwork[kabcd], &ldabcd, &zwork[
		jwork], &i__1, info, (ftnlen)5, (ftnlen)19);
/* Computing MAX */
#line 556 "AG08BZ.f"
	i__3 = jwork;
#line 556 "AG08BZ.f"
	i__1 = wrkopt, i__2 = (integer) zwork[i__3].r + jwork - 1;
#line 556 "AG08BZ.f"
	wrkopt = max(i__1,i__2);

/*        Save Af. */

#line 560 "AG08BZ.f"
	zlacpy_("Full", &nu, &nu, &zwork[kabcd + ldabcd * mu], &ldabcd, &a[
		a_offset], lda, (ftnlen)4);

/*        Compute Ef by applying the saved transformations from previous */
/*        reduction to ( 0  Er ) . */

#line 566 "AG08BZ.f"
	zlaset_("Full", &nu, &mu, &c_b1, &c_b1, &zwork[kabcd], &ldabcd, (
		ftnlen)4);
#line 568 "AG08BZ.f"
	zlacpy_("Full", &nu, &nu, &e[e_offset], lde, &zwork[kabcd + ldabcd * 
		mu], &ldabcd, (ftnlen)4);

#line 571 "AG08BZ.f"
	i__1 = *lzwork - jwork + 1;
#line 571 "AG08BZ.f"
	zunmrz_("Right", "Conjugate transpose", &nu, &numu, &mu, &nu, &zwork[
		ipd], &ldabcd, &zwork[itau], &zwork[kabcd], &ldabcd, &zwork[
		jwork], &i__1, info, (ftnlen)5, (ftnlen)19);

/*        Save Ef. */

#line 577 "AG08BZ.f"
	zlacpy_("Full", &nu, &nu, &zwork[kabcd + ldabcd * mu], &ldabcd, &e[
		e_offset], lde, (ftnlen)4);
#line 579 "AG08BZ.f"
    }

#line 581 "AG08BZ.f"
    *nfz = nu;

/*     Set right Kronecker indices (column indices). */

#line 585 "AG08BZ.f"
    i__1 = *nkror;
#line 585 "AG08BZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 586 "AG08BZ.f"
	iwork[i__] = kronr[i__];
#line 587 "AG08BZ.f"
/* L10: */
#line 587 "AG08BZ.f"
    }

#line 589 "AG08BZ.f"
    j = 0;
#line 590 "AG08BZ.f"
    i__1 = *nkror;
#line 590 "AG08BZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 591 "AG08BZ.f"
	i__2 = j + iwork[i__];
#line 591 "AG08BZ.f"
	for (ii = j + 1; ii <= i__2; ++ii) {
#line 592 "AG08BZ.f"
	    kronr[ii] = i__ - 1;
#line 593 "AG08BZ.f"
/* L20: */
#line 593 "AG08BZ.f"
	}
#line 594 "AG08BZ.f"
	j += iwork[i__];
#line 595 "AG08BZ.f"
/* L30: */
#line 595 "AG08BZ.f"
    }

#line 597 "AG08BZ.f"
    *nkror = j;

/*     Set left Kronecker indices (row indices). */

#line 601 "AG08BZ.f"
    i__1 = *nkrol;
#line 601 "AG08BZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 602 "AG08BZ.f"
	iwork[i__] = kronl[i__];
#line 603 "AG08BZ.f"
/* L40: */
#line 603 "AG08BZ.f"
    }

#line 605 "AG08BZ.f"
    j = 0;
#line 606 "AG08BZ.f"
    i__1 = *nkrol;
#line 606 "AG08BZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 607 "AG08BZ.f"
	i__2 = j + iwork[i__];
#line 607 "AG08BZ.f"
	for (ii = j + 1; ii <= i__2; ++ii) {
#line 608 "AG08BZ.f"
	    kronl[ii] = i__ - 1;
#line 609 "AG08BZ.f"
/* L50: */
#line 609 "AG08BZ.f"
	}
#line 610 "AG08BZ.f"
	j += iwork[i__];
#line 611 "AG08BZ.f"
/* L60: */
#line 611 "AG08BZ.f"
    }

#line 613 "AG08BZ.f"
    *nkrol = j;

/*     Determine the number of simple infinite blocks */
/*     as the difference between the number of infinite blocks */
/*     of order greater than one and the order of Dr. */

#line 619 "AG08BZ.f"
    *ninfe = 0;
#line 620 "AG08BZ.f"
    i__1 = *dinfz;
#line 620 "AG08BZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 621 "AG08BZ.f"
	*ninfe += infz[i__];
#line 622 "AG08BZ.f"
/* L70: */
#line 622 "AG08BZ.f"
    }
#line 623 "AG08BZ.f"
    *ninfe = nsinfe - *ninfe;
#line 624 "AG08BZ.f"
    i__1 = *ninfe;
#line 624 "AG08BZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 625 "AG08BZ.f"
	infe[i__] = 1;
#line 626 "AG08BZ.f"
/* L80: */
#line 626 "AG08BZ.f"
    }

/*     Set the structure of infinite eigenvalues. */

#line 630 "AG08BZ.f"
    i__1 = *dinfz;
#line 630 "AG08BZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 631 "AG08BZ.f"
	i__2 = *ninfe + infz[i__];
#line 631 "AG08BZ.f"
	for (ii = *ninfe + 1; ii <= i__2; ++ii) {
#line 632 "AG08BZ.f"
	    infe[ii] = i__ + 1;
#line 633 "AG08BZ.f"
/* L90: */
#line 633 "AG08BZ.f"
	}
#line 634 "AG08BZ.f"
	*ninfe += infz[i__];
#line 635 "AG08BZ.f"
/* L100: */
#line 635 "AG08BZ.f"
    }

#line 637 "AG08BZ.f"
    iwork[1] = nsinfe;
#line 638 "AG08BZ.f"
    zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
#line 639 "AG08BZ.f"
    return 0;
/* *** Last line of AG08BZ *** */
} /* ag08bz_ */

