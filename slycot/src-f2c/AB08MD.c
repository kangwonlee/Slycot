#line 1 "AB08MD.f"
/* AB08MD.f -- translated by f2c (version 20100827).
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

#line 1 "AB08MD.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c_n1 = -1;

/* Subroutine */ int ab08md_(char *equil, integer *n, integer *m, integer *p, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, integer *rank, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, nm, np, ro, kw, mu, nu;
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer sigma;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ab08nx_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *);
    static integer ninfz, nkrol;
    static doublereal toler;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical lequil;
    static doublereal thresh, svlmax;
    static logical lquery;
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

/*     To compute the normal rank of the transfer-function matrix of a */
/*     state-space model (A,B,C,D). */

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

/*     RANK    (output) INTEGER */
/*             The normal rank of the transfer-function matrix. */

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

/*     IWORK   INTEGER array, dimension (2*N+MAX(M,P)+1) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= (N+P)*(N+M) + */
/*                       MAX( MIN(P,M) + MAX(3*M-1,N), 1, */
/*                            MIN(P,N) + MAX(3*P-1,N+P,N+M) ) */
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

/*     The routine reduces the (N+P)-by-(M+N) compound matrix (B  A) */
/*                                                            (D  C) */

/*     to one with the same invariant zeros and with D of full row rank. */
/*     The normal rank of the transfer-function matrix is the rank of D. */

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

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     Dec. 2003, Jan. 2009, Mar. 2009, Apr. 2009. */

/*     KEYWORDS */

/*     Multivariable system, orthogonal transformation, */
/*     structural invariant. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 186 "AB08MD.f"
    /* Parameter adjustments */
#line 186 "AB08MD.f"
    a_dim1 = *lda;
#line 186 "AB08MD.f"
    a_offset = 1 + a_dim1;
#line 186 "AB08MD.f"
    a -= a_offset;
#line 186 "AB08MD.f"
    b_dim1 = *ldb;
#line 186 "AB08MD.f"
    b_offset = 1 + b_dim1;
#line 186 "AB08MD.f"
    b -= b_offset;
#line 186 "AB08MD.f"
    c_dim1 = *ldc;
#line 186 "AB08MD.f"
    c_offset = 1 + c_dim1;
#line 186 "AB08MD.f"
    c__ -= c_offset;
#line 186 "AB08MD.f"
    d_dim1 = *ldd;
#line 186 "AB08MD.f"
    d_offset = 1 + d_dim1;
#line 186 "AB08MD.f"
    d__ -= d_offset;
#line 186 "AB08MD.f"
    --iwork;
#line 186 "AB08MD.f"
    --dwork;
#line 186 "AB08MD.f"

#line 186 "AB08MD.f"
    /* Function Body */
#line 186 "AB08MD.f"
    np = *n + *p;
#line 187 "AB08MD.f"
    nm = *n + *m;
#line 188 "AB08MD.f"
    *info = 0;
#line 189 "AB08MD.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 190 "AB08MD.f"
    lquery = *ldwork == -1;
#line 191 "AB08MD.f"
    wrkopt = np * nm;

/*     Test the input scalar arguments. */

#line 195 "AB08MD.f"
    if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
#line 196 "AB08MD.f"
	*info = -1;
#line 197 "AB08MD.f"
    } else if (*n < 0) {
#line 198 "AB08MD.f"
	*info = -2;
#line 199 "AB08MD.f"
    } else if (*m < 0) {
#line 200 "AB08MD.f"
	*info = -3;
#line 201 "AB08MD.f"
    } else if (*p < 0) {
#line 202 "AB08MD.f"
	*info = -4;
#line 203 "AB08MD.f"
    } else if (*lda < max(1,*n)) {
#line 204 "AB08MD.f"
	*info = -6;
#line 205 "AB08MD.f"
    } else if (*ldb < max(1,*n)) {
#line 206 "AB08MD.f"
	*info = -8;
#line 207 "AB08MD.f"
    } else if (*ldc < max(1,*p)) {
#line 208 "AB08MD.f"
	*info = -10;
#line 209 "AB08MD.f"
    } else if (*ldd < max(1,*p)) {
#line 210 "AB08MD.f"
	*info = -12;
#line 211 "AB08MD.f"
    } else {
/* Computing MAX */
/* Computing MAX */
#line 212 "AB08MD.f"
	i__3 = *m * 3 - 1;
/* Computing MAX */
#line 212 "AB08MD.f"
	i__4 = *p * 3 - 1, i__4 = max(i__4,np);
#line 212 "AB08MD.f"
	i__1 = min(*p,*m) + max(i__3,*n), i__1 = max(i__1,1), i__2 = min(*p,*
		n) + max(i__4,nm);
#line 212 "AB08MD.f"
	kw = wrkopt + max(i__1,i__2);
#line 214 "AB08MD.f"
	if (lquery) {
#line 215 "AB08MD.f"
	    svlmax = 0.;
#line 216 "AB08MD.f"
	    ninfz = 0;
#line 217 "AB08MD.f"
	    i__1 = max(1,np);
#line 217 "AB08MD.f"
	    ab08nx_(n, m, p, p, &c__0, &svlmax, &dwork[1], &i__1, &ninfz, &
		    iwork[1], &iwork[1], &mu, &nu, &nkrol, tol, &iwork[1], &
		    dwork[1], &c_n1, info);
/* Computing MAX */
#line 220 "AB08MD.f"
	    i__1 = kw, i__2 = wrkopt + (integer) dwork[1];
#line 220 "AB08MD.f"
	    wrkopt = max(i__1,i__2);
#line 221 "AB08MD.f"
	} else if (*ldwork < kw) {
#line 222 "AB08MD.f"
	    *info = -17;
#line 223 "AB08MD.f"
	}
#line 224 "AB08MD.f"
    }

#line 226 "AB08MD.f"
    if (*info != 0) {

/*        Error return. */

#line 230 "AB08MD.f"
	i__1 = -(*info);
#line 230 "AB08MD.f"
	xerbla_("AB08MD", &i__1, (ftnlen)6);
#line 231 "AB08MD.f"
	return 0;
#line 232 "AB08MD.f"
    } else if (lquery) {
#line 233 "AB08MD.f"
	dwork[1] = (doublereal) wrkopt;
#line 234 "AB08MD.f"
	return 0;
#line 235 "AB08MD.f"
    }

/*     Quick return if possible. */

#line 239 "AB08MD.f"
    if (min(*m,*p) == 0) {
#line 240 "AB08MD.f"
	*rank = 0;
#line 241 "AB08MD.f"
	dwork[1] = 1.;
#line 242 "AB08MD.f"
	return 0;
#line 243 "AB08MD.f"
    }

#line 245 "AB08MD.f"
    i__1 = (*n << 1) + 1;
#line 245 "AB08MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 246 "AB08MD.f"
	iwork[i__] = 0;
#line 247 "AB08MD.f"
/* L10: */
#line 247 "AB08MD.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance.) */

/*     Construct the compound matrix  ( B  A ), dimension (N+P)-by-(M+N). */
/*                                    ( D  C ) */
/*     Workspace: need   (N+P)*(N+M). */

#line 257 "AB08MD.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], &np, (ftnlen)4);
#line 258 "AB08MD.f"
    dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[*n + 1], &np, (ftnlen)4)
	    ;
#line 259 "AB08MD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[np * *m + 1], &np, (
	    ftnlen)4);
#line 260 "AB08MD.f"
    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[np * *m + *n + 1], &np, 
	    (ftnlen)4);

/*     If required, balance the compound matrix (default MAXRED). */
/*     Workspace: need   N. */

#line 265 "AB08MD.f"
    kw = wrkopt + 1;
#line 266 "AB08MD.f"
    if (lequil) {
#line 267 "AB08MD.f"
	maxred = 0.;
#line 268 "AB08MD.f"
	tb01id_("A", n, m, p, &maxred, &dwork[np * *m + 1], &np, &dwork[1], &
		np, &dwork[np * *m + *n + 1], &np, &dwork[kw], info, (ftnlen)
		1);
#line 270 "AB08MD.f"
	wrkopt += *n;
#line 271 "AB08MD.f"
    }

/*     If required, set tolerance. */

#line 275 "AB08MD.f"
    thresh = sqrt((doublereal) (np * nm)) * dlamch_("Precision", (ftnlen)9);
#line 276 "AB08MD.f"
    toler = *tol;
#line 277 "AB08MD.f"
    if (toler < thresh) {
#line 277 "AB08MD.f"
	toler = thresh;
#line 277 "AB08MD.f"
    }
#line 278 "AB08MD.f"
    svlmax = dlange_("Frobenius", &np, &nm, &dwork[1], &np, &dwork[kw], (
	    ftnlen)9);

/*     Reduce this system to one with the same invariant zeros and with */
/*     D full row rank MU (the normal rank of the original system). */
/*     Real workspace: need   (N+P)*(N+M) + */
/*                            MAX( 1, MIN(P,M) + MAX(3*M-1,N), */
/*                                    MIN(P,N) + MAX(3*P-1,N+P,N+M) ); */
/*                     prefer larger. */
/*     Integer workspace: 2*N+MAX(M,P)+1. */

#line 288 "AB08MD.f"
    ro = *p;
#line 289 "AB08MD.f"
    sigma = 0;
#line 290 "AB08MD.f"
    ninfz = 0;
#line 291 "AB08MD.f"
    i__1 = *ldwork - kw + 1;
#line 291 "AB08MD.f"
    ab08nx_(n, m, p, &ro, &sigma, &svlmax, &dwork[1], &np, &ninfz, &iwork[1], 
	    &iwork[*n + 1], &mu, &nu, &nkrol, &toler, &iwork[(*n << 1) + 2], &
	    dwork[kw], &i__1, info);
#line 294 "AB08MD.f"
    *rank = mu;

/* Computing MAX */
#line 296 "AB08MD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 296 "AB08MD.f"
    dwork[1] = (doublereal) max(i__1,i__2);
#line 297 "AB08MD.f"
    return 0;
/* *** Last line of AB08MD *** */
} /* ab08md_ */
