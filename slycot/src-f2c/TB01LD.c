#line 1 "TB01LD.f"
/* TB01LD.f -- translated by f2c (version 20100827).
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

#line 1 "TB01LD.f"
/* Table of constant values */

static doublereal c_b13 = 0.;
static doublereal c_b14 = 1.;
static integer c__1 = 1;

/* Subroutine */ int tb01ld_(char *dico, char *stdom, char *joba, integer *n, 
	integer *m, integer *p, doublereal *alpha, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	integer *ndim, doublereal *u, integer *ldu, doublereal *wr, 
	doublereal *wi, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen dico_len, ftnlen stdom_len, ftnlen joba_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, u_dim1, 
	    u_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, sdim, ierr, ldwp;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), mb03qd_(char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static logical ljobg;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static logical discr;
    extern /* Subroutine */ int mb03qx_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical bwork[1];
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    extern logical select_();
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

/*     To reduce the system state matrix A to an ordered upper real */
/*     Schur form by using an orthogonal similarity transformation */
/*     A <-- U'*A*U and to apply the transformation to the matrices */
/*     B and C: B <-- U'*B and C <-- C*U. */
/*     The leading block of the resulting A has eigenvalues in a */
/*     suitably defined domain of interest. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     STDOM   CHARACTER*1 */
/*             Specifies whether the domain of interest is of stability */
/*             type (left part of complex plane or inside of a circle) */
/*             or of instability type (right part of complex plane or */
/*             outside of a circle) as follows: */
/*             = 'S':  stability type domain; */
/*             = 'U':  instability type domain. */

/*     JOBA    CHARACTER*1 */
/*             Specifies the shape of the state dynamics matrix on entry */
/*             as follows: */
/*             = 'S':  A is in an upper real Schur form; */
/*             = 'G':  A is a general square dense matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation, */
/*             i.e. the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs, or of columns of B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs, or of rows of C.  P >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION. */
/*             Specifies the boundary of the domain of interest for the */
/*             eigenvalues of A. For a continuous-time system */
/*             (DICO = 'C'), ALPHA is the boundary value for the real */
/*             parts of eigenvalues, while for a discrete-time system */
/*             (DICO = 'D'), ALPHA >= 0 represents the boundary value */
/*             for the moduli of eigenvalues. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the unreduced state dynamics matrix A. */
/*             If JOBA = 'S' then A must be a matrix in real Schur form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the ordered real Schur matrix U' * A * U with the elements */
/*             below the first subdiagonal set to zero. */
/*             The leading NDIM-by-NDIM part of A has eigenvalues in the */
/*             domain of interest and the trailing (N-NDIM)-by-(N-NDIM) */
/*             part has eigenvalues outside the domain of interest. */
/*             The domain of interest for lambda(A), the eigenvalues */
/*             of A, is defined by the parameters ALPHA, DICO and STDOM */
/*             as follows: */
/*             For a continuous-time system (DICO = 'C'): */
/*               Real(lambda(A)) < ALPHA if STDOM = 'S'; */
/*               Real(lambda(A)) > ALPHA if STDOM = 'U'; */
/*             For a discrete-time system (DICO = 'D'): */
/*               Abs(lambda(A)) < ALPHA if STDOM = 'S'; */
/*               Abs(lambda(A)) > ALPHA if STDOM = 'U'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed input matrix U' * B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the transformed output matrix C * U. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     NDIM    (output) INTEGER */
/*             The number of eigenvalues of A lying inside the domain of */
/*             interest for eigenvalues. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,N) */
/*             The leading N-by-N part of this array contains the */
/*             orthogonal transformation matrix used to reduce A to the */
/*             real Schur form and/or to reorder the diagonal blocks of */
/*             real Schur form of A. The first NDIM columns of U form */
/*             an orthogonal basis for the invariant subspace of A */
/*             corresponding to the first NDIM eigenvalues. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= max(1,N). */

/*     WR, WI  (output) DOUBLE PRECISION arrays, dimension (N) */
/*             WR and WI contain the real and imaginary parts, */
/*             respectively, of the computed eigenvalues of A. The */
/*             eigenvalues will be in the same order that they appear on */
/*             the diagonal of the output real Schur form of A. Complex */
/*             conjugate pairs of eigenvalues will appear consecutively */
/*             with the eigenvalue having the positive imaginary part */
/*             first. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of working array DWORK. */
/*             LDWORK >= MAX(1,N)   if JOBA = 'S'; */
/*             LDWORK >= MAX(1,3*N) if JOBA = 'G'. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0: successful exit; */
/*             < 0: if INFO = -i, the i-th argument had an illegal */
/*                  value; */
/*             = 1: the QR algorithm failed to compute all the */
/*                  eigenvalues of A; */
/*             = 2: a failure occured during the ordering of the real */
/*                  Schur form of A. */

/*     METHOD */

/*     Matrix A is reduced to an ordered upper real Schur form using an */
/*     orthogonal similarity transformation A <-- U'*A*U. This */
/*     transformation is determined so that the leading block of the */
/*     resulting A has eigenvalues in a suitably defined domain of */
/*     interest. Then, the transformation is applied to the matrices B */
/*     and C: B <-- U'*B and C <-- C*U. */

/*     NUMERICAL ASPECTS */
/*                                     3 */
/*     The algorithm requires about 14N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routine SRSFOD. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2001. */

/*     KEYWORDS */

/*     Invariant subspace, orthogonal transformation, real Schur form, */
/*     similarity transformation. */

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

#line 220 "TB01LD.f"
    /* Parameter adjustments */
#line 220 "TB01LD.f"
    a_dim1 = *lda;
#line 220 "TB01LD.f"
    a_offset = 1 + a_dim1;
#line 220 "TB01LD.f"
    a -= a_offset;
#line 220 "TB01LD.f"
    b_dim1 = *ldb;
#line 220 "TB01LD.f"
    b_offset = 1 + b_dim1;
#line 220 "TB01LD.f"
    b -= b_offset;
#line 220 "TB01LD.f"
    c_dim1 = *ldc;
#line 220 "TB01LD.f"
    c_offset = 1 + c_dim1;
#line 220 "TB01LD.f"
    c__ -= c_offset;
#line 220 "TB01LD.f"
    u_dim1 = *ldu;
#line 220 "TB01LD.f"
    u_offset = 1 + u_dim1;
#line 220 "TB01LD.f"
    u -= u_offset;
#line 220 "TB01LD.f"
    --wr;
#line 220 "TB01LD.f"
    --wi;
#line 220 "TB01LD.f"
    --dwork;
#line 220 "TB01LD.f"

#line 220 "TB01LD.f"
    /* Function Body */
#line 220 "TB01LD.f"
    *info = 0;
#line 221 "TB01LD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 222 "TB01LD.f"
    ljobg = lsame_(joba, "G", (ftnlen)1, (ftnlen)1);

/*     Check input scalar arguments. */

#line 226 "TB01LD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 227 "TB01LD.f"
	*info = -1;
#line 228 "TB01LD.f"
    } else if (! (lsame_(stdom, "S", (ftnlen)1, (ftnlen)1) || lsame_(stdom, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 230 "TB01LD.f"
	*info = -2;
#line 231 "TB01LD.f"
    } else if (! (lsame_(joba, "S", (ftnlen)1, (ftnlen)1) || ljobg)) {
#line 232 "TB01LD.f"
	*info = -3;
#line 233 "TB01LD.f"
    } else if (*n < 0) {
#line 234 "TB01LD.f"
	*info = -4;
#line 235 "TB01LD.f"
    } else if (*m < 0) {
#line 236 "TB01LD.f"
	*info = -5;
#line 237 "TB01LD.f"
    } else if (*p < 0) {
#line 238 "TB01LD.f"
	*info = -6;
#line 239 "TB01LD.f"
    } else if (discr && *alpha < 0.) {
#line 240 "TB01LD.f"
	*info = -7;
#line 241 "TB01LD.f"
    } else if (*lda < max(1,*n)) {
#line 242 "TB01LD.f"
	*info = -9;
#line 243 "TB01LD.f"
    } else if (*ldb < max(1,*n)) {
#line 244 "TB01LD.f"
	*info = -11;
#line 245 "TB01LD.f"
    } else if (*ldc < max(1,*p)) {
#line 246 "TB01LD.f"
	*info = -13;
#line 247 "TB01LD.f"
    } else if (*ldu < max(1,*n)) {
#line 248 "TB01LD.f"
	*info = -16;
#line 249 "TB01LD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 249 "TB01LD.f"
	i__1 = 1, i__2 = *n * 3;
#line 249 "TB01LD.f"
	if (*ldwork < max(1,*n) || *ldwork < max(i__1,i__2) && ljobg) {
#line 251 "TB01LD.f"
	    *info = -20;
#line 252 "TB01LD.f"
	}
#line 252 "TB01LD.f"
    }

#line 254 "TB01LD.f"
    if (*info != 0) {

/*        Error return. */

#line 258 "TB01LD.f"
	i__1 = -(*info);
#line 258 "TB01LD.f"
	xerbla_("TB01LD", &i__1, (ftnlen)6);
#line 259 "TB01LD.f"
	return 0;
#line 260 "TB01LD.f"
    }

/*     Quick return if possible. */

#line 264 "TB01LD.f"
    *ndim = 0;
#line 265 "TB01LD.f"
    if (*n == 0) {
#line 265 "TB01LD.f"
	return 0;
#line 265 "TB01LD.f"
    }

#line 268 "TB01LD.f"
    if (lsame_(joba, "G", (ftnlen)1, (ftnlen)1)) {

/*        Reduce A to real Schur form using an orthogonal similarity */
/*        transformation A <- U'*A*U, accumulate the transformation in U */
/*        and compute the eigenvalues of A in (WR,WI). */

/*        Workspace:  need   3*N; */
/*                    prefer larger. */

#line 277 "TB01LD.f"
	dgees_("Vectors", "Not ordered", (L_fp)select_, n, &a[a_offset], lda, 
		&sdim, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], ldwork, 
		bwork, info, (ftnlen)7, (ftnlen)11);
#line 279 "TB01LD.f"
	wrkopt = dwork[1];
#line 280 "TB01LD.f"
	if (*info != 0) {
#line 281 "TB01LD.f"
	    *info = 1;
#line 282 "TB01LD.f"
	    return 0;
#line 283 "TB01LD.f"
	}
#line 284 "TB01LD.f"
    } else {

/*        Initialize U with an identity matrix. */

#line 288 "TB01LD.f"
	dlaset_("Full", n, n, &c_b13, &c_b14, &u[u_offset], ldu, (ftnlen)4);
#line 289 "TB01LD.f"
	wrkopt = 0.;
#line 290 "TB01LD.f"
    }

/*     Separate the spectrum of A. The leading NDIM-by-NDIM submatrix of */
/*     A corresponds to the eigenvalues of interest. */
/*     Workspace:  need   N. */

#line 296 "TB01LD.f"
    mb03qd_(dico, stdom, "Update", n, &c__1, n, alpha, &a[a_offset], lda, &u[
	    u_offset], ldu, ndim, &dwork[1], info, (ftnlen)1, (ftnlen)1, (
	    ftnlen)6);
#line 298 "TB01LD.f"
    if (*info != 0) {
#line 298 "TB01LD.f"
	return 0;
#line 298 "TB01LD.f"
    }

/*     Compute the eigenvalues. */

#line 303 "TB01LD.f"
    mb03qx_(n, &a[a_offset], lda, &wr[1], &wi[1], &ierr);

/*     Apply the transformation: B <-- U'*B. */

#line 307 "TB01LD.f"
    if (*ldwork < *n * *m) {

/*        Not enough working space for using DGEMM. */

#line 311 "TB01LD.f"
	i__1 = *m;
#line 311 "TB01LD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 312 "TB01LD.f"
	    dcopy_(n, &b[i__ * b_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 313 "TB01LD.f"
	    dgemv_("Transpose", n, n, &c_b14, &u[u_offset], ldu, &dwork[1], &
		    c__1, &c_b13, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)9);
#line 315 "TB01LD.f"
/* L10: */
#line 315 "TB01LD.f"
	}

#line 317 "TB01LD.f"
    } else {
#line 318 "TB01LD.f"
	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
#line 319 "TB01LD.f"
	dgemm_("Transpose", "No transpose", n, m, n, &c_b14, &u[u_offset], 
		ldu, &dwork[1], n, &c_b13, &b[b_offset], ldb, (ftnlen)9, (
		ftnlen)12);
/* Computing MAX */
#line 321 "TB01LD.f"
	d__1 = wrkopt, d__2 = (doublereal) (*n * *m);
#line 321 "TB01LD.f"
	wrkopt = max(d__1,d__2);
#line 322 "TB01LD.f"
    }

/*     Apply the transformation: C <-- C*U. */

#line 326 "TB01LD.f"
    if (*ldwork < *n * *p) {

/*        Not enough working space for using DGEMM. */

#line 330 "TB01LD.f"
	i__1 = *p;
#line 330 "TB01LD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 331 "TB01LD.f"
	    dcopy_(n, &c__[i__ + c_dim1], ldc, &dwork[1], &c__1);
#line 332 "TB01LD.f"
	    dgemv_("Transpose", n, n, &c_b14, &u[u_offset], ldu, &dwork[1], &
		    c__1, &c_b13, &c__[i__ + c_dim1], ldc, (ftnlen)9);
#line 334 "TB01LD.f"
/* L20: */
#line 334 "TB01LD.f"
	}

#line 336 "TB01LD.f"
    } else {
#line 337 "TB01LD.f"
	ldwp = max(1,*p);
#line 338 "TB01LD.f"
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[1], &ldwp, (ftnlen)
		4);
#line 339 "TB01LD.f"
	dgemm_("No transpose", "No transpose", p, n, n, &c_b14, &dwork[1], &
		ldwp, &u[u_offset], ldu, &c_b13, &c__[c_offset], ldc, (ftnlen)
		12, (ftnlen)12);
/* Computing MAX */
#line 341 "TB01LD.f"
	d__1 = wrkopt, d__2 = (doublereal) (*n * *p);
#line 341 "TB01LD.f"
	wrkopt = max(d__1,d__2);
#line 342 "TB01LD.f"
    }

#line 344 "TB01LD.f"
    dwork[1] = wrkopt;

#line 346 "TB01LD.f"
    return 0;
/* *** Last line of TB01LD *** */
} /* tb01ld_ */

