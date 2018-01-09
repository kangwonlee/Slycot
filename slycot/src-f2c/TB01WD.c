#line 1 "TB01WD.f"
/* TB01WD.f -- translated by f2c (version 20100827).
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

#line 1 "TB01WD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = 1.;
static doublereal c_b11 = 0.;

/* Subroutine */ int tb01wd_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *u, integer *ldu, doublereal *wr, doublereal 
	*wi, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, u_dim1, 
	    u_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, sdim, ldwp;
    extern /* Subroutine */ int dgees_(char *, char *, L_fp, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen), dgemm_(char *, char *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen), dgemv_(char *, integer *, integer *, doublereal *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical bwork[1];
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
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

/*     To reduce the system state matrix A to an upper real Schur form */
/*     by using an orthogonal similarity transformation A <-- U'*A*U and */
/*     to apply the transformation to the matrices B and C: B <-- U'*B */
/*     and C <-- C*U. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, */
/*             i.e. the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs, or of columns of B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs, or of rows of C.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix U' * A * U in real Schur form. The elements */
/*             below the first subdiagonal are set to zero. */
/*             Note:  A matrix is in real Schur form if it is upper */
/*                    quasi-triangular with 1-by-1 and 2-by-2 blocks. */
/*                    2-by-2 blocks are standardized in the form */
/*                             [  a  b  ] */
/*                             [  c  a  ] */
/*                    where b*c < 0. The eigenvalues of such a block */
/*                    are a +- sqrt(bc). */

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

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,N) */
/*             The leading N-by-N part of this array contains the */
/*             orthogonal transformation matrix used to reduce A to the */
/*             real Schur form. The columns of U are the Schur vectors of */
/*             matrix A. */

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
/*             The dimension of working array DWORK.  LWORK >= 3*N. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, the QR algorithm failed to compute */
/*                   all the eigenvalues; elements i+1:N of WR and WI */
/*                   contain those eigenvalues which have converged; */
/*                   U contains the matrix which reduces A to its */
/*                   partially converged Schur form. */

/*     METHOD */

/*     Matrix A is reduced to a real Schur form using an orthogonal */
/*     similarity transformation A <- U'*A*U. Then, the transformation */
/*     is applied to the matrices B and C: B <-- U'*B and C <-- C*U. */

/*     NUMERICAL ASPECTS */
/*                                     3 */
/*     The algorithm requires about 10N  floating point operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, March 1998. */
/*     Based on the RASP routine SRSFDC. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Orthogonal transformation, real Schur form, similarity */
/*     transformation. */

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

#line 168 "TB01WD.f"
    /* Parameter adjustments */
#line 168 "TB01WD.f"
    a_dim1 = *lda;
#line 168 "TB01WD.f"
    a_offset = 1 + a_dim1;
#line 168 "TB01WD.f"
    a -= a_offset;
#line 168 "TB01WD.f"
    b_dim1 = *ldb;
#line 168 "TB01WD.f"
    b_offset = 1 + b_dim1;
#line 168 "TB01WD.f"
    b -= b_offset;
#line 168 "TB01WD.f"
    c_dim1 = *ldc;
#line 168 "TB01WD.f"
    c_offset = 1 + c_dim1;
#line 168 "TB01WD.f"
    c__ -= c_offset;
#line 168 "TB01WD.f"
    u_dim1 = *ldu;
#line 168 "TB01WD.f"
    u_offset = 1 + u_dim1;
#line 168 "TB01WD.f"
    u -= u_offset;
#line 168 "TB01WD.f"
    --wr;
#line 168 "TB01WD.f"
    --wi;
#line 168 "TB01WD.f"
    --dwork;
#line 168 "TB01WD.f"

#line 168 "TB01WD.f"
    /* Function Body */
#line 168 "TB01WD.f"
    *info = 0;

/*     Check input parameters. */

#line 172 "TB01WD.f"
    if (*n < 0) {
#line 173 "TB01WD.f"
	*info = -1;
#line 174 "TB01WD.f"
    } else if (*m < 0) {
#line 175 "TB01WD.f"
	*info = -2;
#line 176 "TB01WD.f"
    } else if (*p < 0) {
#line 177 "TB01WD.f"
	*info = -3;
#line 178 "TB01WD.f"
    } else if (*lda < max(1,*n)) {
#line 179 "TB01WD.f"
	*info = -5;
#line 180 "TB01WD.f"
    } else if (*ldb < max(1,*n)) {
#line 181 "TB01WD.f"
	*info = -7;
#line 182 "TB01WD.f"
    } else if (*ldc < max(1,*p)) {
#line 183 "TB01WD.f"
	*info = -9;
#line 184 "TB01WD.f"
    } else if (*ldu < max(1,*n)) {
#line 185 "TB01WD.f"
	*info = -11;
#line 186 "TB01WD.f"
    } else if (*ldwork < *n * 3) {
#line 187 "TB01WD.f"
	*info = -15;
#line 188 "TB01WD.f"
    }

#line 190 "TB01WD.f"
    if (*info != 0) {

/*        Error return. */

#line 194 "TB01WD.f"
	i__1 = -(*info);
#line 194 "TB01WD.f"
	xerbla_("TB01WD", &i__1, (ftnlen)6);
#line 195 "TB01WD.f"
	return 0;
#line 196 "TB01WD.f"
    }

/*     Quick return if possible. */

#line 200 "TB01WD.f"
    if (*n == 0) {
#line 200 "TB01WD.f"
	return 0;
#line 200 "TB01WD.f"
    }

/*     Reduce A to real Schur form using an orthogonal similarity */
/*     transformation A <- U'*A*U, accumulate the transformation in U */
/*     and compute the eigenvalues of A in (WR,WI). */

/*     Workspace:  need   3*N; */
/*                 prefer larger. */

#line 210 "TB01WD.f"
    dgees_("Vectors", "Not ordered", (L_fp)select_, n, &a[a_offset], lda, &
	    sdim, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], ldwork, bwork,
	     info, (ftnlen)7, (ftnlen)11);
#line 212 "TB01WD.f"
    wrkopt = dwork[1];
#line 213 "TB01WD.f"
    if (*info != 0) {
#line 213 "TB01WD.f"
	return 0;
#line 213 "TB01WD.f"
    }

/*     Apply the transformation: B <-- U'*B. */

#line 218 "TB01WD.f"
    if (*ldwork < *n * *m) {

/*        Not enough working space for using DGEMM. */

#line 222 "TB01WD.f"
	i__1 = *m;
#line 222 "TB01WD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 223 "TB01WD.f"
	    dcopy_(n, &b[i__ * b_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 224 "TB01WD.f"
	    dgemv_("Transpose", n, n, &c_b9, &u[u_offset], ldu, &dwork[1], &
		    c__1, &c_b11, &b[i__ * b_dim1 + 1], &c__1, (ftnlen)9);
#line 226 "TB01WD.f"
/* L10: */
#line 226 "TB01WD.f"
	}

#line 228 "TB01WD.f"
    } else {
#line 229 "TB01WD.f"
	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
#line 230 "TB01WD.f"
	dgemm_("Transpose", "No transpose", n, m, n, &c_b9, &u[u_offset], ldu,
		 &dwork[1], n, &c_b11, &b[b_offset], ldb, (ftnlen)9, (ftnlen)
		12);
/* Computing MAX */
#line 232 "TB01WD.f"
	d__1 = wrkopt, d__2 = (doublereal) (*n * *m);
#line 232 "TB01WD.f"
	wrkopt = max(d__1,d__2);
#line 233 "TB01WD.f"
    }

/*     Apply the transformation: C <-- C*U. */

#line 237 "TB01WD.f"
    if (*ldwork < *n * *p) {

/*        Not enough working space for using DGEMM. */

#line 241 "TB01WD.f"
	i__1 = *p;
#line 241 "TB01WD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 242 "TB01WD.f"
	    dcopy_(n, &c__[i__ + c_dim1], ldc, &dwork[1], &c__1);
#line 243 "TB01WD.f"
	    dgemv_("Transpose", n, n, &c_b9, &u[u_offset], ldu, &dwork[1], &
		    c__1, &c_b11, &c__[i__ + c_dim1], ldc, (ftnlen)9);
#line 245 "TB01WD.f"
/* L20: */
#line 245 "TB01WD.f"
	}

#line 247 "TB01WD.f"
    } else {
#line 248 "TB01WD.f"
	ldwp = max(1,*p);
#line 249 "TB01WD.f"
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[1], &ldwp, (ftnlen)
		4);
#line 250 "TB01WD.f"
	dgemm_("No transpose", "No transpose", p, n, n, &c_b9, &dwork[1], &
		ldwp, &u[u_offset], ldu, &c_b11, &c__[c_offset], ldc, (ftnlen)
		12, (ftnlen)12);
/* Computing MAX */
#line 252 "TB01WD.f"
	d__1 = wrkopt, d__2 = (doublereal) (*n * *p);
#line 252 "TB01WD.f"
	wrkopt = max(d__1,d__2);
#line 253 "TB01WD.f"
    }

#line 255 "TB01WD.f"
    dwork[1] = wrkopt;

#line 257 "TB01WD.f"
    return 0;
/* *** Last line of TB01WD *** */
} /* tb01wd_ */
