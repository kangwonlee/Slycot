#line 1 "MB02VD.f"
/* MB02VD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02VD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b12 = 1.;
static integer c_n1 = -1;

/* Subroutine */ int mb02vd_(char *trans, integer *m, integer *n, doublereal *
	a, integer *lda, integer *ipiv, doublereal *b, integer *ldb, integer *
	info, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;

    /* Local variables */
    static logical tran;
    extern /* Subroutine */ int ma02gd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dgetrf_(
	    integer *, integer *, doublereal *, integer *, integer *, integer 
	    *), xerbla_(char *, integer *, ftnlen);


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

/*     To compute the solution to a real system of linear equations */
/*        X * op(A) = B, */
/*     where op(A) is either A or its transpose, A is an N-by-N matrix, */
/*     and X and B are M-by-N matrices. */
/*     The LU decomposition with partial pivoting and row interchanges, */
/*     A = P * L * U, is used, where P is a permutation matrix, L is unit */
/*     lower triangular, and U is upper triangular. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of op(A) to be used as follows: */
/*             = 'N':  op(A) = A; */
/*             = 'T':  op(A) = A'; */
/*             = 'C':  op(A) = A'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix B.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix B, and the order of */
/*             the matrix A.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the coefficient matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the factors L and U from the factorization A = P*L*U; */
/*             the unit diagonal elements of L are not stored. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     IPIV    (output) INTEGER array, dimension (N) */
/*             The pivot indices that define the permutation matrix P; */
/*             row i of the matrix was interchanged with row IPIV(i). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the right hand side matrix B. */
/*             On exit, if INFO = 0, the leading M-by-N part of this */
/*             array contains the solution matrix X. */

/*     LDB     (input) INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,M). */

/*     INFO    (output) INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, U(i,i) is exactly zero.  The */
/*                   factorization has been completed, but the factor U */
/*                   is exactly singular, so the solution could not be */
/*                   computed. */

/*     METHOD */

/*     The LU decomposition with partial pivoting and row interchanges is */
/*     used to factor A as */
/*        A = P * L * U, */
/*     where P is a permutation matrix, L is unit lower triangular, and */
/*     U is upper triangular.  The factored form of A is then used to */
/*     solve the system of equations X * A = B or X * A' = B. */

/*     FURTHER COMMENTS */

/*     This routine enables to solve the system X * A = B or X * A' = B */
/*     as easily and efficiently as possible; it is similar to the LAPACK */
/*     Library routine DGESV, which solves A * X = B. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, linear algebra. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the scalar input parameters. */

#line 139 "MB02VD.f"
    /* Parameter adjustments */
#line 139 "MB02VD.f"
    a_dim1 = *lda;
#line 139 "MB02VD.f"
    a_offset = 1 + a_dim1;
#line 139 "MB02VD.f"
    a -= a_offset;
#line 139 "MB02VD.f"
    --ipiv;
#line 139 "MB02VD.f"
    b_dim1 = *ldb;
#line 139 "MB02VD.f"
    b_offset = 1 + b_dim1;
#line 139 "MB02VD.f"
    b -= b_offset;
#line 139 "MB02VD.f"

#line 139 "MB02VD.f"
    /* Function Body */
#line 139 "MB02VD.f"
    *info = 0;
#line 140 "MB02VD.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);

#line 142 "MB02VD.f"
    if (! tran && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 143 "MB02VD.f"
	*info = -1;
#line 144 "MB02VD.f"
    } else if (*m < 0) {
#line 145 "MB02VD.f"
	*info = -2;
#line 146 "MB02VD.f"
    } else if (*n < 0) {
#line 147 "MB02VD.f"
	*info = -3;
#line 148 "MB02VD.f"
    } else if (*lda < max(1,*n)) {
#line 149 "MB02VD.f"
	*info = -5;
#line 150 "MB02VD.f"
    } else if (*ldb < max(1,*m)) {
#line 151 "MB02VD.f"
	*info = -8;
#line 152 "MB02VD.f"
    }

#line 154 "MB02VD.f"
    if (*info != 0) {
#line 155 "MB02VD.f"
	i__1 = -(*info);
#line 155 "MB02VD.f"
	xerbla_("MB02VD", &i__1, (ftnlen)6);
#line 156 "MB02VD.f"
	return 0;
#line 157 "MB02VD.f"
    }

/*     Compute the LU factorization of A. */

#line 161 "MB02VD.f"
    dgetrf_(n, n, &a[a_offset], lda, &ipiv[1], info);

#line 163 "MB02VD.f"
    if (*info == 0) {
#line 164 "MB02VD.f"
	if (tran) {

/*           Compute X = B * A**(-T). */

#line 168 "MB02VD.f"
	    ma02gd_(m, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c__1);
#line 169 "MB02VD.f"
	    dtrsm_("Right", "Lower", "Transpose", "Unit", m, n, &c_b12, &a[
		    a_offset], lda, &b[b_offset], ldb, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)4);
#line 171 "MB02VD.f"
	    dtrsm_("Right", "Upper", "Transpose", "NonUnit", m, n, &c_b12, &a[
		    a_offset], lda, &b[b_offset], ldb, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)7);
#line 173 "MB02VD.f"
	} else {

/*           Compute X = B * A**(-1). */

#line 177 "MB02VD.f"
	    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", m, n, &c_b12, &
		    a[a_offset], lda, &b[b_offset], ldb, (ftnlen)5, (ftnlen)5,
		     (ftnlen)11, (ftnlen)7);
#line 179 "MB02VD.f"
	    dtrsm_("Right", "Lower", "NoTranspose", "Unit", m, n, &c_b12, &a[
		    a_offset], lda, &b[b_offset], ldb, (ftnlen)5, (ftnlen)5, (
		    ftnlen)11, (ftnlen)4);
#line 181 "MB02VD.f"
	    ma02gd_(m, &b[b_offset], ldb, &c__1, n, &ipiv[1], &c_n1);
#line 182 "MB02VD.f"
	}
#line 183 "MB02VD.f"
    }
#line 184 "MB02VD.f"
    return 0;

/* *** Last line of MB02VD *** */
} /* mb02vd_ */

