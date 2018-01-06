#line 1 "MB02RD.f"
/* MB02RD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02RD.f"
/* Table of constant values */

static doublereal c_b11 = 1.;

/* Subroutine */ int mb02rd_(char *trans, integer *n, integer *nrhs, 
	doublereal *h__, integer *ldh, integer *ipiv, doublereal *b, integer *
	ldb, integer *info, ftnlen trans_len)
{
    /* System generated locals */
    integer b_dim1, b_offset, h_dim1, h_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j, jp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), dtrsm_(char *, 
	    char *, char *, char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    static logical notran;


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

/*     To solve a system of linear equations */
/*        H * X = B  or  H' * X = B */
/*     with an upper Hessenberg N-by-N matrix H using the LU */
/*     factorization computed by MB02SD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANS   CHARACTER*1 */
/*             Specifies the form of the system of equations: */
/*             = 'N':  H * X = B  (No transpose) */
/*             = 'T':  H'* X = B  (Transpose) */
/*             = 'C':  H'* X = B  (Conjugate transpose = Transpose) */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix H.  N >= 0. */

/*     NRHS    (input) INTEGER */
/*             The number of right hand sides, i.e., the number of */
/*             columns of the matrix B.  NRHS >= 0. */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,N) */
/*             The factors L and U from the factorization H = P*L*U */
/*             as computed by MB02SD. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

/*     IPIV    (input) INTEGER array, dimension (N) */
/*             The pivot indices from MB02SD; for 1<=i<=N, row i of the */
/*             matrix was interchanged with row IPIV(i). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,NRHS) */
/*             On entry, the right hand side matrix B. */
/*             On exit, the solution matrix X. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     Error Indicator */

/*     INFO    (output) INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine uses the factorization */
/*        H = P * L * U */
/*     where P is a permutation matrix, L is lower triangular with unit */
/*     diagonal elements (and one nonzero subdiagonal), and U is upper */
/*     triangular. */

/*     REFERENCES */

/*     - */

/*     NUMERICAL ASPECTS */
/*                                2 */
/*     The algorithm requires 0( N x NRHS ) operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Hessenberg form, matrix algebra. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 127 "MB02RD.f"
    /* Parameter adjustments */
#line 127 "MB02RD.f"
    h_dim1 = *ldh;
#line 127 "MB02RD.f"
    h_offset = 1 + h_dim1;
#line 127 "MB02RD.f"
    h__ -= h_offset;
#line 127 "MB02RD.f"
    --ipiv;
#line 127 "MB02RD.f"
    b_dim1 = *ldb;
#line 127 "MB02RD.f"
    b_offset = 1 + b_dim1;
#line 127 "MB02RD.f"
    b -= b_offset;
#line 127 "MB02RD.f"

#line 127 "MB02RD.f"
    /* Function Body */
#line 127 "MB02RD.f"
    *info = 0;
#line 128 "MB02RD.f"
    notran = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 129 "MB02RD.f"
    if (! notran && ! lsame_(trans, "T", (ftnlen)1, (ftnlen)1) && ! lsame_(
	    trans, "C", (ftnlen)1, (ftnlen)1)) {
#line 131 "MB02RD.f"
	*info = -1;
#line 132 "MB02RD.f"
    } else if (*n < 0) {
#line 133 "MB02RD.f"
	*info = -2;
#line 134 "MB02RD.f"
    } else if (*nrhs < 0) {
#line 135 "MB02RD.f"
	*info = -3;
#line 136 "MB02RD.f"
    } else if (*ldh < max(1,*n)) {
#line 137 "MB02RD.f"
	*info = -5;
#line 138 "MB02RD.f"
    } else if (*ldb < max(1,*n)) {
#line 139 "MB02RD.f"
	*info = -8;
#line 140 "MB02RD.f"
    }
#line 141 "MB02RD.f"
    if (*info != 0) {
#line 142 "MB02RD.f"
	i__1 = -(*info);
#line 142 "MB02RD.f"
	xerbla_("MB02RD", &i__1, (ftnlen)6);
#line 143 "MB02RD.f"
	return 0;
#line 144 "MB02RD.f"
    }

/*     Quick return if possible. */

#line 148 "MB02RD.f"
    if (*n == 0 || *nrhs == 0) {
#line 148 "MB02RD.f"
	return 0;
#line 148 "MB02RD.f"
    }

#line 151 "MB02RD.f"
    if (notran) {

/*        Solve H * X = B. */

/*        Solve L * X = B, overwriting B with X. */

/*        L is represented as a product of permutations and unit lower */
/*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1), */
/*        where each transformation L(i) is a rank-one modification of */
/*        the identity matrix. */

#line 162 "MB02RD.f"
	i__1 = *n - 1;
#line 162 "MB02RD.f"
	for (j = 1; j <= i__1; ++j) {
#line 163 "MB02RD.f"
	    jp = ipiv[j];
#line 164 "MB02RD.f"
	    if (jp != j) {
#line 164 "MB02RD.f"
		dswap_(nrhs, &b[jp + b_dim1], ldb, &b[j + b_dim1], ldb);
#line 164 "MB02RD.f"
	    }
#line 166 "MB02RD.f"
	    d__1 = -h__[j + 1 + j * h_dim1];
#line 166 "MB02RD.f"
	    daxpy_(nrhs, &d__1, &b[j + b_dim1], ldb, &b[j + 1 + b_dim1], ldb);
#line 168 "MB02RD.f"
/* L10: */
#line 168 "MB02RD.f"
	}

/*        Solve U * X = B, overwriting B with X. */

#line 172 "MB02RD.f"
	dtrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b11, &
		h__[h_offset], ldh, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)12, (ftnlen)8);

#line 175 "MB02RD.f"
    } else {

/*        Solve H' * X = B. */

/*        Solve U' * X = B, overwriting B with X. */

#line 181 "MB02RD.f"
	dtrsm_("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b11, &
		h__[h_offset], ldh, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
		ftnlen)9, (ftnlen)8);

/*        Solve L' * X = B, overwriting B with X. */

#line 186 "MB02RD.f"
	for (j = *n - 1; j >= 1; --j) {
#line 187 "MB02RD.f"
	    d__1 = -h__[j + 1 + j * h_dim1];
#line 187 "MB02RD.f"
	    daxpy_(nrhs, &d__1, &b[j + 1 + b_dim1], ldb, &b[j + b_dim1], ldb);
#line 189 "MB02RD.f"
	    jp = ipiv[j];
#line 190 "MB02RD.f"
	    if (jp != j) {
#line 190 "MB02RD.f"
		dswap_(nrhs, &b[jp + b_dim1], ldb, &b[j + b_dim1], ldb);
#line 190 "MB02RD.f"
	    }
#line 192 "MB02RD.f"
/* L20: */
#line 192 "MB02RD.f"
	}
#line 193 "MB02RD.f"
    }

#line 195 "MB02RD.f"
    return 0;
/* *** Last line of MB02RD *** */
} /* mb02rd_ */

