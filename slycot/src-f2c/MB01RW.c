#line 1 "MB01RW.f"
/* MB01RW.f -- translated by f2c (version 20100827).
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

#line 1 "MB01RW.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = 1.;
static doublereal c_b13 = 0.;

/* Subroutine */ int mb01rw_(char *uplo, char *trans, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *z__, integer *ldz, 
	doublereal *dwork, integer *info, ftnlen uplo_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical nottra;


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

/*     To compute the transformation of the symmetric matrix A by the */
/*     matrix Z in the form */

/*        A := op(Z)*A*op(Z)', */

/*     where op(Z) is either Z or its transpose, Z'. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Specifies whether the upper or lower triangle of A */
/*             is stored: */
/*             = 'U':  Upper triangle of A is stored; */
/*             = 'L':  Lower triangle of A is stored. */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether op(Z) is Z or its transpose Z': */
/*             = 'N':  op(Z) = Z; */
/*             = 'T':  op(Z) = Z'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the resulting symmetric matrix op(Z)*A*op(Z)' */
/*             and the number of rows of the matrix Z, if TRANS = 'N', */
/*             or the number of columns of the matrix Z, if TRANS = 'T'. */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The order of the symmetric matrix A and the number of */
/*             columns of the matrix Z, if TRANS = 'N', or the number of */
/*             rows of the matrix Z, if TRANS = 'T'.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDA,MAX(M,N)) */
/*             On entry, the leading N-by-N upper or lower triangular */
/*             part of this array must contain the upper (UPLO = 'U') */
/*             or lower (UPLO = 'L') triangular part of the symmetric */
/*             matrix A. */
/*             On exit, the leading M-by-M upper or lower triangular */
/*             part of this array contains the upper (UPLO = 'U') or */
/*             lower (UPLO = 'L') triangular part of the symmetric */
/*             matrix op(Z)*A*op(Z)'. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,M,N). */

/*     Z       (input) DOUBLE PRECISION array, dimension (LDQ,K) */
/*             where K = N if TRANS = 'N' and K = M if TRANS = 'T'. */
/*             The leading M-by-N part, if TRANS = 'N', or N-by-M part, */
/*             if TRANS = 'T', of this array contains the matrix Z. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z. */
/*             LDZ >= MAX(1,M) if TRANS = 'N' and */
/*             LDZ >= MAX(1,N) if TRANS = 'T'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     FURTHER COMMENTS */

/*     This is a simpler, BLAS 2 version for MB01RD. */

/*     CONTRIBUTOR */

/*     A. Varga, DLR, Feb. 1995. */

/*     REVISIONS */

/*     April 1998 (T. Penzl). */
/*     Sep. 1998 (V. Sima). */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements */

#line 130 "MB01RW.f"
    /* Parameter adjustments */
#line 130 "MB01RW.f"
    a_dim1 = *lda;
#line 130 "MB01RW.f"
    a_offset = 1 + a_dim1;
#line 130 "MB01RW.f"
    a -= a_offset;
#line 130 "MB01RW.f"
    z_dim1 = *ldz;
#line 130 "MB01RW.f"
    z_offset = 1 + z_dim1;
#line 130 "MB01RW.f"
    z__ -= z_offset;
#line 130 "MB01RW.f"
    --dwork;
#line 130 "MB01RW.f"

#line 130 "MB01RW.f"
    /* Function Body */
#line 130 "MB01RW.f"
    nottra = lsame_(trans, "N", (ftnlen)1, (ftnlen)1);
#line 131 "MB01RW.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 133 "MB01RW.f"
    *info = 0;
#line 134 "MB01RW.f"
    if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 135 "MB01RW.f"
	*info = -1;
#line 136 "MB01RW.f"
    } else if (! (nottra || lsame_(trans, "T", (ftnlen)1, (ftnlen)1))) {
#line 137 "MB01RW.f"
	*info = -2;
#line 138 "MB01RW.f"
    } else if (*m < 0) {
#line 139 "MB01RW.f"
	*info = -3;
#line 140 "MB01RW.f"
    } else if (*n < 0) {
#line 141 "MB01RW.f"
	*info = -4;
#line 142 "MB01RW.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 142 "MB01RW.f"
	i__1 = max(1,*m);
#line 142 "MB01RW.f"
	if (*lda < max(i__1,*n)) {
#line 143 "MB01RW.f"
	    *info = -6;
#line 144 "MB01RW.f"
	} else if (nottra && *ldz < max(1,*m) || ! nottra && *ldz < max(1,*n))
		 {
#line 146 "MB01RW.f"
	    *info = -8;
#line 147 "MB01RW.f"
	}
#line 147 "MB01RW.f"
    }

#line 149 "MB01RW.f"
    if (*info != 0) {
#line 150 "MB01RW.f"
	i__1 = -(*info);
#line 150 "MB01RW.f"
	xerbla_("MB01RW", &i__1, (ftnlen)6);
#line 151 "MB01RW.f"
	return 0;
#line 152 "MB01RW.f"
    }

/*     Quick return if possible. */

#line 156 "MB01RW.f"
    if (*n == 0 || *m == 0) {
#line 156 "MB01RW.f"
	return 0;
#line 156 "MB01RW.f"
    }

#line 159 "MB01RW.f"
    if (nottra) {

/*        Compute Z*A*Z'. */

#line 163 "MB01RW.f"
	if (upper) {

/*           Compute Z*A in A (M-by-N). */

#line 167 "MB01RW.f"
	    i__1 = *n;
#line 167 "MB01RW.f"
	    for (j = 1; j <= i__1; ++j) {
#line 168 "MB01RW.f"
		i__2 = j - 1;
#line 168 "MB01RW.f"
		dcopy_(&i__2, &a[j * a_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 169 "MB01RW.f"
		i__2 = *n - j + 1;
#line 169 "MB01RW.f"
		dcopy_(&i__2, &a[j + j * a_dim1], lda, &dwork[j], &c__1);
#line 170 "MB01RW.f"
		dgemv_(trans, m, n, &c_b11, &z__[z_offset], ldz, &dwork[1], &
			c__1, &c_b13, &a[j * a_dim1 + 1], &c__1, (ftnlen)1);
#line 172 "MB01RW.f"
/* L10: */
#line 172 "MB01RW.f"
	    }

/*           Compute A*Z' in the upper triangular part of A. */

#line 176 "MB01RW.f"
	    i__1 = *m;
#line 176 "MB01RW.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 177 "MB01RW.f"
		dcopy_(n, &a[i__ + a_dim1], lda, &dwork[1], &c__1);
#line 178 "MB01RW.f"
		i__2 = *m - i__ + 1;
#line 178 "MB01RW.f"
		dgemv_(trans, &i__2, n, &c_b11, &z__[i__ + z_dim1], ldz, &
			dwork[1], &c__1, &c_b13, &a[i__ + i__ * a_dim1], lda, 
			(ftnlen)1);
#line 180 "MB01RW.f"
/* L20: */
#line 180 "MB01RW.f"
	    }

#line 182 "MB01RW.f"
	} else {

/*           Compute A*Z' in A (N-by-M). */

#line 186 "MB01RW.f"
	    i__1 = *n;
#line 186 "MB01RW.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 187 "MB01RW.f"
		i__2 = i__ - 1;
#line 187 "MB01RW.f"
		dcopy_(&i__2, &a[i__ + a_dim1], lda, &dwork[1], &c__1);
#line 188 "MB01RW.f"
		i__2 = *n - i__ + 1;
#line 188 "MB01RW.f"
		dcopy_(&i__2, &a[i__ + i__ * a_dim1], &c__1, &dwork[i__], &
			c__1);
#line 189 "MB01RW.f"
		dgemv_(trans, m, n, &c_b11, &z__[z_offset], ldz, &dwork[1], &
			c__1, &c_b13, &a[i__ + a_dim1], lda, (ftnlen)1);
#line 191 "MB01RW.f"
/* L30: */
#line 191 "MB01RW.f"
	    }

/*           Compute Z*A in the lower triangular part of A. */

#line 195 "MB01RW.f"
	    i__1 = *m;
#line 195 "MB01RW.f"
	    for (j = 1; j <= i__1; ++j) {
#line 196 "MB01RW.f"
		dcopy_(n, &a[j * a_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 197 "MB01RW.f"
		i__2 = *m - j + 1;
#line 197 "MB01RW.f"
		dgemv_(trans, &i__2, n, &c_b11, &z__[j + z_dim1], ldz, &dwork[
			1], &c__1, &c_b13, &a[j + j * a_dim1], &c__1, (ftnlen)
			1);
#line 199 "MB01RW.f"
/* L40: */
#line 199 "MB01RW.f"
	    }

#line 201 "MB01RW.f"
	}
#line 202 "MB01RW.f"
    } else {

/*        Compute Z'*A*Z. */

#line 206 "MB01RW.f"
	if (upper) {

/*           Compute Z'*A in A (M-by-N). */

#line 210 "MB01RW.f"
	    i__1 = *n;
#line 210 "MB01RW.f"
	    for (j = 1; j <= i__1; ++j) {
#line 211 "MB01RW.f"
		i__2 = j - 1;
#line 211 "MB01RW.f"
		dcopy_(&i__2, &a[j * a_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 212 "MB01RW.f"
		i__2 = *n - j + 1;
#line 212 "MB01RW.f"
		dcopy_(&i__2, &a[j + j * a_dim1], lda, &dwork[j], &c__1);
#line 213 "MB01RW.f"
		dgemv_(trans, n, m, &c_b11, &z__[z_offset], ldz, &dwork[1], &
			c__1, &c_b13, &a[j * a_dim1 + 1], &c__1, (ftnlen)1);
#line 215 "MB01RW.f"
/* L50: */
#line 215 "MB01RW.f"
	    }

/*           Compute A*Z in the upper triangular part of A. */

#line 219 "MB01RW.f"
	    i__1 = *m;
#line 219 "MB01RW.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 220 "MB01RW.f"
		dcopy_(n, &a[i__ + a_dim1], lda, &dwork[1], &c__1);
#line 221 "MB01RW.f"
		i__2 = *m - i__ + 1;
#line 221 "MB01RW.f"
		dgemv_(trans, n, &i__2, &c_b11, &z__[i__ * z_dim1 + 1], ldz, &
			dwork[1], &c__1, &c_b13, &a[i__ + i__ * a_dim1], lda, 
			(ftnlen)1);
#line 223 "MB01RW.f"
/* L60: */
#line 223 "MB01RW.f"
	    }

#line 225 "MB01RW.f"
	} else {

/*           Compute A*Z in A (N-by-M). */

#line 229 "MB01RW.f"
	    i__1 = *n;
#line 229 "MB01RW.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 230 "MB01RW.f"
		i__2 = i__ - 1;
#line 230 "MB01RW.f"
		dcopy_(&i__2, &a[i__ + a_dim1], lda, &dwork[1], &c__1);
#line 231 "MB01RW.f"
		i__2 = *n - i__ + 1;
#line 231 "MB01RW.f"
		dcopy_(&i__2, &a[i__ + i__ * a_dim1], &c__1, &dwork[i__], &
			c__1);
#line 232 "MB01RW.f"
		dgemv_(trans, n, m, &c_b11, &z__[z_offset], ldz, &dwork[1], &
			c__1, &c_b13, &a[i__ + a_dim1], lda, (ftnlen)1);
#line 234 "MB01RW.f"
/* L70: */
#line 234 "MB01RW.f"
	    }

/*           Compute Z'*A in the lower triangular part of A. */

#line 238 "MB01RW.f"
	    i__1 = *m;
#line 238 "MB01RW.f"
	    for (j = 1; j <= i__1; ++j) {
#line 239 "MB01RW.f"
		dcopy_(n, &a[j * a_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 240 "MB01RW.f"
		i__2 = *m - j + 1;
#line 240 "MB01RW.f"
		dgemv_(trans, n, &i__2, &c_b11, &z__[j * z_dim1 + 1], ldz, &
			dwork[1], &c__1, &c_b13, &a[j + j * a_dim1], &c__1, (
			ftnlen)1);
#line 242 "MB01RW.f"
/* L80: */
#line 242 "MB01RW.f"
	    }

#line 244 "MB01RW.f"
	}
#line 245 "MB01RW.f"
    }

#line 247 "MB01RW.f"
    return 0;
/* *** Last line of MB01RW *** */
} /* mb01rw_ */

