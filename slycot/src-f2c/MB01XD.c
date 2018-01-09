#line 1 "MB01XD.f"
/* MB01XD.f -- translated by f2c (version 20100827).
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

#line 1 "MB01XD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b15 = 1.;

/* Subroutine */ int mb01xd_(char *uplo, integer *n, doublereal *a, integer *
	lda, integer *info, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, ib, nb, ii;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb01xy_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dtrmm_(char *, char *, char *, 
	    char *, integer *, integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical upper;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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

/*     To compute the matrix product U' * U or L * L', where U and L are */
/*     upper and lower triangular matrices, respectively, stored in the */
/*     corresponding upper or lower triangular part of the array A. */

/*     If UPLO = 'U' then the upper triangle of the result is stored, */
/*     overwriting the matrix U in A. */
/*     If UPLO = 'L' then the lower triangle of the result is stored, */
/*     overwriting the matrix L in A. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle (U or L) is given in the array A, */
/*             as follows: */
/*             = 'U':  the upper triangular part U is given; */
/*             = 'L':  the lower triangular part L is given. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the triangular matrices U or L.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, if UPLO = 'U', the leading N-by-N upper */
/*             triangular part of this array must contain the upper */
/*             triangular matrix U. */
/*             On entry, if UPLO = 'L', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangular matrix L. */
/*             On exit, if UPLO = 'U', the leading N-by-N upper */
/*             triangular part of this array contains the upper */
/*             triangular part of the product U' * U. The strictly lower */
/*             triangular part is not referenced. */
/*             On exit, if UPLO = 'L', the leading N-by-N lower */
/*             triangular part of this array contains the lower */
/*             triangular part of the product L * L'. The strictly upper */
/*             triangular part is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= max(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix product U' * U or L * L' is computed using BLAS 3 */
/*     operations as much as possible (a block algorithm). */

/*     FURTHER COMMENTS */

/*     This routine is a counterpart of LAPACK Library routine DLAUUM, */
/*     which computes the matrix product U * U' or L' * L. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Nov. 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
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

/*     Test the input scalar arguments. */

#line 127 "MB01XD.f"
    /* Parameter adjustments */
#line 127 "MB01XD.f"
    a_dim1 = *lda;
#line 127 "MB01XD.f"
    a_offset = 1 + a_dim1;
#line 127 "MB01XD.f"
    a -= a_offset;
#line 127 "MB01XD.f"

#line 127 "MB01XD.f"
    /* Function Body */
#line 127 "MB01XD.f"
    *info = 0;
#line 128 "MB01XD.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 129 "MB01XD.f"
    if (! upper && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 130 "MB01XD.f"
	*info = -1;
#line 131 "MB01XD.f"
    } else if (*n < 0) {
#line 132 "MB01XD.f"
	*info = -2;
#line 133 "MB01XD.f"
    } else if (*lda < max(1,*n)) {
#line 134 "MB01XD.f"
	*info = -4;
#line 135 "MB01XD.f"
    }

#line 137 "MB01XD.f"
    if (*info != 0) {

/*        Error return. */

#line 141 "MB01XD.f"
	i__1 = -(*info);
#line 141 "MB01XD.f"
	xerbla_("MB01XD", &i__1, (ftnlen)6);
#line 142 "MB01XD.f"
	return 0;
#line 143 "MB01XD.f"
    }

/*     Quick return, if possible. */

#line 147 "MB01XD.f"
    if (*n == 0) {
#line 147 "MB01XD.f"
	return 0;
#line 147 "MB01XD.f"
    }

/*     Determine the block size for this environment (as for DLAUUM). */

#line 152 "MB01XD.f"
    nb = ilaenv_(&c__1, "DLAUUM", uplo, n, &c_n1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

#line 154 "MB01XD.f"
    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code. */

#line 158 "MB01XD.f"
	mb01xy_(uplo, n, &a[a_offset], lda, info, (ftnlen)1);
#line 159 "MB01XD.f"
    } else {

/*        Use blocked code. */

#line 163 "MB01XD.f"
	if (upper) {

/*           Compute the product U' * U. */

#line 167 "MB01XD.f"
	    i__1 = -nb;
#line 167 "MB01XD.f"
	    for (i__ = *n; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
#line 168 "MB01XD.f"
		ib = min(nb,i__);
#line 169 "MB01XD.f"
		ii = i__ - ib + 1;
#line 170 "MB01XD.f"
		if (i__ < *n) {
#line 171 "MB01XD.f"
		    i__2 = *n - i__;
#line 171 "MB01XD.f"
		    dtrmm_("Left", "Upper", "Transpose", "Non-unit", &ib, &
			    i__2, &c_b15, &a[ii + ii * a_dim1], lda, &a[ii + (
			    ii + ib) * a_dim1], lda, (ftnlen)4, (ftnlen)5, (
			    ftnlen)9, (ftnlen)8);
#line 174 "MB01XD.f"
		    i__2 = *n - i__;
#line 174 "MB01XD.f"
		    i__3 = i__ - ib;
#line 174 "MB01XD.f"
		    dgemm_("Transpose", "No transpose", &ib, &i__2, &i__3, &
			    c_b15, &a[ii * a_dim1 + 1], lda, &a[(ii + ib) * 
			    a_dim1 + 1], lda, &c_b15, &a[ii + (ii + ib) * 
			    a_dim1], lda, (ftnlen)9, (ftnlen)12);
#line 177 "MB01XD.f"
		}
#line 178 "MB01XD.f"
		mb01xy_("Upper", &ib, &a[ii + ii * a_dim1], lda, info, (
			ftnlen)5);
#line 179 "MB01XD.f"
		i__2 = ii - 1;
#line 179 "MB01XD.f"
		dsyrk_("Upper", "Transpose", &ib, &i__2, &c_b15, &a[ii * 
			a_dim1 + 1], lda, &c_b15, &a[ii + ii * a_dim1], lda, (
			ftnlen)5, (ftnlen)9);
#line 181 "MB01XD.f"
/* L10: */
#line 181 "MB01XD.f"
	    }
#line 182 "MB01XD.f"
	} else {

/*           Compute the product L * L'. */

#line 186 "MB01XD.f"
	    i__1 = -nb;
#line 186 "MB01XD.f"
	    for (i__ = *n; i__1 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__1) {
#line 187 "MB01XD.f"
		ib = min(nb,i__);
#line 188 "MB01XD.f"
		ii = i__ - ib + 1;
#line 189 "MB01XD.f"
		if (i__ < *n) {
#line 190 "MB01XD.f"
		    i__2 = *n - i__;
#line 190 "MB01XD.f"
		    dtrmm_("Right", "Lower", "Transpose", "Non-unit", &i__2, &
			    ib, &c_b15, &a[ii + ii * a_dim1], lda, &a[ii + ib 
			    + ii * a_dim1], lda, (ftnlen)5, (ftnlen)5, (
			    ftnlen)9, (ftnlen)8);
#line 193 "MB01XD.f"
		    i__2 = *n - i__;
#line 193 "MB01XD.f"
		    i__3 = i__ - ib;
#line 193 "MB01XD.f"
		    dgemm_("No transpose", "Transpose", &i__2, &ib, &i__3, &
			    c_b15, &a[ii + ib + a_dim1], lda, &a[ii + a_dim1],
			     lda, &c_b15, &a[ii + ib + ii * a_dim1], lda, (
			    ftnlen)12, (ftnlen)9);
#line 196 "MB01XD.f"
		}
#line 197 "MB01XD.f"
		mb01xy_("Lower", &ib, &a[ii + ii * a_dim1], lda, info, (
			ftnlen)5);
#line 198 "MB01XD.f"
		i__2 = ii - 1;
#line 198 "MB01XD.f"
		dsyrk_("Lower", "No Transpose", &ib, &i__2, &c_b15, &a[ii + 
			a_dim1], lda, &c_b15, &a[ii + ii * a_dim1], lda, (
			ftnlen)5, (ftnlen)12);
#line 200 "MB01XD.f"
/* L20: */
#line 200 "MB01XD.f"
	    }
#line 201 "MB01XD.f"
	}
#line 202 "MB01XD.f"
    }

#line 204 "MB01XD.f"
    return 0;

/* *** Last line of MB01XD *** */
} /* mb01xd_ */
