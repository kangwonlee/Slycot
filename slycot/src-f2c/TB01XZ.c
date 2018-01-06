#line 1 "TB01XZ.f"
/* TB01XZ.f -- translated by f2c (version 20100827).
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

#line 1 "TB01XZ.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int tb01xz_(char *jobd, integer *n, integer *m, integer *p, 
	integer *kl, integer *ku, doublecomplex *a, integer *lda, 
	doublecomplex *b, integer *ldb, doublecomplex *c__, integer *ldc, 
	doublecomplex *d__, integer *ldd, integer *info, ftnlen jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, j1, nm1, lda1;
    static logical ljobd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer minmp, maxmp;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), xerbla_(char *, integer *,
	     ftnlen);


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

/*     To apply a special transformation to a system given as a triple */
/*     (A,B,C), */

/*        A <-- P * A' * P,  B <-- P * C',  C <-- B' * P, */

/*     where P is a matrix with 1 on the secondary diagonal, and with 0 */
/*     in the other entries. Matrix A can be specified as a band matrix. */
/*     Optionally, matrix D of the system can be transposed. This */
/*     transformation is actually a special similarity transformation of */
/*     the dual system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A, the number of rows of matrix B */
/*             and the number of columns of matrix C. */
/*             N represents the dimension of the state vector.  N >= 0. */

/*     M       (input) INTEGER. */
/*             The number of columns of matrix B. */
/*             M represents the dimension of input vector.  M >= 0. */

/*     P       (input) INTEGER. */
/*             The number of rows of matrix C. */
/*             P represents the dimension of output vector.  P >= 0. */

/*     KL      (input) INTEGER */
/*             The number of subdiagonals of A to be transformed. */
/*             MAX( 0, N-1 ) >= KL >= 0. */

/*     KU      (input) INTEGER */
/*             The number of superdiagonals of A to be transformed. */
/*             MAX( 0, N-1 ) >= KU >= 0. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed (pertransposed) matrix P*A'*P. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) COMPLEX*16 array, dimension (LDB,MAX(M,P)) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */
/*             On exit, the leading N-by-P part of this array contains */
/*             the dual input/state matrix P*C'. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N) if M > 0 or  P > 0. */
/*             LDB >= 1        if M = 0 and P = 0. */

/*     C       (input/output) COMPLEX*16 array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the dual state/output matrix B'*P. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,M,P) if N > 0. */
/*             LDC >= 1          if N = 0. */

/*     D       (input/output) COMPLEX*16 array, dimension (LDD,MAX(M,P)) */
/*             On entry, if JOBD = 'D', the leading P-by-M part of this */
/*             array must contain the original direct transmission */
/*             matrix D. */
/*             On exit, if JOBD = 'D', the leading M-by-P part of this */
/*             array contains the transposed direct transmission matrix */
/*             D'. The array D is not referenced if JOBD = 'Z'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,M,P) if JOBD = 'D'. */
/*             LDD >= 1          if JOBD = 'Z'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The rows and/or columns of the matrices of the triplet (A,B,C) */
/*     and, optionally, of the matrix D are swapped in a special way. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Matrix algebra, matrix operations, similarity transformation. */

/*  ********************************************************************* */

/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the scalar input arguments. */

#line 171 "TB01XZ.f"
    /* Parameter adjustments */
#line 171 "TB01XZ.f"
    a_dim1 = *lda;
#line 171 "TB01XZ.f"
    a_offset = 1 + a_dim1;
#line 171 "TB01XZ.f"
    a -= a_offset;
#line 171 "TB01XZ.f"
    b_dim1 = *ldb;
#line 171 "TB01XZ.f"
    b_offset = 1 + b_dim1;
#line 171 "TB01XZ.f"
    b -= b_offset;
#line 171 "TB01XZ.f"
    c_dim1 = *ldc;
#line 171 "TB01XZ.f"
    c_offset = 1 + c_dim1;
#line 171 "TB01XZ.f"
    c__ -= c_offset;
#line 171 "TB01XZ.f"
    d_dim1 = *ldd;
#line 171 "TB01XZ.f"
    d_offset = 1 + d_dim1;
#line 171 "TB01XZ.f"
    d__ -= d_offset;
#line 171 "TB01XZ.f"

#line 171 "TB01XZ.f"
    /* Function Body */
#line 171 "TB01XZ.f"
    *info = 0;
#line 172 "TB01XZ.f"
    ljobd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
#line 173 "TB01XZ.f"
    maxmp = max(*m,*p);
#line 174 "TB01XZ.f"
    minmp = min(*m,*p);
#line 175 "TB01XZ.f"
    nm1 = *n - 1;

#line 177 "TB01XZ.f"
    if (! ljobd && ! lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1)) {
#line 178 "TB01XZ.f"
	*info = -1;
#line 179 "TB01XZ.f"
    } else if (*n < 0) {
#line 180 "TB01XZ.f"
	*info = -2;
#line 181 "TB01XZ.f"
    } else if (*m < 0) {
#line 182 "TB01XZ.f"
	*info = -3;
#line 183 "TB01XZ.f"
    } else if (*p < 0) {
#line 184 "TB01XZ.f"
	*info = -4;
#line 185 "TB01XZ.f"
    } else if (*kl < 0 || *kl > max(0,nm1)) {
#line 186 "TB01XZ.f"
	*info = -5;
#line 187 "TB01XZ.f"
    } else if (*ku < 0 || *ku > max(0,nm1)) {
#line 188 "TB01XZ.f"
	*info = -6;
#line 189 "TB01XZ.f"
    } else if (*lda < max(1,*n)) {
#line 190 "TB01XZ.f"
	*info = -8;
#line 191 "TB01XZ.f"
    } else if (maxmp > 0 && *ldb < max(1,*n) || minmp == 0 && *ldb < 1) {
#line 193 "TB01XZ.f"
	*info = -10;
#line 194 "TB01XZ.f"
    } else if (*ldc < 1 || *n > 0 && *ldc < maxmp) {
#line 195 "TB01XZ.f"
	*info = -12;
#line 196 "TB01XZ.f"
    } else if (*ldd < 1 || ljobd && *ldd < maxmp) {
#line 197 "TB01XZ.f"
	*info = -14;
#line 198 "TB01XZ.f"
    }

#line 200 "TB01XZ.f"
    if (*info != 0) {

/*        Error return. */

#line 204 "TB01XZ.f"
	i__1 = -(*info);
#line 204 "TB01XZ.f"
	xerbla_("TB01XZ", &i__1, (ftnlen)6);
#line 205 "TB01XZ.f"
	return 0;
#line 206 "TB01XZ.f"
    }

/*     Quick return if possible. */

#line 210 "TB01XZ.f"
    if (ljobd) {

/*        Replace D by D', if non-scalar. */

#line 214 "TB01XZ.f"
	i__1 = maxmp;
#line 214 "TB01XZ.f"
	for (j = 1; j <= i__1; ++j) {
#line 215 "TB01XZ.f"
	    if (j < minmp) {
#line 216 "TB01XZ.f"
		i__2 = minmp - j;
#line 216 "TB01XZ.f"
		zswap_(&i__2, &d__[j + 1 + j * d_dim1], &c__1, &d__[j + (j + 
			1) * d_dim1], ldd);
#line 217 "TB01XZ.f"
	    } else if (j > *p) {
#line 218 "TB01XZ.f"
		zcopy_(p, &d__[j * d_dim1 + 1], &c__1, &d__[j + d_dim1], ldd);
#line 219 "TB01XZ.f"
	    } else if (j > *m) {
#line 220 "TB01XZ.f"
		zcopy_(m, &d__[j + d_dim1], ldd, &d__[j * d_dim1 + 1], &c__1);
#line 221 "TB01XZ.f"
	    }
#line 222 "TB01XZ.f"
/* L5: */
#line 222 "TB01XZ.f"
	}

#line 224 "TB01XZ.f"
    }

#line 226 "TB01XZ.f"
    if (*n == 0) {
#line 226 "TB01XZ.f"
	return 0;
#line 226 "TB01XZ.f"
    }

/*     Replace matrix A by P*A'*P. */

#line 231 "TB01XZ.f"
    if (*kl == nm1 && *ku == nm1) {

/*        Full matrix A. */

#line 235 "TB01XZ.f"
	i__1 = nm1;
#line 235 "TB01XZ.f"
	for (j = 1; j <= i__1; ++j) {
#line 236 "TB01XZ.f"
	    i__2 = *n - j;
#line 236 "TB01XZ.f"
	    i__3 = -(*lda);
#line 236 "TB01XZ.f"
	    zswap_(&i__2, &a[j * a_dim1 + 1], &c__1, &a[*n - j + 1 + (j + 1) *
		     a_dim1], &i__3);
#line 237 "TB01XZ.f"
/* L10: */
#line 237 "TB01XZ.f"
	}

#line 239 "TB01XZ.f"
    } else {

/*        Band matrix A. */

#line 243 "TB01XZ.f"
	lda1 = *lda + 1;

/*        Pertranspose the KL subdiagonals. */

/* Computing MIN */
#line 247 "TB01XZ.f"
	i__2 = *kl, i__3 = *n - 2;
#line 247 "TB01XZ.f"
	i__1 = min(i__2,i__3);
#line 247 "TB01XZ.f"
	for (j = 1; j <= i__1; ++j) {
#line 248 "TB01XZ.f"
	    j1 = (*n - j) / 2;
#line 249 "TB01XZ.f"
	    i__2 = -lda1;
#line 249 "TB01XZ.f"
	    zswap_(&j1, &a[j + 1 + a_dim1], &lda1, &a[*n - j1 + 1 + (*n - j1 
		    + 1 - j) * a_dim1], &i__2);
#line 250 "TB01XZ.f"
/* L20: */
#line 250 "TB01XZ.f"
	}

/*        Pertranspose the KU superdiagonals. */

/* Computing MIN */
#line 254 "TB01XZ.f"
	i__2 = *ku, i__3 = *n - 2;
#line 254 "TB01XZ.f"
	i__1 = min(i__2,i__3);
#line 254 "TB01XZ.f"
	for (j = 1; j <= i__1; ++j) {
#line 255 "TB01XZ.f"
	    j1 = (*n - j) / 2;
#line 256 "TB01XZ.f"
	    i__2 = -lda1;
#line 256 "TB01XZ.f"
	    zswap_(&j1, &a[(j + 1) * a_dim1 + 1], &lda1, &a[*n - j1 + 1 - j + 
		    (*n - j1 + 1) * a_dim1], &i__2);
#line 257 "TB01XZ.f"
/* L30: */
#line 257 "TB01XZ.f"
	}

/*        Pertranspose the diagonal. */

#line 261 "TB01XZ.f"
	j1 = *n / 2;
#line 262 "TB01XZ.f"
	i__1 = -lda1;
#line 262 "TB01XZ.f"
	zswap_(&j1, &a[a_dim1 + 1], &lda1, &a[*n - j1 + 1 + (*n - j1 + 1) * 
		a_dim1], &i__1);

#line 264 "TB01XZ.f"
    }

/*     Replace matrix B by P*C' and matrix C by B'*P. */

#line 268 "TB01XZ.f"
    i__1 = maxmp;
#line 268 "TB01XZ.f"
    for (j = 1; j <= i__1; ++j) {
#line 269 "TB01XZ.f"
	if (j <= minmp) {
#line 270 "TB01XZ.f"
	    i__2 = -(*ldc);
#line 270 "TB01XZ.f"
	    zswap_(n, &b[j * b_dim1 + 1], &c__1, &c__[j + c_dim1], &i__2);
#line 271 "TB01XZ.f"
	} else if (j > *p) {
#line 272 "TB01XZ.f"
	    i__2 = -(*ldc);
#line 272 "TB01XZ.f"
	    zcopy_(n, &b[j * b_dim1 + 1], &c__1, &c__[j + c_dim1], &i__2);
#line 273 "TB01XZ.f"
	} else {
#line 274 "TB01XZ.f"
	    i__2 = -(*ldc);
#line 274 "TB01XZ.f"
	    zcopy_(n, &c__[j + c_dim1], &i__2, &b[j * b_dim1 + 1], &c__1);
#line 275 "TB01XZ.f"
	}
#line 276 "TB01XZ.f"
/* L40: */
#line 276 "TB01XZ.f"
    }

#line 278 "TB01XZ.f"
    return 0;
/* *** Last line of TB01XZ *** */
} /* tb01xz_ */

