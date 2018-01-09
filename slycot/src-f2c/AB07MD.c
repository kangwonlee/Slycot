#line 1 "AB07MD.f"
/* AB07MD.f -- translated by f2c (version 20100827).
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

#line 1 "AB07MD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int ab07md_(char *jobd, integer *n, integer *m, integer *p, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, integer *info, 
	ftnlen jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2;

    /* Local variables */
    static integer j;
    static logical ljobd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer mplim;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer minmp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);


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

/*     To find the dual of a given state-space representation. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the dual state dynamics matrix A'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,MAX(M,P)) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */
/*             On exit, the leading N-by-P part of this array contains */
/*             the dual input/state matrix C'. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the dual state/output matrix B'. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,M,P) if N > 0. */
/*             LDC >= 1 if N = 0. */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDD,MAX(M,P)) */
/*             On entry, if JOBD = 'D', the leading P-by-M part of this */
/*             array must contain the original direct transmission */
/*             matrix D. */
/*             On exit, if JOBD = 'D', the leading M-by-P part of this */
/*             array contains the dual direct transmission matrix D'. */
/*             The array D is not referenced if JOBD = 'Z'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,M,P) if JOBD = 'D'. */
/*             LDD >= 1 if JOBD = 'Z'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If the given state-space representation is the M-input/P-output */
/*     (A,B,C,D), its dual is simply the P-input/M-output (A',C',B',D'). */

/*     REFERENCES */

/*     None */

/*     NUMERICAL ASPECTS */

/*     None */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine AB07AD by T.W.C.Williams, Kingston */
/*     Polytechnic, United Kingdom, March 1982. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

/*     KEYWORDS */

/*     Dual system, state-space model, state-space representation. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External functions .. */
/*     .. External subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 144 "AB07MD.f"
    /* Parameter adjustments */
#line 144 "AB07MD.f"
    a_dim1 = *lda;
#line 144 "AB07MD.f"
    a_offset = 1 + a_dim1;
#line 144 "AB07MD.f"
    a -= a_offset;
#line 144 "AB07MD.f"
    b_dim1 = *ldb;
#line 144 "AB07MD.f"
    b_offset = 1 + b_dim1;
#line 144 "AB07MD.f"
    b -= b_offset;
#line 144 "AB07MD.f"
    c_dim1 = *ldc;
#line 144 "AB07MD.f"
    c_offset = 1 + c_dim1;
#line 144 "AB07MD.f"
    c__ -= c_offset;
#line 144 "AB07MD.f"
    d_dim1 = *ldd;
#line 144 "AB07MD.f"
    d_offset = 1 + d_dim1;
#line 144 "AB07MD.f"
    d__ -= d_offset;
#line 144 "AB07MD.f"

#line 144 "AB07MD.f"
    /* Function Body */
#line 144 "AB07MD.f"
    *info = 0;
#line 145 "AB07MD.f"
    ljobd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
#line 146 "AB07MD.f"
    mplim = max(*m,*p);
#line 147 "AB07MD.f"
    minmp = min(*m,*p);

/*     Test the input scalar arguments. */

#line 151 "AB07MD.f"
    if (! ljobd && ! lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1)) {
#line 152 "AB07MD.f"
	*info = -1;
#line 153 "AB07MD.f"
    } else if (*n < 0) {
#line 154 "AB07MD.f"
	*info = -2;
#line 155 "AB07MD.f"
    } else if (*m < 0) {
#line 156 "AB07MD.f"
	*info = -3;
#line 157 "AB07MD.f"
    } else if (*p < 0) {
#line 158 "AB07MD.f"
	*info = -4;
#line 159 "AB07MD.f"
    } else if (*lda < max(1,*n)) {
#line 160 "AB07MD.f"
	*info = -6;
#line 161 "AB07MD.f"
    } else if (*ldb < max(1,*n)) {
#line 162 "AB07MD.f"
	*info = -8;
#line 163 "AB07MD.f"
    } else if (*n > 0 && *ldc < max(1,mplim) || *n == 0 && *ldc < 1) {
#line 165 "AB07MD.f"
	*info = -10;
#line 166 "AB07MD.f"
    } else if (ljobd && *ldd < max(1,mplim) || ! ljobd && *ldd < 1) {
#line 168 "AB07MD.f"
	*info = -12;
#line 169 "AB07MD.f"
    }

#line 171 "AB07MD.f"
    if (*info != 0) {

/*        Error return. */

#line 175 "AB07MD.f"
	i__1 = -(*info);
#line 175 "AB07MD.f"
	xerbla_("AB07MD", &i__1, (ftnlen)6);
#line 176 "AB07MD.f"
	return 0;
#line 177 "AB07MD.f"
    }

/*     Quick return if possible. */

#line 181 "AB07MD.f"
    if (max(*n,minmp) == 0) {
#line 181 "AB07MD.f"
	return 0;
#line 181 "AB07MD.f"
    }

#line 184 "AB07MD.f"
    if (*n > 0) {

/*        Transpose A, if non-scalar. */

#line 188 "AB07MD.f"
	i__1 = *n - 1;
#line 188 "AB07MD.f"
	for (j = 1; j <= i__1; ++j) {
#line 189 "AB07MD.f"
	    i__2 = *n - j;
#line 189 "AB07MD.f"
	    dswap_(&i__2, &a[j + 1 + j * a_dim1], &c__1, &a[j + (j + 1) * 
		    a_dim1], lda);
#line 190 "AB07MD.f"
/* L10: */
#line 190 "AB07MD.f"
	}

/*        Replace B by C' and C by B'. */

#line 194 "AB07MD.f"
	i__1 = mplim;
#line 194 "AB07MD.f"
	for (j = 1; j <= i__1; ++j) {
#line 195 "AB07MD.f"
	    if (j <= minmp) {
#line 196 "AB07MD.f"
		dswap_(n, &b[j * b_dim1 + 1], &c__1, &c__[j + c_dim1], ldc);
#line 197 "AB07MD.f"
	    } else if (j > *p) {
#line 198 "AB07MD.f"
		dcopy_(n, &b[j * b_dim1 + 1], &c__1, &c__[j + c_dim1], ldc);
#line 199 "AB07MD.f"
	    } else {
#line 200 "AB07MD.f"
		dcopy_(n, &c__[j + c_dim1], ldc, &b[j * b_dim1 + 1], &c__1);
#line 201 "AB07MD.f"
	    }
#line 202 "AB07MD.f"
/* L20: */
#line 202 "AB07MD.f"
	}

#line 204 "AB07MD.f"
    }

#line 206 "AB07MD.f"
    if (ljobd && minmp > 0) {

/*        Transpose D, if non-scalar. */

#line 210 "AB07MD.f"
	i__1 = mplim;
#line 210 "AB07MD.f"
	for (j = 1; j <= i__1; ++j) {
#line 211 "AB07MD.f"
	    if (j < minmp) {
#line 212 "AB07MD.f"
		i__2 = minmp - j;
#line 212 "AB07MD.f"
		dswap_(&i__2, &d__[j + 1 + j * d_dim1], &c__1, &d__[j + (j + 
			1) * d_dim1], ldd);
#line 213 "AB07MD.f"
	    } else if (j > *p) {
#line 214 "AB07MD.f"
		dcopy_(p, &d__[j * d_dim1 + 1], &c__1, &d__[j + d_dim1], ldd);
#line 215 "AB07MD.f"
	    } else if (j > *m) {
#line 216 "AB07MD.f"
		dcopy_(m, &d__[j + d_dim1], ldd, &d__[j * d_dim1 + 1], &c__1);
#line 217 "AB07MD.f"
	    }
#line 218 "AB07MD.f"
/* L30: */
#line 218 "AB07MD.f"
	}

#line 220 "AB07MD.f"
    }

#line 222 "AB07MD.f"
    return 0;
/* *** Last line of AB07MD *** */
} /* ab07md_ */
