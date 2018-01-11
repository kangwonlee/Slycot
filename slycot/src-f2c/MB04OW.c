#line 1 "MB04OW.f"
/* MB04OW.f -- translated by f2c (version 20100827).
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

#line 1 "MB04OW.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb04ow_(integer *m, integer *n, integer *p, doublereal *
	a, integer *lda, doublereal *t, integer *ldt, doublereal *x, integer *
	incx, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *incd)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, t_dim1, 
	    t_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal ci;
    static integer mn;
    static doublereal si;
    static integer ix;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlartg_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);


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

/*     To perform the QR factorization */

/*        ( U  ) = Q*( R ),  where  U = ( U1  U2 ),  R = ( R1  R2 ), */
/*        ( x' )     ( 0 )              ( 0   T  )       ( 0   R3 ) */

/*     where U and R are (m+n)-by-(m+n) upper triangular matrices, x is */
/*     an m+n element vector, U1 is m-by-m, T is n-by-n, stored */
/*     separately, and Q is an (m+n+1)-by-(m+n+1) orthogonal matrix. */

/*     The matrix ( U1 U2 ) must be supplied in the m-by-(m+n) upper */
/*     trapezoidal part of the array A and this is overwritten by the */
/*     corresponding part ( R1 R2 ) of R. The remaining upper triangular */
/*     part of R, R3, is overwritten on the array T. */

/*     The transformations performed are also applied to the (m+n+1)-by-p */
/*     matrix ( B' C' d )' (' denotes transposition), where B, C, and d' */
/*     are m-by-p, n-by-p, and 1-by-p matrices, respectively. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M      (input) INTEGER */
/*            The number of rows of the matrix ( U1  U2 ).  M >= 0. */

/*     N      (input) INTEGER */
/*            The order of the matrix T.  N >= 0. */

/*     P      (input) INTEGER */
/*            The number of columns of the matrices B and C.  P >= 0. */

/*     A      (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*            On entry, the leading M-by-(M+N) upper trapezoidal part of */
/*            this array must contain the upper trapezoidal matrix */
/*            ( U1 U2 ). */
/*            On exit, the leading M-by-(M+N) upper trapezoidal part of */
/*            this array contains the upper trapezoidal matrix ( R1 R2 ). */
/*            The strict lower triangle of A is not referenced. */

/*     LDA    INTEGER */
/*            The leading dimension of the array A.  LDA >= max(1,M). */

/*     T      (input/output) DOUBLE PRECISION array, dimension (LDT,N) */
/*            On entry, the leading N-by-N upper triangular part of this */
/*            array must contain the upper triangular matrix T. */
/*            On exit, the leading N-by-N upper triangular part of this */
/*            array contains the upper triangular matrix R3. */
/*            The strict lower triangle of T is not referenced. */

/*     LDT    INTEGER */
/*            The leading dimension of the array T.  LDT >= max(1,N). */

/*     X      (input/output) DOUBLE PRECISION array, dimension */
/*            (1+(M+N-1)*INCX), if M+N > 0, or dimension (0), if M+N = 0. */
/*            On entry, the incremented array X must contain the */
/*            vector x. On exit, the content of X is changed. */

/*     INCX   (input) INTEGER */
/*            Specifies the increment for the elements of X.  INCX > 0. */

/*     B      (input/output) DOUBLE PRECISION array, dimension (LDB,P) */
/*            On entry, the leading M-by-P part of this array must */
/*            contain the matrix B. */
/*            On exit, the leading M-by-P part of this array contains */
/*            the transformed matrix B. */
/*            If M = 0 or P = 0, this array is not referenced. */

/*     LDB    INTEGER */
/*            The leading dimension of the array B. */
/*            LDB >= max(1,M), if P > 0; */
/*            LDB >= 1,        if P = 0. */

/*     C      (input/output) DOUBLE PRECISION array, dimension (LDC,P) */
/*            On entry, the leading N-by-P part of this array must */
/*            contain the matrix C. */
/*            On exit, the leading N-by-P part of this array contains */
/*            the transformed matrix C. */
/*            If N = 0 or P = 0, this array is not referenced. */

/*     LDC    INTEGER */
/*            The leading dimension of the array C. */
/*            LDC >= max(1,N), if P > 0; */
/*            LDC >= 1,        if P = 0. */

/*     D      (input/output) DOUBLE PRECISION array, dimension */
/*            (1+(P-1)*INCD), if P > 0, or dimension (0), if P = 0. */
/*            On entry, the incremented array D must contain the */
/*            vector d. */
/*            On exit, this incremented array contains the transformed */
/*            vector d. */
/*            If P = 0, this array is not referenced. */

/*     INCD   (input) INTEGER */
/*            Specifies the increment for the elements of D.  INCD > 0. */

/*     METHOD */

/*     Let q = m+n. The matrix Q is formed as a sequence of plane */
/*     rotations in planes (1, q+1), (2, q+1), ..., (q, q+1), the */
/*     rotation in the (j, q+1)th plane, Q(j), being chosen to */
/*     annihilate the jth element of x. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0((M+N)*(M+N+P)) operations and is backward */
/*     stable. */

/*     FURTHER COMMENTS */

/*     For P = 0, this routine produces the same result as SLICOT Library */
/*     routine MB04OX, but matrix T may not be stored in the array A. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Matrix operations, plane rotations. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

/*     For efficiency reasons, the parameters are not checked. */

#line 165 "MB04OW.f"
    /* Parameter adjustments */
#line 165 "MB04OW.f"
    a_dim1 = *lda;
#line 165 "MB04OW.f"
    a_offset = 1 + a_dim1;
#line 165 "MB04OW.f"
    a -= a_offset;
#line 165 "MB04OW.f"
    t_dim1 = *ldt;
#line 165 "MB04OW.f"
    t_offset = 1 + t_dim1;
#line 165 "MB04OW.f"
    t -= t_offset;
#line 165 "MB04OW.f"
    --x;
#line 165 "MB04OW.f"
    b_dim1 = *ldb;
#line 165 "MB04OW.f"
    b_offset = 1 + b_dim1;
#line 165 "MB04OW.f"
    b -= b_offset;
#line 165 "MB04OW.f"
    c_dim1 = *ldc;
#line 165 "MB04OW.f"
    c_offset = 1 + c_dim1;
#line 165 "MB04OW.f"
    c__ -= c_offset;
#line 165 "MB04OW.f"
    --d__;
#line 165 "MB04OW.f"

#line 165 "MB04OW.f"
    /* Function Body */
#line 165 "MB04OW.f"
    mn = *m + *n;
#line 166 "MB04OW.f"
    if (*incx > 1) {

/*        Code for increment INCX > 1. */

#line 170 "MB04OW.f"
	ix = 1;
#line 171 "MB04OW.f"
	if (*m > 0) {

#line 173 "MB04OW.f"
	    i__1 = *m - 1;
#line 173 "MB04OW.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 174 "MB04OW.f"
		dlartg_(&a[i__ + i__ * a_dim1], &x[ix], &ci, &si, &temp);
#line 175 "MB04OW.f"
		a[i__ + i__ * a_dim1] = temp;
#line 176 "MB04OW.f"
		ix += *incx;
#line 177 "MB04OW.f"
		i__2 = mn - i__;
#line 177 "MB04OW.f"
		drot_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda, &x[ix], incx, 
			&ci, &si);
#line 178 "MB04OW.f"
		if (*p > 0) {
#line 178 "MB04OW.f"
		    drot_(p, &b[i__ + b_dim1], ldb, &d__[1], incd, &ci, &si);
#line 178 "MB04OW.f"
		}
#line 180 "MB04OW.f"
/* L10: */
#line 180 "MB04OW.f"
	    }

#line 182 "MB04OW.f"
	    dlartg_(&a[*m + *m * a_dim1], &x[ix], &ci, &si, &temp);
#line 183 "MB04OW.f"
	    a[*m + *m * a_dim1] = temp;
#line 184 "MB04OW.f"
	    ix += *incx;
#line 185 "MB04OW.f"
	    if (*n > 0) {
#line 185 "MB04OW.f"
		drot_(n, &a[*m + (*m + 1) * a_dim1], lda, &x[ix], incx, &ci, &
			si);
#line 185 "MB04OW.f"
	    }
#line 187 "MB04OW.f"
	    if (*p > 0) {
#line 187 "MB04OW.f"
		drot_(p, &b[*m + b_dim1], ldb, &d__[1], incd, &ci, &si);
#line 187 "MB04OW.f"
	    }
#line 189 "MB04OW.f"
	}

#line 191 "MB04OW.f"
	if (*n > 0) {

#line 193 "MB04OW.f"
	    i__1 = *n - 1;
#line 193 "MB04OW.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 194 "MB04OW.f"
		dlartg_(&t[i__ + i__ * t_dim1], &x[ix], &ci, &si, &temp);
#line 195 "MB04OW.f"
		t[i__ + i__ * t_dim1] = temp;
#line 196 "MB04OW.f"
		ix += *incx;
#line 197 "MB04OW.f"
		i__2 = *n - i__;
#line 197 "MB04OW.f"
		drot_(&i__2, &t[i__ + (i__ + 1) * t_dim1], ldt, &x[ix], incx, 
			&ci, &si);
#line 198 "MB04OW.f"
		if (*p > 0) {
#line 198 "MB04OW.f"
		    drot_(p, &c__[i__ + c_dim1], ldc, &d__[1], incd, &ci, &si)
			    ;
#line 198 "MB04OW.f"
		}
#line 200 "MB04OW.f"
/* L20: */
#line 200 "MB04OW.f"
	    }

#line 202 "MB04OW.f"
	    dlartg_(&t[*n + *n * t_dim1], &x[ix], &ci, &si, &temp);
#line 203 "MB04OW.f"
	    t[*n + *n * t_dim1] = temp;
#line 204 "MB04OW.f"
	    if (*p > 0) {
#line 204 "MB04OW.f"
		drot_(p, &c__[*n + c_dim1], ldc, &d__[1], incd, &ci, &si);
#line 204 "MB04OW.f"
	    }
#line 206 "MB04OW.f"
	}

#line 208 "MB04OW.f"
    } else if (*incx == 1) {

/*        Code for increment INCX = 1. */

#line 212 "MB04OW.f"
	if (*m > 0) {

#line 214 "MB04OW.f"
	    i__1 = *m - 1;
#line 214 "MB04OW.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 215 "MB04OW.f"
		dlartg_(&a[i__ + i__ * a_dim1], &x[i__], &ci, &si, &temp);
#line 216 "MB04OW.f"
		a[i__ + i__ * a_dim1] = temp;
#line 217 "MB04OW.f"
		i__2 = mn - i__;
#line 217 "MB04OW.f"
		drot_(&i__2, &a[i__ + (i__ + 1) * a_dim1], lda, &x[i__ + 1], &
			c__1, &ci, &si);
#line 218 "MB04OW.f"
		if (*p > 0) {
#line 218 "MB04OW.f"
		    drot_(p, &b[i__ + b_dim1], ldb, &d__[1], incd, &ci, &si);
#line 218 "MB04OW.f"
		}
#line 220 "MB04OW.f"
/* L30: */
#line 220 "MB04OW.f"
	    }

#line 222 "MB04OW.f"
	    dlartg_(&a[*m + *m * a_dim1], &x[*m], &ci, &si, &temp);
#line 223 "MB04OW.f"
	    a[*m + *m * a_dim1] = temp;
#line 224 "MB04OW.f"
	    if (*n > 0) {
#line 224 "MB04OW.f"
		drot_(n, &a[*m + (*m + 1) * a_dim1], lda, &x[*m + 1], &c__1, &
			ci, &si);
#line 224 "MB04OW.f"
	    }
#line 226 "MB04OW.f"
	    if (*p > 0) {
#line 226 "MB04OW.f"
		drot_(p, &b[*m + b_dim1], ldb, &d__[1], incd, &ci, &si);
#line 226 "MB04OW.f"
	    }
#line 228 "MB04OW.f"
	}

#line 230 "MB04OW.f"
	if (*n > 0) {
#line 231 "MB04OW.f"
	    ix = *m + 1;

#line 233 "MB04OW.f"
	    i__1 = *n - 1;
#line 233 "MB04OW.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 234 "MB04OW.f"
		dlartg_(&t[i__ + i__ * t_dim1], &x[ix], &ci, &si, &temp);
#line 235 "MB04OW.f"
		t[i__ + i__ * t_dim1] = temp;
#line 236 "MB04OW.f"
		++ix;
#line 237 "MB04OW.f"
		i__2 = *n - i__;
#line 237 "MB04OW.f"
		drot_(&i__2, &t[i__ + (i__ + 1) * t_dim1], ldt, &x[ix], &c__1,
			 &ci, &si);
#line 238 "MB04OW.f"
		if (*p > 0) {
#line 238 "MB04OW.f"
		    drot_(p, &c__[i__ + c_dim1], ldc, &d__[1], incd, &ci, &si)
			    ;
#line 238 "MB04OW.f"
		}
#line 240 "MB04OW.f"
/* L40: */
#line 240 "MB04OW.f"
	    }

#line 242 "MB04OW.f"
	    dlartg_(&t[*n + *n * t_dim1], &x[ix], &ci, &si, &temp);
#line 243 "MB04OW.f"
	    t[*n + *n * t_dim1] = temp;
#line 244 "MB04OW.f"
	    if (*p > 0) {
#line 244 "MB04OW.f"
		drot_(p, &c__[*n + c_dim1], ldc, &d__[1], incd, &ci, &si);
#line 244 "MB04OW.f"
	    }
#line 246 "MB04OW.f"
	}
#line 247 "MB04OW.f"
    }

#line 249 "MB04OW.f"
    return 0;
/* *** Last line of MB04OW *** */
} /* mb04ow_ */

