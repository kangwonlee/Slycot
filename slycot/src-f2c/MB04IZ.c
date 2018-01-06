#line 1 "MB04IZ.f"
/* MB04IZ.f -- translated by f2c (version 20100827).
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

#line 1 "MB04IZ.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int mb04iz_(integer *n, integer *m, integer *p, integer *l, 
	doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, 
	doublecomplex *tau, doublecomplex *zwork, integer *lzwork, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2, i__3, i__4;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    static integer i__, nb;
    extern /* Subroutine */ int zlarf_(char *, integer *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, doublecomplex *, ftnlen);
    static doublecomplex first;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int zlarfg_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *), zgeqrf_(integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, integer *);
    static logical lquery;
    static integer wrkopt;
    extern /* Subroutine */ int zunmqr_(char *, char *, integer *, integer *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *, integer *,
	     ftnlen, ftnlen);


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

/*     To compute a QR factorization of an n-by-m matrix A (A = Q * R), */
/*     having a p-by-min(p,m) zero triangle in the lower left-hand side */
/*     corner, as shown below, for n = 8, m = 7, and p = 2: */

/*            [ x x x x x x x ] */
/*            [ x x x x x x x ] */
/*            [ x x x x x x x ] */
/*            [ x x x x x x x ] */
/*        A = [ x x x x x x x ], */
/*            [ x x x x x x x ] */
/*            [ 0 x x x x x x ] */
/*            [ 0 0 x x x x x ] */

/*     and optionally apply the transformations to an n-by-l matrix B */
/*     (from the left). The problem structure is exploited. This */
/*     computation is useful, for instance, in combined measurement and */
/*     time update of one iteration of the time-invariant Kalman filter */
/*     (square root information filter). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix A.  M >= 0. */

/*     P       (input) INTEGER */
/*             The order of the zero triagle.  P >= 0. */

/*     L       (input) INTEGER */
/*             The number of columns of the matrix B.  L >= 0. */

/*     A       (input/output) COMPLEX*16 array, dimension (LDA,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix A. The elements corresponding to the */
/*             zero P-by-MIN(P,M) lower trapezoidal/triangular part */
/*             (if P > 0) are not referenced. */
/*             On exit, the elements on and above the diagonal of this */
/*             array contain the MIN(N,M)-by-M upper trapezoidal matrix */
/*             R (R is upper triangular, if N >= M) of the QR */
/*             factorization, and the relevant elements below the */
/*             diagonal contain the trailing components (the vectors v, */
/*             see Method) of the elementary reflectors used in the */
/*             factorization. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) COMPLEX*16 array, dimension (LDB,L) */
/*             On entry, the leading N-by-L part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading N-by-L part of this array contains */
/*             the updated matrix B. */
/*             If L = 0, this array is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,N) if L > 0; */
/*             LDB >= 1        if L = 0. */

/*     TAU     (output) COMPLEX*16 array, dimension MIN(N,M) */
/*             The scalar factors of the elementary reflectors used. */

/*     Workspace */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */
/*             On exit, if INFO = 0, ZWORK(1) returns the optimal value */
/*             of LZWORK. */

/*     LZWORK  The length of the array ZWORK. */
/*             LZWORK >= MAX(1,M-1,M-P,L). */
/*             For optimum performance LZWORK should be larger. */

/*             If LZWORK = -1, then a workspace query is assumed; */
/*             the routine only calculates the optimal size of the */
/*             ZWORK array, returns this value as the first entry of */
/*             the ZWORK array, and no error message related to LZWORK */
/*             is issued by XERBLA. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine uses min(N,M) Householder transformations exploiting */
/*     the zero pattern of the matrix.  A Householder matrix has the form */

/*                                     ( 1 ), */
/*        H  = I - tau *u *u',    u  = ( v ) */
/*         i          i  i  i      i   (  i) */

/*     where v  is an (N-P+I-2)-vector.  The components of v  are stored */
/*            i                                             i */
/*     in the i-th column of A, beginning from the location i+1, and */
/*     tau  is stored in TAU(i). */
/*        i */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Nov. 2008. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2009. */

/*     KEYWORDS */

/*     Elementary reflector, QR factorization, unitary transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 171 "MB04IZ.f"
    /* Parameter adjustments */
#line 171 "MB04IZ.f"
    a_dim1 = *lda;
#line 171 "MB04IZ.f"
    a_offset = 1 + a_dim1;
#line 171 "MB04IZ.f"
    a -= a_offset;
#line 171 "MB04IZ.f"
    b_dim1 = *ldb;
#line 171 "MB04IZ.f"
    b_offset = 1 + b_dim1;
#line 171 "MB04IZ.f"
    b -= b_offset;
#line 171 "MB04IZ.f"
    --tau;
#line 171 "MB04IZ.f"
    --zwork;
#line 171 "MB04IZ.f"

#line 171 "MB04IZ.f"
    /* Function Body */
#line 171 "MB04IZ.f"
    *info = 0;
#line 172 "MB04IZ.f"
    lquery = *lzwork == -1;
#line 173 "MB04IZ.f"
    if (*n < 0) {
#line 174 "MB04IZ.f"
	*info = -1;
#line 175 "MB04IZ.f"
    } else if (*m < 0) {
#line 176 "MB04IZ.f"
	*info = -2;
#line 177 "MB04IZ.f"
    } else if (*p < 0) {
#line 178 "MB04IZ.f"
	*info = -3;
#line 179 "MB04IZ.f"
    } else if (*l < 0) {
#line 180 "MB04IZ.f"
	*info = -4;
#line 181 "MB04IZ.f"
    } else if (*lda < max(1,*n)) {
#line 182 "MB04IZ.f"
	*info = -6;
#line 183 "MB04IZ.f"
    } else if (*ldb < 1 || *l > 0 && *ldb < *n) {
#line 184 "MB04IZ.f"
	*info = -8;
#line 185 "MB04IZ.f"
    } else {
/* Computing MAX */
#line 186 "MB04IZ.f"
	i__1 = 1, i__2 = *m - 1, i__1 = max(i__1,i__2), i__2 = *m - *p, i__1 =
		 max(i__1,i__2);
#line 186 "MB04IZ.f"
	i__ = max(i__1,*l);
#line 187 "MB04IZ.f"
	if (lquery) {
#line 188 "MB04IZ.f"
	    if (*m > *p) {
#line 189 "MB04IZ.f"
		i__1 = *n - *p;
#line 189 "MB04IZ.f"
		i__2 = *m - *p;
#line 189 "MB04IZ.f"
		nb = ilaenv_(&c__1, "ZGEQRF", " ", &i__1, &i__2, &c_n1, &c_n1,
			 (ftnlen)6, (ftnlen)1);
/* Computing MAX */
#line 190 "MB04IZ.f"
		i__1 = i__, i__2 = (*m - *p) * nb;
#line 190 "MB04IZ.f"
		wrkopt = max(i__1,i__2);
#line 191 "MB04IZ.f"
		if (*l > 0) {
/* Computing MIN */
#line 192 "MB04IZ.f"
		    i__3 = *n - *p;
#line 192 "MB04IZ.f"
		    i__4 = min(*n,*m) - *p;
#line 192 "MB04IZ.f"
		    i__1 = 64, i__2 = ilaenv_(&c__1, "ZUNMQR", "LC", &i__3, l,
			     &i__4, &c_n1, (ftnlen)6, (ftnlen)2);
#line 192 "MB04IZ.f"
		    nb = min(i__1,i__2);
/* Computing MAX */
#line 194 "MB04IZ.f"
		    i__1 = wrkopt, i__2 = max(1,*l) * nb;
#line 194 "MB04IZ.f"
		    wrkopt = max(i__1,i__2);
#line 195 "MB04IZ.f"
		}
#line 196 "MB04IZ.f"
	    }
#line 197 "MB04IZ.f"
	} else if (*lzwork < i__) {
#line 198 "MB04IZ.f"
	    *info = -11;
#line 199 "MB04IZ.f"
	}
#line 200 "MB04IZ.f"
    }

#line 202 "MB04IZ.f"
    if (*info != 0) {

/*        Error return. */

#line 206 "MB04IZ.f"
	i__1 = -(*info);
#line 206 "MB04IZ.f"
	xerbla_("MB04IZ", &i__1, (ftnlen)6);
#line 207 "MB04IZ.f"
	return 0;
#line 208 "MB04IZ.f"
    } else if (lquery) {
#line 209 "MB04IZ.f"
	zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
#line 210 "MB04IZ.f"
	return 0;
#line 211 "MB04IZ.f"
    }

/*     Quick return if possible. */

#line 215 "MB04IZ.f"
    if (min(*m,*n) == 0) {
#line 216 "MB04IZ.f"
	zwork[1].r = 1., zwork[1].i = 0.;
#line 217 "MB04IZ.f"
	return 0;
#line 218 "MB04IZ.f"
    } else if (*n <= *p + 1) {
#line 219 "MB04IZ.f"
	i__1 = min(*n,*m);
#line 219 "MB04IZ.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 220 "MB04IZ.f"
	    i__2 = i__;
#line 220 "MB04IZ.f"
	    tau[i__2].r = 0., tau[i__2].i = 0.;
#line 221 "MB04IZ.f"
/* L5: */
#line 221 "MB04IZ.f"
	}
#line 222 "MB04IZ.f"
	zwork[1].r = 1., zwork[1].i = 0.;
#line 223 "MB04IZ.f"
	return 0;
#line 224 "MB04IZ.f"
    }

/*     Annihilate the subdiagonal elements of A and apply the */
/*     transformations to B, if L > 0. */
/*     Workspace: need MAX(M-1,L). */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of complex workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 236 "MB04IZ.f"
    i__1 = min(*p,*m);
#line 236 "MB04IZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Exploit the structure of the I-th column of A. */

#line 240 "MB04IZ.f"
	i__2 = *n - *p;
#line 240 "MB04IZ.f"
	zlarfg_(&i__2, &a[i__ + i__ * a_dim1], &a[i__ + 1 + i__ * a_dim1], &
		c__1, &tau[i__]);
#line 241 "MB04IZ.f"
	i__2 = i__;
#line 241 "MB04IZ.f"
	if (tau[i__2].r != 0. || tau[i__2].i != 0.) {

#line 243 "MB04IZ.f"
	    i__2 = i__ + i__ * a_dim1;
#line 243 "MB04IZ.f"
	    first.r = a[i__2].r, first.i = a[i__2].i;
#line 244 "MB04IZ.f"
	    i__2 = i__ + i__ * a_dim1;
#line 244 "MB04IZ.f"
	    a[i__2].r = 1., a[i__2].i = 0.;

#line 246 "MB04IZ.f"
	    if (i__ < *m) {
#line 246 "MB04IZ.f"
		i__2 = *n - *p;
#line 246 "MB04IZ.f"
		i__3 = *m - i__;
#line 246 "MB04IZ.f"
		d_cnjg(&z__1, &tau[i__]);
#line 246 "MB04IZ.f"
		zlarf_("Left", &i__2, &i__3, &a[i__ + i__ * a_dim1], &c__1, &
			z__1, &a[i__ + (i__ + 1) * a_dim1], lda, &zwork[1], (
			ftnlen)4);
#line 246 "MB04IZ.f"
	    }
#line 249 "MB04IZ.f"
	    if (*l > 0) {
#line 249 "MB04IZ.f"
		i__2 = *n - *p;
#line 249 "MB04IZ.f"
		d_cnjg(&z__1, &tau[i__]);
#line 249 "MB04IZ.f"
		zlarf_("Left", &i__2, l, &a[i__ + i__ * a_dim1], &c__1, &z__1,
			 &b[i__ + b_dim1], ldb, &zwork[1], (ftnlen)4);
#line 249 "MB04IZ.f"
	    }

#line 253 "MB04IZ.f"
	    i__2 = i__ + i__ * a_dim1;
#line 253 "MB04IZ.f"
	    a[i__2].r = first.r, a[i__2].i = first.i;
#line 254 "MB04IZ.f"
	}
#line 255 "MB04IZ.f"
/* L10: */
#line 255 "MB04IZ.f"
    }

/* Computing MAX */
#line 257 "MB04IZ.f"
    i__1 = 1, i__2 = *m - 1, i__1 = max(i__1,i__2);
#line 257 "MB04IZ.f"
    wrkopt = max(i__1,*l);

/*     Fast QR factorization of the remaining right submatrix, if any. */
/*     Workspace: need M-P;  prefer (M-P)*NB. */

#line 262 "MB04IZ.f"
    if (*m > *p) {
#line 263 "MB04IZ.f"
	i__1 = *n - *p;
#line 263 "MB04IZ.f"
	i__2 = *m - *p;
#line 263 "MB04IZ.f"
	zgeqrf_(&i__1, &i__2, &a[*p + 1 + (*p + 1) * a_dim1], lda, &tau[*p + 
		1], &zwork[1], lzwork, info);
/* Computing MAX */
#line 265 "MB04IZ.f"
	i__1 = wrkopt, i__2 = (integer) zwork[1].r;
#line 265 "MB04IZ.f"
	wrkopt = max(i__1,i__2);

#line 267 "MB04IZ.f"
	if (*l > 0) {

/*           Apply the transformations to B. */
/*           Workspace: need L;  prefer L*NB. */

#line 272 "MB04IZ.f"
	    i__1 = *n - *p;
#line 272 "MB04IZ.f"
	    i__2 = min(*n,*m) - *p;
#line 272 "MB04IZ.f"
	    zunmqr_("Left", "Conjugate", &i__1, l, &i__2, &a[*p + 1 + (*p + 1)
		     * a_dim1], lda, &tau[*p + 1], &b[*p + 1 + b_dim1], ldb, &
		    zwork[1], lzwork, info, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 275 "MB04IZ.f"
	    i__1 = wrkopt, i__2 = (integer) zwork[1].r;
#line 275 "MB04IZ.f"
	    wrkopt = max(i__1,i__2);
#line 276 "MB04IZ.f"
	}
#line 277 "MB04IZ.f"
    }

#line 279 "MB04IZ.f"
    zwork[1].r = (doublereal) wrkopt, zwork[1].i = 0.;
#line 280 "MB04IZ.f"
    return 0;
/* *** Last line of MB04IZ *** */
} /* mb04iz_ */

