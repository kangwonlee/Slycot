#line 1 "MB04QU.f"
/* MB04QU.f -- translated by f2c (version 20100827).
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

#line 1 "MB04QU.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb04qu_(char *tranc, char *trand, char *tranq, char *
	storev, char *storew, integer *m, integer *n, integer *k, doublereal *
	v, integer *ldv, doublereal *w, integer *ldw, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *cs, 
	doublereal *tau, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen tranc_len, ftnlen trand_len, ftnlen tranq_len, ftnlen 
	storev_len, ftnlen storew_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, d_dim1, d_offset, v_dim1, v_offset, w_dim1, 
	    w_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal nu;
    static logical ltrc, ltrd;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static logical ltrq;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lcolv, lcolw;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To overwrite general real m-by-n matrices C and D, or their */
/*     transposes, with */

/*               [ op(C) ] */
/*         Q  *  [       ]   if TRANQ = 'N', or */
/*               [ op(D) ] */

/*          T    [ op(C) ] */
/*         Q  *  [       ]   if TRANQ = 'T', */
/*               [ op(D) ] */

/*     where Q is defined as the product of symplectic reflectors and */
/*     Givens rotators, */

/*         Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) ) */
/*             diag( H(2),H(2) ) G(2) diag( F(2),F(2) ) */
/*                               .... */
/*             diag( H(k),H(k) ) G(k) diag( F(k),F(k) ). */

/*     Unblocked version. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANC   CHARACTER*1 */
/*             Specifies the form of op( C ) as follows: */
/*             = 'N':  op( C ) = C; */
/*             = 'T':  op( C ) = C'; */
/*             = 'C':  op( C ) = C'. */

/*     STOREV  CHARACTER*1 */
/*             Specifies how the vectors which define the concatenated */
/*             Householder reflectors contained in V are stored: */
/*             = 'C':  columnwise; */
/*             = 'R':  rowwise. */

/*     STOREW  CHARACTER*1 */
/*             Specifies how the vectors which define the concatenated */
/*             Householder reflectors contained in W are stored: */
/*             = 'C':  columnwise; */
/*             = 'R':  rowwise. */

/*     TRAND   CHARACTER*1 */
/*             Specifies the form of op( D ) as follows: */
/*             = 'N':  op( D ) = D; */
/*             = 'T':  op( D ) = D'; */
/*             = 'C':  op( D ) = D'. */

/*     TRANQ   CHARACTER*1 */
/*             = 'N':  apply Q; */
/*             = 'T':  apply Q'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices op(C) and op(D). */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrices op(C) and op(D). */
/*             N >= 0. */

/*     K       (input) INTEGER */
/*             The number of elementary reflectors whose product defines */
/*             the matrix Q.  M >= K >= 0. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*                     (LDV,K) if STOREV = 'C', */
/*                     (LDV,M) if STOREV = 'R' */
/*             On entry with STOREV = 'C', the leading M-by-K part of */
/*             this array must contain in its columns the vectors which */
/*             define the elementary reflectors F(i). */
/*             On entry with STOREV = 'R', the leading K-by-M part of */
/*             this array must contain in its rows the vectors which */
/*             define the elementary reflectors F(i). */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= MAX(1,M),  if STOREV = 'C'; */
/*             LDV >= MAX(1,K),  if STOREV = 'R'. */

/*     W       (input) DOUBLE PRECISION array, dimension */
/*                     (LDW,K) if STOREW = 'C', */
/*                     (LDW,M) if STOREW = 'R' */
/*             On entry with STOREW = 'C', the leading M-by-K part of */
/*             this array must contain in its columns the vectors which */
/*             define the elementary reflectors H(i). */
/*             On entry with STOREW = 'R', the leading K-by-M part of */
/*             this array must contain in its rows the vectors which */
/*             define the elementary reflectors H(i). */

/*     LDW     INTEGER */
/*             The leading dimension of the array W. */
/*             LDW >= MAX(1,M),  if STOREW = 'C'; */
/*             LDW >= MAX(1,K),  if STOREW = 'R'. */

/*     C       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDC,N) if TRANC = 'N', */
/*                     (LDC,M) if TRANC = 'T' or TRANC = 'C' */
/*             On entry with TRANC = 'N', the leading M-by-N part of */
/*             this array must contain the matrix C. */
/*             On entry with TRANC = 'C' or TRANC = 'T', the leading */
/*             N-by-M part of this array must contain the transpose of */
/*             the matrix C. */
/*             On exit with TRANC = 'N', the leading M-by-N part of */
/*             this array contains the updated matrix C. */
/*             On exit with TRANC = 'C' or TRANC = 'T', the leading */
/*             N-by-M part of this array contains the transpose of the */
/*             updated matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= MAX(1,M),  if TRANC = 'N'; */
/*             LDC >= MAX(1,N),  if TRANC = 'T' or TRANC = 'C'. */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDD,N) if TRAND = 'N', */
/*                     (LDD,M) if TRAND = 'T' or TRAND = 'C' */
/*             On entry with TRAND = 'N', the leading M-by-N part of */
/*             this array must contain the matrix D. */
/*             On entry with TRAND = 'C' or TRAND = 'T', the leading */
/*             N-by-M part of this array must contain the transpose of */
/*             the matrix D. */
/*             On exit with TRAND = 'N', the leading M-by-N part of */
/*             this array contains the updated matrix D. */
/*             On exit with TRAND = 'C' or TRAND = 'T', the leading */
/*             N-by-M part of this array contains the transpose of the */
/*             updated matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D. */
/*             LDD >= MAX(1,M),  if TRAND = 'N'; */
/*             LDD >= MAX(1,N),  if TRAND = 'T' or TRAND = 'C'. */

/*     CS      (input) DOUBLE PRECISION array, dimension (2*K) */
/*             On entry, the first 2*K elements of this array must */
/*             contain the cosines and sines of the symplectic Givens */
/*             rotators G(i). */

/*     TAU     (input) DOUBLE PRECISION array, dimension (K) */
/*             On entry, the first K elements of this array must */
/*             contain the scalar factors of the elementary reflectors */
/*             F(i). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -20,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSMSQ). */

/*     KEYWORDS */

/*     Elementary matrix operations, orthogonal symplectic matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 228 "MB04QU.f"
    /* Parameter adjustments */
#line 228 "MB04QU.f"
    v_dim1 = *ldv;
#line 228 "MB04QU.f"
    v_offset = 1 + v_dim1;
#line 228 "MB04QU.f"
    v -= v_offset;
#line 228 "MB04QU.f"
    w_dim1 = *ldw;
#line 228 "MB04QU.f"
    w_offset = 1 + w_dim1;
#line 228 "MB04QU.f"
    w -= w_offset;
#line 228 "MB04QU.f"
    c_dim1 = *ldc;
#line 228 "MB04QU.f"
    c_offset = 1 + c_dim1;
#line 228 "MB04QU.f"
    c__ -= c_offset;
#line 228 "MB04QU.f"
    d_dim1 = *ldd;
#line 228 "MB04QU.f"
    d_offset = 1 + d_dim1;
#line 228 "MB04QU.f"
    d__ -= d_offset;
#line 228 "MB04QU.f"
    --cs;
#line 228 "MB04QU.f"
    --tau;
#line 228 "MB04QU.f"
    --dwork;
#line 228 "MB04QU.f"

#line 228 "MB04QU.f"
    /* Function Body */
#line 228 "MB04QU.f"
    *info = 0;
#line 229 "MB04QU.f"
    lcolv = lsame_(storev, "C", (ftnlen)1, (ftnlen)1);
#line 230 "MB04QU.f"
    lcolw = lsame_(storew, "C", (ftnlen)1, (ftnlen)1);
#line 231 "MB04QU.f"
    ltrc = lsame_(tranc, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranc, "C", (
	    ftnlen)1, (ftnlen)1);
#line 232 "MB04QU.f"
    ltrd = lsame_(trand, "T", (ftnlen)1, (ftnlen)1) || lsame_(trand, "C", (
	    ftnlen)1, (ftnlen)1);
#line 233 "MB04QU.f"
    ltrq = lsame_(tranq, "T", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 237 "MB04QU.f"
    if (! (ltrc || lsame_(tranc, "N", (ftnlen)1, (ftnlen)1))) {
#line 238 "MB04QU.f"
	*info = -1;
#line 239 "MB04QU.f"
    } else if (! (ltrd || lsame_(trand, "N", (ftnlen)1, (ftnlen)1))) {
#line 240 "MB04QU.f"
	*info = -2;
#line 241 "MB04QU.f"
    } else if (! (ltrq || lsame_(tranq, "N", (ftnlen)1, (ftnlen)1))) {
#line 242 "MB04QU.f"
	*info = -3;
#line 243 "MB04QU.f"
    } else if (! (lcolv || lsame_(storev, "R", (ftnlen)1, (ftnlen)1))) {
#line 244 "MB04QU.f"
	*info = -4;
#line 245 "MB04QU.f"
    } else if (! (lcolw || lsame_(storew, "R", (ftnlen)1, (ftnlen)1))) {
#line 246 "MB04QU.f"
	*info = -5;
#line 247 "MB04QU.f"
    } else if (*m < 0) {
#line 248 "MB04QU.f"
	*info = -6;
#line 249 "MB04QU.f"
    } else if (*n < 0) {
#line 250 "MB04QU.f"
	*info = -7;
#line 251 "MB04QU.f"
    } else if (*k < 0 || *k > *m) {
#line 252 "MB04QU.f"
	*info = -8;
#line 253 "MB04QU.f"
    } else if (lcolv && *ldv < max(1,*m) || ! lcolv && *ldv < max(1,*k)) {
#line 255 "MB04QU.f"
	*info = -10;
#line 256 "MB04QU.f"
    } else if (lcolw && *ldw < max(1,*m) || ! lcolw && *ldw < max(1,*k)) {
#line 258 "MB04QU.f"
	*info = -12;
#line 259 "MB04QU.f"
    } else if (ltrc && *ldc < max(1,*n) || ! ltrc && *ldc < max(1,*m)) {
#line 261 "MB04QU.f"
	*info = -14;
#line 262 "MB04QU.f"
    } else if (ltrd && *ldd < max(1,*n) || ! ltrd && *ldd < max(1,*m)) {
#line 264 "MB04QU.f"
	*info = -16;
#line 265 "MB04QU.f"
    } else if (*ldwork < max(1,*n)) {
#line 266 "MB04QU.f"
	dwork[1] = (doublereal) max(1,*n);
#line 267 "MB04QU.f"
	*info = -20;
#line 268 "MB04QU.f"
    }

/*     Return if there were illegal values. */

#line 272 "MB04QU.f"
    if (*info != 0) {
#line 273 "MB04QU.f"
	i__1 = -(*info);
#line 273 "MB04QU.f"
	xerbla_("MB04QU", &i__1, (ftnlen)6);
#line 274 "MB04QU.f"
	return 0;
#line 275 "MB04QU.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 279 "MB04QU.f"
    i__1 = min(*k,*m);
#line 279 "MB04QU.f"
    if (min(i__1,*n) == 0) {
#line 280 "MB04QU.f"
	dwork[1] = 1.;
#line 281 "MB04QU.f"
	return 0;
#line 282 "MB04QU.f"
    }

#line 284 "MB04QU.f"
    if (ltrq) {
#line 285 "MB04QU.f"
	i__1 = *k;
#line 285 "MB04QU.f"
	for (i__ = 1; i__ <= i__1; ++i__) {

/*           Apply H(I) to C(I:M,:) and D(I:M,:) from the left. */

#line 289 "MB04QU.f"
	    nu = w[i__ + i__ * w_dim1];
#line 290 "MB04QU.f"
	    w[i__ + i__ * w_dim1] = 1.;
#line 291 "MB04QU.f"
	    if (lcolw) {
#line 292 "MB04QU.f"
		if (ltrc) {
#line 293 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 293 "MB04QU.f"
		    dlarf_("Right", n, &i__2, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (
			    ftnlen)5);
#line 295 "MB04QU.f"
		} else {
#line 296 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 296 "MB04QU.f"
		    dlarf_("Left", &i__2, n, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)4)
			    ;
#line 298 "MB04QU.f"
		}
#line 299 "MB04QU.f"
		if (ltrd) {
#line 300 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 300 "MB04QU.f"
		    dlarf_("Right", n, &i__2, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &d__[i__ * d_dim1 + 1], ldd, &dwork[1], (
			    ftnlen)5);
#line 302 "MB04QU.f"
		} else {
#line 303 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 303 "MB04QU.f"
		    dlarf_("Left", &i__2, n, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &d__[i__ + d_dim1], ldd, &dwork[1], (ftnlen)4)
			    ;
#line 305 "MB04QU.f"
		}
#line 306 "MB04QU.f"
	    } else {
#line 307 "MB04QU.f"
		if (ltrc) {
#line 308 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 308 "MB04QU.f"
		    dlarf_("Right", n, &i__2, &w[i__ + i__ * w_dim1], ldw, &
			    nu, &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (
			    ftnlen)5);
#line 310 "MB04QU.f"
		} else {
#line 311 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 311 "MB04QU.f"
		    dlarf_("Left", &i__2, n, &w[i__ + i__ * w_dim1], ldw, &nu,
			     &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)4);
#line 313 "MB04QU.f"
		}
#line 314 "MB04QU.f"
		if (ltrd) {
#line 315 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 315 "MB04QU.f"
		    dlarf_("Right", n, &i__2, &w[i__ + i__ * w_dim1], ldw, &
			    nu, &d__[i__ * d_dim1 + 1], ldd, &dwork[1], (
			    ftnlen)5);
#line 317 "MB04QU.f"
		} else {
#line 318 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 318 "MB04QU.f"
		    dlarf_("Left", &i__2, n, &w[i__ + i__ * w_dim1], ldw, &nu,
			     &d__[i__ + d_dim1], ldd, &dwork[1], (ftnlen)4);
#line 320 "MB04QU.f"
		}
#line 321 "MB04QU.f"
	    }
#line 322 "MB04QU.f"
	    w[i__ + i__ * w_dim1] = nu;

/*           Apply G(i) to C(I,:) and D(I,:) from the left. */

#line 326 "MB04QU.f"
	    if (ltrc && ltrd) {
#line 327 "MB04QU.f"
		drot_(n, &c__[i__ * c_dim1 + 1], &c__1, &d__[i__ * d_dim1 + 1]
			, &c__1, &cs[(i__ << 1) - 1], &cs[i__ * 2]);
#line 328 "MB04QU.f"
	    } else if (ltrc) {
#line 329 "MB04QU.f"
		drot_(n, &c__[i__ * c_dim1 + 1], &c__1, &d__[i__ + d_dim1], 
			ldd, &cs[(i__ << 1) - 1], &cs[i__ * 2]);
#line 331 "MB04QU.f"
	    } else if (ltrd) {
#line 332 "MB04QU.f"
		drot_(n, &c__[i__ + c_dim1], ldc, &d__[i__ * d_dim1 + 1], &
			c__1, &cs[(i__ << 1) - 1], &cs[i__ * 2]);
#line 334 "MB04QU.f"
	    } else {
#line 335 "MB04QU.f"
		drot_(n, &c__[i__ + c_dim1], ldc, &d__[i__ + d_dim1], ldd, &
			cs[(i__ << 1) - 1], &cs[i__ * 2]);
#line 337 "MB04QU.f"
	    }

/*           Apply F(I) to C(I:M,:) and D(I:M,:) from the left. */

#line 341 "MB04QU.f"
	    nu = v[i__ + i__ * v_dim1];
#line 342 "MB04QU.f"
	    v[i__ + i__ * v_dim1] = 1.;
#line 343 "MB04QU.f"
	    if (lcolv) {
#line 344 "MB04QU.f"
		if (ltrc) {
#line 345 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 345 "MB04QU.f"
		    dlarf_("Right", n, &i__2, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], 
			    (ftnlen)5);
#line 347 "MB04QU.f"
		} else {
#line 348 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 348 "MB04QU.f"
		    dlarf_("Left", &i__2, n, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &c__[i__ + c_dim1], ldc, &dwork[1], (
			    ftnlen)4);
#line 350 "MB04QU.f"
		}
#line 351 "MB04QU.f"
		if (ltrd) {
#line 352 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 352 "MB04QU.f"
		    dlarf_("Right", n, &i__2, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &d__[i__ * d_dim1 + 1], ldd, &dwork[1], 
			    (ftnlen)5);
#line 354 "MB04QU.f"
		} else {
#line 355 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 355 "MB04QU.f"
		    dlarf_("Left", &i__2, n, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &d__[i__ + d_dim1], ldd, &dwork[1], (
			    ftnlen)4);
#line 357 "MB04QU.f"
		}
#line 358 "MB04QU.f"
	    } else {
#line 359 "MB04QU.f"
		if (ltrc) {
#line 360 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 360 "MB04QU.f"
		    dlarf_("Right", n, &i__2, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], 
			    (ftnlen)5);
#line 362 "MB04QU.f"
		} else {
#line 363 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 363 "MB04QU.f"
		    dlarf_("Left", &i__2, n, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &c__[i__ + c_dim1], ldc, &dwork[1], (
			    ftnlen)4);
#line 365 "MB04QU.f"
		}
#line 366 "MB04QU.f"
		if (ltrd) {
#line 367 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 367 "MB04QU.f"
		    dlarf_("Right", n, &i__2, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &d__[i__ * d_dim1 + 1], ldd, &dwork[1], 
			    (ftnlen)5);
#line 369 "MB04QU.f"
		} else {
#line 370 "MB04QU.f"
		    i__2 = *m - i__ + 1;
#line 370 "MB04QU.f"
		    dlarf_("Left", &i__2, n, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &d__[i__ + d_dim1], ldd, &dwork[1], (
			    ftnlen)4);
#line 372 "MB04QU.f"
		}
#line 373 "MB04QU.f"
	    }
#line 374 "MB04QU.f"
	    v[i__ + i__ * v_dim1] = nu;
#line 375 "MB04QU.f"
/* L10: */
#line 375 "MB04QU.f"
	}
#line 376 "MB04QU.f"
    } else {
#line 377 "MB04QU.f"
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to C(I:M,:) and D(I:M,:) from the left. */

#line 381 "MB04QU.f"
	    nu = v[i__ + i__ * v_dim1];
#line 382 "MB04QU.f"
	    v[i__ + i__ * v_dim1] = 1.;
#line 383 "MB04QU.f"
	    if (lcolv) {
#line 384 "MB04QU.f"
		if (ltrc) {
#line 385 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 385 "MB04QU.f"
		    dlarf_("Right", n, &i__1, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], 
			    (ftnlen)5);
#line 387 "MB04QU.f"
		} else {
#line 388 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 388 "MB04QU.f"
		    dlarf_("Left", &i__1, n, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &c__[i__ + c_dim1], ldc, &dwork[1], (
			    ftnlen)4);
#line 390 "MB04QU.f"
		}
#line 391 "MB04QU.f"
		if (ltrd) {
#line 392 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 392 "MB04QU.f"
		    dlarf_("Right", n, &i__1, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &d__[i__ * d_dim1 + 1], ldd, &dwork[1], 
			    (ftnlen)5);
#line 394 "MB04QU.f"
		} else {
#line 395 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 395 "MB04QU.f"
		    dlarf_("Left", &i__1, n, &v[i__ + i__ * v_dim1], &c__1, &
			    tau[i__], &d__[i__ + d_dim1], ldd, &dwork[1], (
			    ftnlen)4);
#line 397 "MB04QU.f"
		}
#line 398 "MB04QU.f"
	    } else {
#line 399 "MB04QU.f"
		if (ltrc) {
#line 400 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 400 "MB04QU.f"
		    dlarf_("Right", n, &i__1, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], 
			    (ftnlen)5);
#line 402 "MB04QU.f"
		} else {
#line 403 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 403 "MB04QU.f"
		    dlarf_("Left", &i__1, n, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &c__[i__ + c_dim1], ldc, &dwork[1], (
			    ftnlen)4);
#line 405 "MB04QU.f"
		}
#line 406 "MB04QU.f"
		if (ltrd) {
#line 407 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 407 "MB04QU.f"
		    dlarf_("Right", n, &i__1, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &d__[i__ * d_dim1 + 1], ldd, &dwork[1], 
			    (ftnlen)5);
#line 409 "MB04QU.f"
		} else {
#line 410 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 410 "MB04QU.f"
		    dlarf_("Left", &i__1, n, &v[i__ + i__ * v_dim1], ldv, &
			    tau[i__], &d__[i__ + d_dim1], ldd, &dwork[1], (
			    ftnlen)4);
#line 412 "MB04QU.f"
		}
#line 413 "MB04QU.f"
	    }
#line 414 "MB04QU.f"
	    v[i__ + i__ * v_dim1] = nu;

/*           Apply G(i) to C(I,:) and D(I,:) from the left. */

#line 418 "MB04QU.f"
	    if (ltrc && ltrd) {
#line 419 "MB04QU.f"
		d__1 = -cs[i__ * 2];
#line 419 "MB04QU.f"
		drot_(n, &c__[i__ * c_dim1 + 1], &c__1, &d__[i__ * d_dim1 + 1]
			, &c__1, &cs[(i__ << 1) - 1], &d__1);
#line 420 "MB04QU.f"
	    } else if (ltrc) {
#line 421 "MB04QU.f"
		d__1 = -cs[i__ * 2];
#line 421 "MB04QU.f"
		drot_(n, &c__[i__ * c_dim1 + 1], &c__1, &d__[i__ + d_dim1], 
			ldd, &cs[(i__ << 1) - 1], &d__1);
#line 423 "MB04QU.f"
	    } else if (ltrd) {
#line 424 "MB04QU.f"
		d__1 = -cs[i__ * 2];
#line 424 "MB04QU.f"
		drot_(n, &c__[i__ + c_dim1], ldc, &d__[i__ * d_dim1 + 1], &
			c__1, &cs[(i__ << 1) - 1], &d__1);
#line 426 "MB04QU.f"
	    } else {
#line 427 "MB04QU.f"
		d__1 = -cs[i__ * 2];
#line 427 "MB04QU.f"
		drot_(n, &c__[i__ + c_dim1], ldc, &d__[i__ + d_dim1], ldd, &
			cs[(i__ << 1) - 1], &d__1);
#line 429 "MB04QU.f"
	    }

/*           Apply H(I) to C(I:M,:) and D(I:M,:) from the left. */

#line 433 "MB04QU.f"
	    nu = w[i__ + i__ * w_dim1];
#line 434 "MB04QU.f"
	    w[i__ + i__ * w_dim1] = 1.;
#line 435 "MB04QU.f"
	    if (lcolw) {
#line 436 "MB04QU.f"
		if (ltrc) {
#line 437 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 437 "MB04QU.f"
		    dlarf_("Right", n, &i__1, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (
			    ftnlen)5);
#line 439 "MB04QU.f"
		} else {
#line 440 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 440 "MB04QU.f"
		    dlarf_("Left", &i__1, n, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)4)
			    ;
#line 442 "MB04QU.f"
		}
#line 443 "MB04QU.f"
		if (ltrd) {
#line 444 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 444 "MB04QU.f"
		    dlarf_("Right", n, &i__1, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &d__[i__ * d_dim1 + 1], ldd, &dwork[1], (
			    ftnlen)5);
#line 446 "MB04QU.f"
		} else {
#line 447 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 447 "MB04QU.f"
		    dlarf_("Left", &i__1, n, &w[i__ + i__ * w_dim1], &c__1, &
			    nu, &d__[i__ + d_dim1], ldd, &dwork[1], (ftnlen)4)
			    ;
#line 449 "MB04QU.f"
		}
#line 450 "MB04QU.f"
	    } else {
#line 451 "MB04QU.f"
		if (ltrc) {
#line 452 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 452 "MB04QU.f"
		    dlarf_("Right", n, &i__1, &w[i__ + i__ * w_dim1], ldw, &
			    nu, &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (
			    ftnlen)5);
#line 454 "MB04QU.f"
		} else {
#line 455 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 455 "MB04QU.f"
		    dlarf_("Left", &i__1, n, &w[i__ + i__ * w_dim1], ldw, &nu,
			     &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)4);
#line 457 "MB04QU.f"
		}
#line 458 "MB04QU.f"
		if (ltrd) {
#line 459 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 459 "MB04QU.f"
		    dlarf_("Right", n, &i__1, &w[i__ + i__ * w_dim1], ldw, &
			    nu, &d__[i__ * d_dim1 + 1], ldd, &dwork[1], (
			    ftnlen)5);
#line 461 "MB04QU.f"
		} else {
#line 462 "MB04QU.f"
		    i__1 = *m - i__ + 1;
#line 462 "MB04QU.f"
		    dlarf_("Left", &i__1, n, &w[i__ + i__ * w_dim1], ldw, &nu,
			     &d__[i__ + d_dim1], ldd, &dwork[1], (ftnlen)4);
#line 464 "MB04QU.f"
		}
#line 465 "MB04QU.f"
	    }
#line 466 "MB04QU.f"
	    w[i__ + i__ * w_dim1] = nu;
#line 467 "MB04QU.f"
/* L20: */
#line 467 "MB04QU.f"
	}
#line 468 "MB04QU.f"
    }

#line 470 "MB04QU.f"
    dwork[1] = (doublereal) max(1,*n);
/* *** Last line of MB04QU *** */
#line 472 "MB04QU.f"
    return 0;
} /* mb04qu_ */

