#line 1 "MB04QB.f"
/* MB04QB.f -- translated by f2c (version 20100827).
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

#line 1 "MB04QB.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__2 = 2;

/* Subroutine */ int mb04qb_(char *tranc, char *trand, char *tranq, char *
	storev, char *storew, integer *m, integer *n, integer *k, doublereal *
	v, integer *ldv, doublereal *w, integer *ldw, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *cs, 
	doublereal *tau, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen tranc_len, ftnlen trand_len, ftnlen tranq_len, ftnlen 
	storev_len, ftnlen storew_len)
{
    /* System generated locals */
    address a__1[3];
    integer c_dim1, c_offset, d_dim1, d_offset, v_dim1, v_offset, w_dim1, 
	    w_offset, i__1, i__2[3], i__3, i__4;
    char ch__1[3];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ib, ic, id, jc, jd, nb, ki, kk, nx, pdt, pdw, ierr;
    static logical ltrc, ltrd;
    static integer pdrs;
    static logical ltrq;
    extern /* Subroutine */ int mb04qc_(char *, char *, char *, char *, char *
	    , char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, 
	    ftnlen, ftnlen), mb04qf_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    extern integer ue01md_(integer *, char *, char *, integer *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nbmin;
    extern /* Subroutine */ int mb04qu_(char *, char *, char *, char *, char *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    static logical lcolv, lcolw;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer wrkopt;


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

/*     Blocked version. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANC   CHARACTER*1 */
/*             Specifies the form of op( C ) as follows: */
/*             = 'N':  op( C ) = C; */
/*             = 'T':  op( C ) = C'; */
/*             = 'C':  op( C ) = C'. */

/*     TRAND   CHARACTER*1 */
/*             Specifies the form of op( D ) as follows: */
/*             = 'N':  op( D ) = D; */
/*             = 'T':  op( D ) = D'; */
/*             = 'C':  op( D ) = D'. */

/*     TRANQ   CHARACTER*1 */
/*             = 'N':  apply Q; */
/*             = 'T':  apply Q'. */

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

/*     REFERENCES */

/*     [1] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSMSB). */

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

#line 235 "MB04QB.f"
    /* Parameter adjustments */
#line 235 "MB04QB.f"
    v_dim1 = *ldv;
#line 235 "MB04QB.f"
    v_offset = 1 + v_dim1;
#line 235 "MB04QB.f"
    v -= v_offset;
#line 235 "MB04QB.f"
    w_dim1 = *ldw;
#line 235 "MB04QB.f"
    w_offset = 1 + w_dim1;
#line 235 "MB04QB.f"
    w -= w_offset;
#line 235 "MB04QB.f"
    c_dim1 = *ldc;
#line 235 "MB04QB.f"
    c_offset = 1 + c_dim1;
#line 235 "MB04QB.f"
    c__ -= c_offset;
#line 235 "MB04QB.f"
    d_dim1 = *ldd;
#line 235 "MB04QB.f"
    d_offset = 1 + d_dim1;
#line 235 "MB04QB.f"
    d__ -= d_offset;
#line 235 "MB04QB.f"
    --cs;
#line 235 "MB04QB.f"
    --tau;
#line 235 "MB04QB.f"
    --dwork;
#line 235 "MB04QB.f"

#line 235 "MB04QB.f"
    /* Function Body */
#line 235 "MB04QB.f"
    *info = 0;
#line 236 "MB04QB.f"
    lcolv = lsame_(storev, "C", (ftnlen)1, (ftnlen)1);
#line 237 "MB04QB.f"
    lcolw = lsame_(storew, "C", (ftnlen)1, (ftnlen)1);
#line 238 "MB04QB.f"
    ltrc = lsame_(tranc, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranc, "C", (
	    ftnlen)1, (ftnlen)1);
#line 239 "MB04QB.f"
    ltrd = lsame_(trand, "T", (ftnlen)1, (ftnlen)1) || lsame_(trand, "C", (
	    ftnlen)1, (ftnlen)1);
#line 240 "MB04QB.f"
    ltrq = lsame_(tranq, "T", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 244 "MB04QB.f"
    if (! (ltrc || lsame_(tranc, "N", (ftnlen)1, (ftnlen)1))) {
#line 245 "MB04QB.f"
	*info = -1;
#line 246 "MB04QB.f"
    } else if (! (ltrd || lsame_(trand, "N", (ftnlen)1, (ftnlen)1))) {
#line 247 "MB04QB.f"
	*info = -2;
#line 248 "MB04QB.f"
    } else if (! (ltrq || lsame_(tranq, "N", (ftnlen)1, (ftnlen)1))) {
#line 249 "MB04QB.f"
	*info = -3;
#line 250 "MB04QB.f"
    } else if (! (lcolv || lsame_(storev, "R", (ftnlen)1, (ftnlen)1))) {
#line 251 "MB04QB.f"
	*info = -4;
#line 252 "MB04QB.f"
    } else if (! (lcolw || lsame_(storew, "R", (ftnlen)1, (ftnlen)1))) {
#line 253 "MB04QB.f"
	*info = -5;
#line 254 "MB04QB.f"
    } else if (*m < 0) {
#line 255 "MB04QB.f"
	*info = -6;
#line 256 "MB04QB.f"
    } else if (*n < 0) {
#line 257 "MB04QB.f"
	*info = -7;
#line 258 "MB04QB.f"
    } else if (*k < 0 || *k > *m) {
#line 259 "MB04QB.f"
	*info = -8;
#line 260 "MB04QB.f"
    } else if (lcolv && *ldv < max(1,*m) || ! lcolv && *ldv < max(1,*k)) {
#line 262 "MB04QB.f"
	*info = -10;
#line 263 "MB04QB.f"
    } else if (lcolw && *ldw < max(1,*m) || ! lcolw && *ldw < max(1,*k)) {
#line 265 "MB04QB.f"
	*info = -12;
#line 266 "MB04QB.f"
    } else if (ltrc && *ldc < max(1,*n) || ! ltrc && *ldc < max(1,*m)) {
#line 268 "MB04QB.f"
	*info = -14;
#line 269 "MB04QB.f"
    } else if (ltrd && *ldd < max(1,*n) || ! ltrd && *ldd < max(1,*m)) {
#line 271 "MB04QB.f"
	*info = -16;
#line 272 "MB04QB.f"
    } else if (*ldwork < max(1,*n)) {
#line 273 "MB04QB.f"
	dwork[1] = (doublereal) max(1,*n);
#line 274 "MB04QB.f"
	*info = -20;
#line 275 "MB04QB.f"
    }

/*     Return if there were illegal values. */

#line 279 "MB04QB.f"
    if (*info != 0) {
#line 280 "MB04QB.f"
	i__1 = -(*info);
#line 280 "MB04QB.f"
	xerbla_("MB04QB", &i__1, (ftnlen)6);
#line 281 "MB04QB.f"
	return 0;
#line 282 "MB04QB.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 286 "MB04QB.f"
    i__1 = min(*k,*m);
#line 286 "MB04QB.f"
    if (min(i__1,*n) == 0) {
#line 287 "MB04QB.f"
	dwork[1] = 1.;
#line 288 "MB04QB.f"
	return 0;
#line 289 "MB04QB.f"
    }

#line 291 "MB04QB.f"
    nbmin = 2;
#line 292 "MB04QB.f"
    nx = 0;
#line 293 "MB04QB.f"
    wrkopt = *n;
/* Writing concatenation */
#line 294 "MB04QB.f"
    i__2[0] = 1, a__1[0] = tranc;
#line 294 "MB04QB.f"
    i__2[1] = 1, a__1[1] = trand;
#line 294 "MB04QB.f"
    i__2[2] = 1, a__1[2] = tranq;
#line 294 "MB04QB.f"
    s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)3);
#line 294 "MB04QB.f"
    nb = ue01md_(&c__1, "MB04QB", ch__1, m, n, k, (ftnlen)6, (ftnlen)3);
#line 295 "MB04QB.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
/* Writing concatenation */
#line 299 "MB04QB.f"
	i__2[0] = 1, a__1[0] = tranc;
#line 299 "MB04QB.f"
	i__2[1] = 1, a__1[1] = trand;
#line 299 "MB04QB.f"
	i__2[2] = 1, a__1[2] = tranq;
#line 299 "MB04QB.f"
	s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)3);
#line 299 "MB04QB.f"
	i__1 = 0, i__3 = ue01md_(&c__3, "MB04QB", ch__1, m, n, k, (ftnlen)6, (
		ftnlen)3);
#line 299 "MB04QB.f"
	nx = max(i__1,i__3);
#line 301 "MB04QB.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

/* Computing MAX */
#line 305 "MB04QB.f"
	    i__1 = wrkopt, i__3 = *n * 9 * nb + nb * 15 * nb;
#line 305 "MB04QB.f"
	    wrkopt = max(i__1,i__3);
#line 306 "MB04QB.f"
	    if (*ldwork < wrkopt) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 311 "MB04QB.f"
		nb = (integer) ((sqrt((doublereal) (*n * 81 * *n + *ldwork * 
			60)) - (doublereal) (*n * 9)) / 30.);
/* Computing MAX */
/* Writing concatenation */
#line 313 "MB04QB.f"
		i__2[0] = 1, a__1[0] = tranc;
#line 313 "MB04QB.f"
		i__2[1] = 1, a__1[1] = trand;
#line 313 "MB04QB.f"
		i__2[2] = 1, a__1[2] = tranq;
#line 313 "MB04QB.f"
		s_cat(ch__1, a__1, i__2, &c__3, (ftnlen)3);
#line 313 "MB04QB.f"
		i__1 = 2, i__3 = ue01md_(&c__2, "MB04QB", ch__1, m, n, k, (
			ftnlen)6, (ftnlen)3);
#line 313 "MB04QB.f"
		nbmin = max(i__1,i__3);
#line 315 "MB04QB.f"
	    }
#line 316 "MB04QB.f"
	}
#line 317 "MB04QB.f"
    }

#line 319 "MB04QB.f"
    pdrs = 1;
#line 320 "MB04QB.f"
    pdt = pdrs + nb * 6 * nb;
#line 321 "MB04QB.f"
    pdw = pdt + nb * 9 * nb;
#line 322 "MB04QB.f"
    ic = 1;
#line 323 "MB04QB.f"
    jc = 1;
#line 324 "MB04QB.f"
    id = 1;
#line 325 "MB04QB.f"
    jd = 1;

#line 327 "MB04QB.f"
    if (ltrq) {

/*        Use blocked code initially. */

#line 331 "MB04QB.f"
	if (nb >= nbmin && nb < *k && nx < *k) {
#line 332 "MB04QB.f"
	    i__1 = *k - nx;
#line 332 "MB04QB.f"
	    i__3 = nb;
#line 332 "MB04QB.f"
	    for (i__ = 1; i__3 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__3) {
/* Computing MIN */
#line 333 "MB04QB.f"
		i__4 = *k - i__ + 1;
#line 333 "MB04QB.f"
		ib = min(i__4,nb);

/*              Form the triangular factors of the symplectic block */
/*              reflector SH. */

#line 338 "MB04QB.f"
		i__4 = *m - i__ + 1;
#line 338 "MB04QB.f"
		mb04qf_("Forward", storev, storew, &i__4, &ib, &v[i__ + i__ * 
			v_dim1], ldv, &w[i__ + i__ * w_dim1], ldw, &cs[(i__ <<
			 1) - 1], &tau[i__], &dwork[pdrs], &nb, &dwork[pdt], &
			nb, &dwork[pdw], (ftnlen)7, (ftnlen)1, (ftnlen)1);

/*              Apply SH' to [ op(C)(i:m,:); op(D)(i:m,:) ] from the */
/*              left. */

#line 346 "MB04QB.f"
		if (ltrc) {
#line 347 "MB04QB.f"
		    jc = i__;
#line 348 "MB04QB.f"
		} else {
#line 349 "MB04QB.f"
		    ic = i__;
#line 350 "MB04QB.f"
		}
#line 351 "MB04QB.f"
		if (ltrd) {
#line 352 "MB04QB.f"
		    jd = i__;
#line 353 "MB04QB.f"
		} else {
#line 354 "MB04QB.f"
		    id = i__;
#line 355 "MB04QB.f"
		}
#line 356 "MB04QB.f"
		i__4 = *m - i__ + 1;
#line 356 "MB04QB.f"
		mb04qc_("No Structure", tranc, trand, tranq, "Forward", 
			storev, storew, &i__4, n, &ib, &v[i__ + i__ * v_dim1],
			 ldv, &w[i__ + i__ * w_dim1], ldw, &dwork[pdrs], &nb, 
			&dwork[pdt], &nb, &c__[ic + jc * c_dim1], ldc, &d__[
			id + jd * d_dim1], ldd, &dwork[pdw], (ftnlen)12, (
			ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)7, (ftnlen)1, 
			(ftnlen)1);
#line 361 "MB04QB.f"
/* L10: */
#line 361 "MB04QB.f"
	    }
#line 362 "MB04QB.f"
	} else {
#line 363 "MB04QB.f"
	    i__ = 1;
#line 364 "MB04QB.f"
	}

/*        Use unblocked code to update last or only block. */

#line 368 "MB04QB.f"
	if (i__ <= *k) {
#line 369 "MB04QB.f"
	    if (ltrc) {
#line 370 "MB04QB.f"
		jc = i__;
#line 371 "MB04QB.f"
	    } else {
#line 372 "MB04QB.f"
		ic = i__;
#line 373 "MB04QB.f"
	    }
#line 374 "MB04QB.f"
	    if (ltrd) {
#line 375 "MB04QB.f"
		jd = i__;
#line 376 "MB04QB.f"
	    } else {
#line 377 "MB04QB.f"
		id = i__;
#line 378 "MB04QB.f"
	    }
#line 379 "MB04QB.f"
	    i__3 = *m - i__ + 1;
#line 379 "MB04QB.f"
	    i__1 = *k - i__ + 1;
#line 379 "MB04QB.f"
	    mb04qu_(tranc, trand, tranq, storev, storew, &i__3, n, &i__1, &v[
		    i__ + i__ * v_dim1], ldv, &w[i__ + i__ * w_dim1], ldw, &
		    c__[ic + jc * c_dim1], ldc, &d__[id + jd * d_dim1], ldd, &
		    cs[(i__ << 1) - 1], &tau[i__], &dwork[1], ldwork, &ierr, (
		    ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 383 "MB04QB.f"
	}
#line 384 "MB04QB.f"
    } else {
#line 385 "MB04QB.f"
	if (nb >= nbmin && nb < *k && nx < *k) {

/*           Use blocked code after the last block. */
/*           The first kk columns are handled by the block method. */

#line 390 "MB04QB.f"
	    ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
#line 391 "MB04QB.f"
	    i__3 = *k, i__1 = ki + nb;
#line 391 "MB04QB.f"
	    kk = min(i__3,i__1);
#line 392 "MB04QB.f"
	} else {
#line 393 "MB04QB.f"
	    kk = 0;
#line 394 "MB04QB.f"
	}

/*        Use unblocked code for the last or only block. */

#line 398 "MB04QB.f"
	if (kk < *k) {
#line 399 "MB04QB.f"
	    if (ltrc) {
#line 400 "MB04QB.f"
		jc = kk + 1;
#line 401 "MB04QB.f"
	    } else {
#line 402 "MB04QB.f"
		ic = kk + 1;
#line 403 "MB04QB.f"
	    }
#line 404 "MB04QB.f"
	    if (ltrd) {
#line 405 "MB04QB.f"
		jd = kk + 1;
#line 406 "MB04QB.f"
	    } else {
#line 407 "MB04QB.f"
		id = kk + 1;
#line 408 "MB04QB.f"
	    }
#line 409 "MB04QB.f"
	    i__3 = *m - kk;
#line 409 "MB04QB.f"
	    i__1 = *k - kk;
#line 409 "MB04QB.f"
	    mb04qu_(tranc, trand, tranq, storev, storew, &i__3, n, &i__1, &v[
		    kk + 1 + (kk + 1) * v_dim1], ldv, &w[kk + 1 + (kk + 1) * 
		    w_dim1], ldw, &c__[ic + jc * c_dim1], ldc, &d__[id + jd * 
		    d_dim1], ldd, &cs[(kk << 1) + 1], &tau[kk + 1], &dwork[1],
		     ldwork, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)
		    1, (ftnlen)1);
#line 413 "MB04QB.f"
	}

/*        Blocked code. */

#line 417 "MB04QB.f"
	if (kk > 0) {
#line 418 "MB04QB.f"
	    i__3 = -nb;
#line 418 "MB04QB.f"
	    for (i__ = ki + 1; i__3 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__3) {
/* Computing MIN */
#line 419 "MB04QB.f"
		i__1 = nb, i__4 = *k - i__ + 1;
#line 419 "MB04QB.f"
		ib = min(i__1,i__4);

/*              Form the triangular factors of the symplectic block */
/*              reflector SH. */

#line 424 "MB04QB.f"
		i__1 = *m - i__ + 1;
#line 424 "MB04QB.f"
		mb04qf_("Forward", storev, storew, &i__1, &ib, &v[i__ + i__ * 
			v_dim1], ldv, &w[i__ + i__ * w_dim1], ldw, &cs[(i__ <<
			 1) - 1], &tau[i__], &dwork[pdrs], &nb, &dwork[pdt], &
			nb, &dwork[pdw], (ftnlen)7, (ftnlen)1, (ftnlen)1);

/*              Apply SH to [ op(C)(i:m,:); op(D)(i:m,:) ] from */
/*              the left. */

#line 432 "MB04QB.f"
		if (ltrc) {
#line 433 "MB04QB.f"
		    jc = i__;
#line 434 "MB04QB.f"
		} else {
#line 435 "MB04QB.f"
		    ic = i__;
#line 436 "MB04QB.f"
		}
#line 437 "MB04QB.f"
		if (ltrd) {
#line 438 "MB04QB.f"
		    jd = i__;
#line 439 "MB04QB.f"
		} else {
#line 440 "MB04QB.f"
		    id = i__;
#line 441 "MB04QB.f"
		}
#line 442 "MB04QB.f"
		i__1 = *m - i__ + 1;
#line 442 "MB04QB.f"
		mb04qc_("No Structure", tranc, trand, tranq, "Forward", 
			storev, storew, &i__1, n, &ib, &v[i__ + i__ * v_dim1],
			 ldv, &w[i__ + i__ * w_dim1], ldw, &dwork[pdrs], &nb, 
			&dwork[pdt], &nb, &c__[ic + jc * c_dim1], ldc, &d__[
			id + jd * d_dim1], ldd, &dwork[pdw], (ftnlen)12, (
			ftnlen)1, (ftnlen)1, (ftnlen)1, (ftnlen)7, (ftnlen)1, 
			(ftnlen)1);
#line 447 "MB04QB.f"
/* L20: */
#line 447 "MB04QB.f"
	    }
#line 448 "MB04QB.f"
	}
#line 449 "MB04QB.f"
    }
#line 450 "MB04QB.f"
    dwork[1] = (doublereal) wrkopt;

#line 452 "MB04QB.f"
    return 0;
/* *** Last line of MB04QB *** */
} /* mb04qb_ */

