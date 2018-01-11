#line 1 "MB04IY.f"
/* MB04IY.f -- translated by f2c (version 20100827).
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

#line 1 "MB04IY.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb04iy_(char *side, char *trans, integer *n, integer *m, 
	integer *k, integer *p, doublereal *a, integer *lda, doublereal *tau, 
	doublereal *c__, integer *ldc, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen side_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, i__1, i__2;

    /* Local variables */
    static integer i__;
    static doublereal aii;
    static logical left, tran;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dormqr_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);
    static doublereal wrkopt;


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

/*     To overwrite the real n-by-m matrix  C  with  Q' * C,  Q * C, */
/*     C * Q',  or  C * Q,  according to the following table */

/*                     SIDE = 'L'     SIDE = 'R' */
/*     TRANS = 'N':      Q * C          C * Q */
/*     TRANS = 'T':      Q'* C          C * Q' */

/*     where  Q  is a real orthogonal matrix defined as the product of */
/*     k elementary reflectors */

/*        Q = H(1) H(2) . . . H(k) */

/*     as returned by SLICOT Library routine MB04ID.  Q  is of order n */
/*     if  SIDE = 'L'  and of order m if  SIDE = 'R'. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Specify if  Q  or  Q'  is applied from the left or right, */
/*             as follows: */
/*             = 'L':  apply  Q  or  Q'  from the left; */
/*             = 'R':  apply  Q  or  Q'  from the right. */

/*     TRANS   CHARACTER*1 */
/*             Specify if  Q  or  Q'  is to be applied, as follows: */
/*             = 'N':  apply  Q   (No transpose); */
/*             = 'T':  apply  Q'  (Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix C.  M >= 0. */

/*     K       (input) INTEGER */
/*             The number of elementary reflectors whose product defines */
/*             the matrix Q. */
/*             N >= K >= 0,  if  SIDE = 'L'; */
/*             M >= K >= 0,  if  SIDE = 'R'. */

/*     P       (input) INTEGER */
/*             The order of the zero triagle (or the number of rows of */
/*             the zero trapezoid) in the matrix triangularized by SLICOT */
/*             Library routine MB04ID.  P >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,K) */
/*             On input, the elements in the rows  i+1:min(n,n-p-1+i)  of */
/*             the  i-th  column, and  TAU(i),  represent the orthogonal */
/*             reflector  H(i),  so that matrix  Q  is the product of */
/*             elementary reflectors:  Q = H(1) H(2) . . . H(k). */
/*             A is modified by the routine but restored on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of the array  A. */
/*             LDA >= max(1,N),  if  SIDE = 'L'; */
/*             LDA >= max(1,M),  if  SIDE = 'R'. */

/*     TAU     (input) DOUBLE PRECISION array, dimension (K) */
/*             The scalar factors of the elementary reflectors. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix  C. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the updated matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array  C.  LDC >= max(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,M),  if  SIDE = 'L'; */
/*             LDWORK >= MAX(1,N),  if  SIDE = 'R'. */
/*             For optimum performance LDWORK >= M*NB if SIDE = 'L', */
/*             or LDWORK >= N*NB if SIDE = 'R', where NB is the optimal */
/*             block size. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If  SIDE = 'L',  each elementary reflector  H(i)  modifies */
/*     n-p  elements of each column of  C,  for  i = 1:p+1,  and */
/*     n-i+1  elements, for  i = p+2:k. */
/*     If  SIDE = 'R',  each elementary reflector  H(i)  modifies */
/*     m-p  elements of each row of  C,  for  i = 1:p+1,  and */
/*     m-i+1  elements, for  i = p+2:k. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Matrix operations, QR decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Check the scalar input arguments. */

#line 168 "MB04IY.f"
    /* Parameter adjustments */
#line 168 "MB04IY.f"
    a_dim1 = *lda;
#line 168 "MB04IY.f"
    a_offset = 1 + a_dim1;
#line 168 "MB04IY.f"
    a -= a_offset;
#line 168 "MB04IY.f"
    --tau;
#line 168 "MB04IY.f"
    c_dim1 = *ldc;
#line 168 "MB04IY.f"
    c_offset = 1 + c_dim1;
#line 168 "MB04IY.f"
    c__ -= c_offset;
#line 168 "MB04IY.f"
    --dwork;
#line 168 "MB04IY.f"

#line 168 "MB04IY.f"
    /* Function Body */
#line 168 "MB04IY.f"
    *info = 0;
#line 169 "MB04IY.f"
    left = lsame_(side, "L", (ftnlen)1, (ftnlen)1);
#line 170 "MB04IY.f"
    tran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);

#line 172 "MB04IY.f"
    if (! left && ! lsame_(side, "R", (ftnlen)1, (ftnlen)1)) {
#line 173 "MB04IY.f"
	*info = -1;
#line 174 "MB04IY.f"
    } else if (! tran && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 175 "MB04IY.f"
	*info = -2;
#line 176 "MB04IY.f"
    } else if (*n < 0) {
#line 177 "MB04IY.f"
	*info = -3;
#line 178 "MB04IY.f"
    } else if (*m < 0) {
#line 179 "MB04IY.f"
	*info = -4;
#line 180 "MB04IY.f"
    } else if (*k < 0 || left && *k > *n || ! left && *k > *m) {
#line 182 "MB04IY.f"
	*info = -5;
#line 183 "MB04IY.f"
    } else if (*p < 0) {
#line 184 "MB04IY.f"
	*info = -6;
#line 185 "MB04IY.f"
    } else if (left && *lda < max(1,*n) || ! left && *lda < max(1,*m)) {
#line 187 "MB04IY.f"
	*info = -8;
#line 188 "MB04IY.f"
    } else if (*ldc < max(1,*n)) {
#line 189 "MB04IY.f"
	*info = -11;
#line 190 "MB04IY.f"
    } else if (left && *ldwork < max(1,*m) || ! left && *ldwork < max(1,*n)) {
#line 192 "MB04IY.f"
	*info = -13;
#line 193 "MB04IY.f"
    }

#line 195 "MB04IY.f"
    if (*info != 0) {
#line 196 "MB04IY.f"
	i__1 = -(*info);
#line 196 "MB04IY.f"
	xerbla_("MB04IY", &i__1, (ftnlen)6);
#line 197 "MB04IY.f"
	return 0;
#line 198 "MB04IY.f"
    }

/*     Quick return if possible. */

#line 202 "MB04IY.f"
    if (*m == 0 || *n == 0 || *k == 0 || left && *n < *p || ! left && *m < *p)
	     {
#line 204 "MB04IY.f"
	dwork[1] = 1.;
#line 205 "MB04IY.f"
	return 0;
#line 206 "MB04IY.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 214 "MB04IY.f"
    if (left) {
#line 215 "MB04IY.f"
	wrkopt = (doublereal) (*m);
#line 216 "MB04IY.f"
	if (tran) {

#line 218 "MB04IY.f"
	    i__1 = min(*k,*p);
#line 218 "MB04IY.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Apply H(i) to C(i:i+n-p-1,1:m), from the left. */
/*              Workspace: need M. */

#line 223 "MB04IY.f"
		aii = a[i__ + i__ * a_dim1];
#line 224 "MB04IY.f"
		a[i__ + i__ * a_dim1] = 1.;
#line 225 "MB04IY.f"
		i__2 = *n - *p;
#line 225 "MB04IY.f"
		dlarf_(side, &i__2, m, &a[i__ + i__ * a_dim1], &c__1, &tau[
			i__], &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)1);
#line 227 "MB04IY.f"
		a[i__ + i__ * a_dim1] = aii;
#line 228 "MB04IY.f"
/* L10: */
#line 228 "MB04IY.f"
	    }

#line 230 "MB04IY.f"
	    if (*p <= min(*n,*k)) {

/*              Apply H(i) to C, i = p+1:k, from the left. */
/*              Workspace: need M;  prefer M*NB. */

#line 235 "MB04IY.f"
		i__1 = *n - *p;
#line 235 "MB04IY.f"
		i__2 = *k - *p;
#line 235 "MB04IY.f"
		dormqr_(side, trans, &i__1, m, &i__2, &a[*p + 1 + (*p + 1) * 
			a_dim1], lda, &tau[*p + 1], &c__[*p + 1 + c_dim1], 
			ldc, &dwork[1], ldwork, &i__, (ftnlen)1, (ftnlen)1);
#line 238 "MB04IY.f"
		wrkopt = max(wrkopt,dwork[1]);
#line 239 "MB04IY.f"
	    }

#line 241 "MB04IY.f"
	} else {

#line 243 "MB04IY.f"
	    if (*p <= min(*n,*k)) {

/*              Apply H(i) to C, i = k:p+1:-1, from the left. */
/*              Workspace: need M;  prefer M*NB. */

#line 248 "MB04IY.f"
		i__1 = *n - *p;
#line 248 "MB04IY.f"
		i__2 = *k - *p;
#line 248 "MB04IY.f"
		dormqr_(side, trans, &i__1, m, &i__2, &a[*p + 1 + (*p + 1) * 
			a_dim1], lda, &tau[*p + 1], &c__[*p + 1 + c_dim1], 
			ldc, &dwork[1], ldwork, &i__, (ftnlen)1, (ftnlen)1);
#line 251 "MB04IY.f"
		wrkopt = max(wrkopt,dwork[1]);
#line 252 "MB04IY.f"
	    }

#line 254 "MB04IY.f"
	    for (i__ = min(*k,*p); i__ >= 1; --i__) {

/*              Apply H(i) to C(i:i+n-p-1,1:m), from the left. */
/*              Workspace: need M. */

#line 259 "MB04IY.f"
		aii = a[i__ + i__ * a_dim1];
#line 260 "MB04IY.f"
		a[i__ + i__ * a_dim1] = 1.;
#line 261 "MB04IY.f"
		i__1 = *n - *p;
#line 261 "MB04IY.f"
		dlarf_(side, &i__1, m, &a[i__ + i__ * a_dim1], &c__1, &tau[
			i__], &c__[i__ + c_dim1], ldc, &dwork[1], (ftnlen)1);
#line 263 "MB04IY.f"
		a[i__ + i__ * a_dim1] = aii;
#line 264 "MB04IY.f"
/* L20: */
#line 264 "MB04IY.f"
	    }
#line 265 "MB04IY.f"
	}

#line 267 "MB04IY.f"
    } else {

#line 269 "MB04IY.f"
	wrkopt = (doublereal) (*n);
#line 270 "MB04IY.f"
	if (tran) {

#line 272 "MB04IY.f"
	    if (*p <= min(*m,*k)) {

/*              Apply H(i) to C, i = k:p+1:-1, from the right. */
/*              Workspace: need N;  prefer N*NB. */

#line 277 "MB04IY.f"
		i__1 = *m - *p;
#line 277 "MB04IY.f"
		i__2 = *k - *p;
#line 277 "MB04IY.f"
		dormqr_(side, trans, n, &i__1, &i__2, &a[*p + 1 + (*p + 1) * 
			a_dim1], lda, &tau[*p + 1], &c__[(*p + 1) * c_dim1 + 
			1], ldc, &dwork[1], ldwork, &i__, (ftnlen)1, (ftnlen)
			1);
#line 280 "MB04IY.f"
		wrkopt = max(wrkopt,dwork[1]);
#line 281 "MB04IY.f"
	    }

#line 283 "MB04IY.f"
	    for (i__ = min(*k,*p); i__ >= 1; --i__) {

/*              Apply H(i) to C(1:n,i:i+m-p-1), from the right. */
/*              Workspace: need N. */

#line 288 "MB04IY.f"
		aii = a[i__ + i__ * a_dim1];
#line 289 "MB04IY.f"
		a[i__ + i__ * a_dim1] = 1.;
#line 290 "MB04IY.f"
		i__1 = *m - *p;
#line 290 "MB04IY.f"
		dlarf_(side, n, &i__1, &a[i__ + i__ * a_dim1], &c__1, &tau[
			i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (ftnlen)
			1);
#line 292 "MB04IY.f"
		a[i__ + i__ * a_dim1] = aii;
#line 293 "MB04IY.f"
/* L30: */
#line 293 "MB04IY.f"
	    }

#line 295 "MB04IY.f"
	} else {

#line 297 "MB04IY.f"
	    i__1 = min(*k,*p);
#line 297 "MB04IY.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Apply H(i) to C(1:n,i:i+m-p-1), from the right. */
/*              Workspace: need N. */

#line 302 "MB04IY.f"
		aii = a[i__ + i__ * a_dim1];
#line 303 "MB04IY.f"
		a[i__ + i__ * a_dim1] = 1.;
#line 304 "MB04IY.f"
		i__2 = *m - *p;
#line 304 "MB04IY.f"
		dlarf_(side, n, &i__2, &a[i__ + i__ * a_dim1], &c__1, &tau[
			i__], &c__[i__ * c_dim1 + 1], ldc, &dwork[1], (ftnlen)
			1);
#line 306 "MB04IY.f"
		a[i__ + i__ * a_dim1] = aii;
#line 307 "MB04IY.f"
/* L40: */
#line 307 "MB04IY.f"
	    }

#line 309 "MB04IY.f"
	    if (*p <= min(*m,*k)) {

/*              Apply H(i) to C, i = p+1:k, from the right. */
/*              Workspace: need N;  prefer N*NB. */

#line 314 "MB04IY.f"
		i__1 = *m - *p;
#line 314 "MB04IY.f"
		i__2 = *k - *p;
#line 314 "MB04IY.f"
		dormqr_(side, trans, n, &i__1, &i__2, &a[*p + 1 + (*p + 1) * 
			a_dim1], lda, &tau[*p + 1], &c__[(*p + 1) * c_dim1 + 
			1], ldc, &dwork[1], ldwork, &i__, (ftnlen)1, (ftnlen)
			1);
#line 317 "MB04IY.f"
		wrkopt = max(wrkopt,dwork[1]);
#line 318 "MB04IY.f"
	    }

#line 320 "MB04IY.f"
	}
#line 321 "MB04IY.f"
    }

#line 323 "MB04IY.f"
    dwork[1] = wrkopt;
#line 324 "MB04IY.f"
    return 0;

/* *** Last line of MB04IY *** */
} /* mb04iy_ */

