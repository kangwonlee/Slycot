#line 1 "MB04WU.f"
/* MB04WU.f -- translated by f2c (version 20100827).
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

#line 1 "MB04WU.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static integer c__1 = 1;

/* Subroutine */ int mb04wu_(char *tranq1, char *tranq2, integer *m, integer *
	n, integer *k, doublereal *q1, integer *ldq1, doublereal *q2, integer 
	*ldq2, doublereal *cs, doublereal *tau, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen tranq1_len, ftnlen tranq2_len)
{
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, j;
    static doublereal nu;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static logical ltrq1, ltrq2;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlarf_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

/*     To generate a matrix Q with orthogonal columns (spanning an */
/*     isotropic subspace), which is defined as the first n columns */
/*     of a product of symplectic reflectors and Givens rotators, */

/*         Q = diag( H(1),H(1) ) G(1) diag( F(1),F(1) ) */
/*             diag( H(2),H(2) ) G(2) diag( F(2),F(2) ) */
/*                               .... */
/*             diag( H(k),H(k) ) G(k) diag( F(k),F(k) ). */

/*     The matrix Q is returned in terms of its first 2*M rows */

/*                      [  op( Q1 )   op( Q2 ) ] */
/*                  Q = [                      ]. */
/*                      [ -op( Q2 )   op( Q1 ) ] */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANQ1  CHARACTER*1 */
/*             Specifies the form of op( Q1 ) as follows: */
/*             = 'N':  op( Q1 ) = Q1; */
/*             = 'T':  op( Q1 ) = Q1'; */
/*             = 'C':  op( Q1 ) = Q1'. */

/*     TRANQ2  CHARACTER*1 */
/*             Specifies the form of op( Q2 ) as follows: */
/*             = 'N':  op( Q2 ) = Q2; */
/*             = 'T':  op( Q2 ) = Q2'; */
/*             = 'C':  op( Q2 ) = Q2'. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices Q1 and Q2. M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrices Q1 and Q2. */
/*             M >= N >= 0. */

/*     K       (input) INTEGER */
/*             The number of symplectic Givens rotators whose product */
/*             partly defines the matrix Q. N >= K >= 0. */

/*     Q1      (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDQ1,N) if TRANQ1 = 'N', */
/*                     (LDQ1,M) if TRANQ1 = 'T' or TRANQ1 = 'C' */
/*             On entry with TRANQ1 = 'N', the leading M-by-K part of */
/*             this array must contain in its i-th column the vector */
/*             which defines the elementary reflector F(i). */
/*             On entry with TRANQ1 = 'T' or TRANQ1 = 'C', the leading */
/*             K-by-M part of this array must contain in its i-th row */
/*             the vector which defines the elementary reflector F(i). */
/*             On exit with TRANQ1 = 'N', the leading M-by-N part of this */
/*             array contains the matrix Q1. */
/*             On exit with TRANQ1 = 'T' or TRANQ1 = 'C', the leading */
/*             N-by-M part of this array contains the matrix Q1'. */

/*     LDQ1    INTEGER */
/*             The leading dimension of the array Q1. */
/*             LDQ1 >= MAX(1,M),  if TRANQ1 = 'N'; */
/*             LDQ1 >= MAX(1,N),  if TRANQ1 = 'T' or TRANQ1 = 'C'. */

/*     Q2      (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDQ2,N) if TRANQ2 = 'N', */
/*                     (LDQ2,M) if TRANQ2 = 'T' or TRANQ2 = 'C' */
/*             On entry with TRANQ2 = 'N', the leading M-by-K part of */
/*             this array must contain in its i-th column the vector */
/*             which defines the elementary reflector H(i) and, on the */
/*             diagonal, the scalar factor of H(i). */
/*             On entry with TRANQ2 = 'T' or TRANQ2 = 'C', the leading */
/*             K-by-M part of this array must contain in its i-th row the */
/*             vector which defines the elementary reflector H(i) and, on */
/*             the diagonal, the scalar factor of H(i). */
/*             On exit with TRANQ2 = 'N', the leading M-by-N part of this */
/*             array contains the matrix Q2. */
/*             On exit with TRANQ2 = 'T' or TRANQ2 = 'C', the leading */
/*             N-by-M part of this array contains the matrix Q2'. */

/*     LDQ2    INTEGER */
/*             The leading dimension of the array Q2. */
/*             LDQ2 >= MAX(1,M),  if TRANQ2 = 'N'; */
/*             LDQ2 >= MAX(1,N),  if TRANQ2 = 'T' or TRANQ2 = 'C'. */

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
/*             On exit, if  INFO = -13,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,M+N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     REFERENCES */

/*     [1] Bunse-Gerstner, A. */
/*         Matrix factorizations for symplectic QR-like methods. */
/*         Linear Algebra Appl., 83, pp. 49-77, 1986. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSGSQ). */

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

#line 181 "MB04WU.f"
    /* Parameter adjustments */
#line 181 "MB04WU.f"
    q1_dim1 = *ldq1;
#line 181 "MB04WU.f"
    q1_offset = 1 + q1_dim1;
#line 181 "MB04WU.f"
    q1 -= q1_offset;
#line 181 "MB04WU.f"
    q2_dim1 = *ldq2;
#line 181 "MB04WU.f"
    q2_offset = 1 + q2_dim1;
#line 181 "MB04WU.f"
    q2 -= q2_offset;
#line 181 "MB04WU.f"
    --cs;
#line 181 "MB04WU.f"
    --tau;
#line 181 "MB04WU.f"
    --dwork;
#line 181 "MB04WU.f"

#line 181 "MB04WU.f"
    /* Function Body */
#line 181 "MB04WU.f"
    *info = 0;
#line 182 "MB04WU.f"
    ltrq1 = lsame_(tranq1, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranq1, "C", (
	    ftnlen)1, (ftnlen)1);
#line 183 "MB04WU.f"
    ltrq2 = lsame_(tranq2, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranq2, "C", (
	    ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 187 "MB04WU.f"
    if (! (ltrq1 || lsame_(tranq1, "N", (ftnlen)1, (ftnlen)1))) {
#line 188 "MB04WU.f"
	*info = -1;
#line 189 "MB04WU.f"
    } else if (! (ltrq2 || lsame_(tranq2, "N", (ftnlen)1, (ftnlen)1))) {
#line 190 "MB04WU.f"
	*info = -2;
#line 191 "MB04WU.f"
    } else if (*m < 0) {
#line 192 "MB04WU.f"
	*info = -3;
#line 193 "MB04WU.f"
    } else if (*n < 0 || *n > *m) {
#line 194 "MB04WU.f"
	*info = -4;
#line 195 "MB04WU.f"
    } else if (*k < 0 || *k > *n) {
#line 196 "MB04WU.f"
	*info = -5;
#line 197 "MB04WU.f"
    } else if (ltrq1 && *ldq1 < max(1,*n) || ! ltrq1 && *ldq1 < max(1,*m)) {
#line 199 "MB04WU.f"
	*info = -7;
#line 200 "MB04WU.f"
    } else if (ltrq2 && *ldq2 < max(1,*n) || ! ltrq2 && *ldq2 < max(1,*m)) {
#line 202 "MB04WU.f"
	*info = -9;
#line 203 "MB04WU.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 203 "MB04WU.f"
	i__1 = 1, i__2 = *m + *n;
#line 203 "MB04WU.f"
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
#line 204 "MB04WU.f"
	    i__1 = 1, i__2 = *m + *n;
#line 204 "MB04WU.f"
	    dwork[1] = (doublereal) max(i__1,i__2);
#line 205 "MB04WU.f"
	    *info = -13;
#line 206 "MB04WU.f"
	}
#line 206 "MB04WU.f"
    }

/*     Return if there were illegal values. */

#line 210 "MB04WU.f"
    if (*info != 0) {
#line 211 "MB04WU.f"
	i__1 = -(*info);
#line 211 "MB04WU.f"
	xerbla_("MB04WU", &i__1, (ftnlen)6);
#line 212 "MB04WU.f"
	return 0;
#line 213 "MB04WU.f"
    }

/*     Quick return if possible. */

#line 217 "MB04WU.f"
    if (*n == 0) {
#line 218 "MB04WU.f"
	dwork[1] = 1.;
#line 219 "MB04WU.f"
	return 0;
#line 220 "MB04WU.f"
    }

/*     Initialize columns K+1:N to columns of the unit matrix. */

#line 224 "MB04WU.f"
    i__1 = *n;
#line 224 "MB04WU.f"
    for (j = *k + 1; j <= i__1; ++j) {
#line 225 "MB04WU.f"
	i__2 = *m;
#line 225 "MB04WU.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 226 "MB04WU.f"
	    q1[i__ + j * q1_dim1] = 0.;
#line 227 "MB04WU.f"
/* L10: */
#line 227 "MB04WU.f"
	}
#line 228 "MB04WU.f"
	q1[j + j * q1_dim1] = 1.;
#line 229 "MB04WU.f"
/* L20: */
#line 229 "MB04WU.f"
    }
#line 230 "MB04WU.f"
    i__1 = *n - *k;
#line 230 "MB04WU.f"
    dlaset_("All", m, &i__1, &c_b12, &c_b12, &q2[(*k + 1) * q2_dim1 + 1], 
	    ldq2, (ftnlen)3);

#line 232 "MB04WU.f"
    if (ltrq1 && ltrq2) {
#line 233 "MB04WU.f"
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to Q1(I+1:N,I:M) and Q2(I+1:N,I:M) from the */
/*           right. */

#line 238 "MB04WU.f"
	    i__1 = *m - i__ + 1;
#line 238 "MB04WU.f"
	    dcopy_(&i__1, &q2[i__ + i__ * q2_dim1], ldq2, &dwork[1], &c__1);
#line 239 "MB04WU.f"
	    if (i__ < *n) {
#line 240 "MB04WU.f"
		q1[i__ + i__ * q1_dim1] = 1.;
#line 241 "MB04WU.f"
		i__1 = *n - i__;
#line 241 "MB04WU.f"
		i__2 = *m - i__ + 1;
#line 241 "MB04WU.f"
		dlarf_("Right", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], ldq1, 
			&tau[i__], &q1[i__ + 1 + i__ * q1_dim1], ldq1, &dwork[
			*m + 1], (ftnlen)5);
#line 243 "MB04WU.f"
		i__1 = *n - i__;
#line 243 "MB04WU.f"
		i__2 = *m - i__ + 1;
#line 243 "MB04WU.f"
		dlarf_("Right", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], ldq1, 
			&tau[i__], &q2[i__ + 1 + i__ * q2_dim1], ldq2, &dwork[
			*m + 1], (ftnlen)5);
#line 245 "MB04WU.f"
	    }
#line 246 "MB04WU.f"
	    if (i__ < *m) {
#line 246 "MB04WU.f"
		i__1 = *m - i__;
#line 246 "MB04WU.f"
		d__1 = -tau[i__];
#line 246 "MB04WU.f"
		dscal_(&i__1, &d__1, &q1[i__ + (i__ + 1) * q1_dim1], ldq1);
#line 246 "MB04WU.f"
	    }
#line 248 "MB04WU.f"
	    q1[i__ + i__ * q1_dim1] = 1. - tau[i__];

/*           Set Q1(I,1:I-1) and Q2(I,1:M) to zero. */

#line 252 "MB04WU.f"
	    i__1 = i__ - 1;
#line 252 "MB04WU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 253 "MB04WU.f"
		q1[i__ + j * q1_dim1] = 0.;
#line 254 "MB04WU.f"
/* L30: */
#line 254 "MB04WU.f"
	    }
#line 255 "MB04WU.f"
	    i__1 = *m;
#line 255 "MB04WU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 256 "MB04WU.f"
		q2[i__ + j * q2_dim1] = 0.;
#line 257 "MB04WU.f"
/* L40: */
#line 257 "MB04WU.f"
	    }

/*           Apply G(I) to Q1(I:N,I) and Q2(I:N,I) from the right. */

#line 261 "MB04WU.f"
	    i__1 = *n - i__ + 1;
#line 261 "MB04WU.f"
	    drot_(&i__1, &q1[i__ + i__ * q1_dim1], &c__1, &q2[i__ + i__ * 
		    q2_dim1], &c__1, &cs[(i__ << 1) - 1], &cs[i__ * 2]);

/*           Apply H(I) to Q1(I:N,I:M) and Q2(I:N,I:M) from the right. */

#line 266 "MB04WU.f"
	    nu = dwork[1];
#line 267 "MB04WU.f"
	    dwork[1] = 1.;
#line 268 "MB04WU.f"
	    i__1 = *n - i__ + 1;
#line 268 "MB04WU.f"
	    i__2 = *m - i__ + 1;
#line 268 "MB04WU.f"
	    dlarf_("Right", &i__1, &i__2, &dwork[1], &c__1, &nu, &q1[i__ + 
		    i__ * q1_dim1], ldq1, &dwork[*m + 1], (ftnlen)5);
#line 270 "MB04WU.f"
	    i__1 = *n - i__ + 1;
#line 270 "MB04WU.f"
	    i__2 = *m - i__ + 1;
#line 270 "MB04WU.f"
	    dlarf_("Right", &i__1, &i__2, &dwork[1], &c__1, &nu, &q2[i__ + 
		    i__ * q2_dim1], ldq2, &dwork[*m + 1], (ftnlen)5);
#line 272 "MB04WU.f"
/* L50: */
#line 272 "MB04WU.f"
	}
#line 273 "MB04WU.f"
    } else if (ltrq1) {
#line 274 "MB04WU.f"
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to Q1(I+1:N,I:M) from the right and to */
/*           Q2(I:M,I+1:N) from the left. */

#line 279 "MB04WU.f"
	    i__1 = *m - i__ + 1;
#line 279 "MB04WU.f"
	    dcopy_(&i__1, &q2[i__ + i__ * q2_dim1], &c__1, &dwork[1], &c__1);
#line 280 "MB04WU.f"
	    if (i__ < *n) {
#line 281 "MB04WU.f"
		q1[i__ + i__ * q1_dim1] = 1.;
#line 282 "MB04WU.f"
		i__1 = *n - i__;
#line 282 "MB04WU.f"
		i__2 = *m - i__ + 1;
#line 282 "MB04WU.f"
		dlarf_("Right", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], ldq1, 
			&tau[i__], &q1[i__ + 1 + i__ * q1_dim1], ldq1, &dwork[
			*m + 1], (ftnlen)5);
#line 284 "MB04WU.f"
		i__1 = *m - i__ + 1;
#line 284 "MB04WU.f"
		i__2 = *n - i__;
#line 284 "MB04WU.f"
		dlarf_("Left", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], ldq1, &
			tau[i__], &q2[i__ + (i__ + 1) * q2_dim1], ldq2, &
			dwork[*m + 1], (ftnlen)4);
#line 286 "MB04WU.f"
	    }
#line 287 "MB04WU.f"
	    if (i__ < *m) {
#line 287 "MB04WU.f"
		i__1 = *m - i__;
#line 287 "MB04WU.f"
		d__1 = -tau[i__];
#line 287 "MB04WU.f"
		dscal_(&i__1, &d__1, &q1[i__ + (i__ + 1) * q1_dim1], ldq1);
#line 287 "MB04WU.f"
	    }
#line 289 "MB04WU.f"
	    q1[i__ + i__ * q1_dim1] = 1. - tau[i__];

/*           Set Q1(I,1:I-1) and Q2(1:M,I) to zero. */

#line 293 "MB04WU.f"
	    i__1 = i__ - 1;
#line 293 "MB04WU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 294 "MB04WU.f"
		q1[i__ + j * q1_dim1] = 0.;
#line 295 "MB04WU.f"
/* L60: */
#line 295 "MB04WU.f"
	    }
#line 296 "MB04WU.f"
	    i__1 = *m;
#line 296 "MB04WU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 297 "MB04WU.f"
		q2[j + i__ * q2_dim1] = 0.;
#line 298 "MB04WU.f"
/* L70: */
#line 298 "MB04WU.f"
	    }

/*           Apply G(I) to Q1(I:N,I) from the right and to Q2(I,I:N) */
/*           from the left. */

#line 303 "MB04WU.f"
	    i__1 = *n - i__ + 1;
#line 303 "MB04WU.f"
	    drot_(&i__1, &q1[i__ + i__ * q1_dim1], &c__1, &q2[i__ + i__ * 
		    q2_dim1], ldq2, &cs[(i__ << 1) - 1], &cs[i__ * 2]);

/*           Apply H(I) to Q1(I:N,I:M) from the right and to Q2(I:M,I:N) */
/*           from the left. */

#line 309 "MB04WU.f"
	    nu = dwork[1];
#line 310 "MB04WU.f"
	    dwork[1] = 1.;
#line 311 "MB04WU.f"
	    i__1 = *n - i__ + 1;
#line 311 "MB04WU.f"
	    i__2 = *m - i__ + 1;
#line 311 "MB04WU.f"
	    dlarf_("Right", &i__1, &i__2, &dwork[1], &c__1, &nu, &q1[i__ + 
		    i__ * q1_dim1], ldq1, &dwork[*m + 1], (ftnlen)5);
#line 313 "MB04WU.f"
	    i__1 = *m - i__ + 1;
#line 313 "MB04WU.f"
	    i__2 = *n - i__ + 1;
#line 313 "MB04WU.f"
	    dlarf_("Left", &i__1, &i__2, &dwork[1], &c__1, &nu, &q2[i__ + i__ 
		    * q2_dim1], ldq2, &dwork[*m + 1], (ftnlen)4);
#line 315 "MB04WU.f"
/* L80: */
#line 315 "MB04WU.f"
	}
#line 316 "MB04WU.f"
    } else if (ltrq2) {
#line 317 "MB04WU.f"
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to Q1(I:M,I+1:N) from the left and to */
/*           Q2(I+1:N,I:M) from the right. */

#line 322 "MB04WU.f"
	    i__1 = *m - i__ + 1;
#line 322 "MB04WU.f"
	    dcopy_(&i__1, &q2[i__ + i__ * q2_dim1], ldq2, &dwork[1], &c__1);
#line 323 "MB04WU.f"
	    if (i__ < *n) {
#line 324 "MB04WU.f"
		q1[i__ + i__ * q1_dim1] = 1.;
#line 325 "MB04WU.f"
		i__1 = *m - i__ + 1;
#line 325 "MB04WU.f"
		i__2 = *n - i__;
#line 325 "MB04WU.f"
		dlarf_("Left", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], &c__1, 
			&tau[i__], &q1[i__ + (i__ + 1) * q1_dim1], ldq1, &
			dwork[*m + 1], (ftnlen)4);
#line 327 "MB04WU.f"
		i__1 = *n - i__;
#line 327 "MB04WU.f"
		i__2 = *m - i__ + 1;
#line 327 "MB04WU.f"
		dlarf_("Right", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], &c__1,
			 &tau[i__], &q2[i__ + 1 + i__ * q2_dim1], ldq2, &
			dwork[*m + 1], (ftnlen)5);
#line 329 "MB04WU.f"
	    }
#line 330 "MB04WU.f"
	    if (i__ < *m) {
#line 330 "MB04WU.f"
		i__1 = *m - i__;
#line 330 "MB04WU.f"
		d__1 = -tau[i__];
#line 330 "MB04WU.f"
		dscal_(&i__1, &d__1, &q1[i__ + 1 + i__ * q1_dim1], &c__1);
#line 330 "MB04WU.f"
	    }
#line 332 "MB04WU.f"
	    q1[i__ + i__ * q1_dim1] = 1. - tau[i__];

/*           Set Q1(1:I-1,I) and Q2(I,1:M) to zero. */

#line 336 "MB04WU.f"
	    i__1 = i__ - 1;
#line 336 "MB04WU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 337 "MB04WU.f"
		q1[j + i__ * q1_dim1] = 0.;
#line 338 "MB04WU.f"
/* L90: */
#line 338 "MB04WU.f"
	    }
#line 339 "MB04WU.f"
	    i__1 = *m;
#line 339 "MB04WU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 340 "MB04WU.f"
		q2[i__ + j * q2_dim1] = 0.;
#line 341 "MB04WU.f"
/* L100: */
#line 341 "MB04WU.f"
	    }

/*           Apply G(I) to Q1(I,I:N) from the left and to Q2(I:N,I) */
/*           from the right. */

#line 346 "MB04WU.f"
	    i__1 = *n - i__ + 1;
#line 346 "MB04WU.f"
	    drot_(&i__1, &q1[i__ + i__ * q1_dim1], ldq1, &q2[i__ + i__ * 
		    q2_dim1], &c__1, &cs[(i__ << 1) - 1], &cs[i__ * 2]);

/*           Apply H(I) to Q1(I:M,I:N) from the left and to Q2(I:N,I:M) */
/*           from the left. */

#line 352 "MB04WU.f"
	    nu = dwork[1];
#line 353 "MB04WU.f"
	    dwork[1] = 1.;
#line 354 "MB04WU.f"
	    i__1 = *m - i__ + 1;
#line 354 "MB04WU.f"
	    i__2 = *n - i__ + 1;
#line 354 "MB04WU.f"
	    dlarf_("Left", &i__1, &i__2, &dwork[1], &c__1, &nu, &q1[i__ + i__ 
		    * q1_dim1], ldq1, &dwork[*m + 1], (ftnlen)4);
#line 356 "MB04WU.f"
	    i__1 = *n - i__ + 1;
#line 356 "MB04WU.f"
	    i__2 = *m - i__ + 1;
#line 356 "MB04WU.f"
	    dlarf_("Right", &i__1, &i__2, &dwork[1], &c__1, &nu, &q2[i__ + 
		    i__ * q2_dim1], ldq2, &dwork[*m + 1], (ftnlen)5);
#line 358 "MB04WU.f"
/* L110: */
#line 358 "MB04WU.f"
	}
#line 359 "MB04WU.f"
    } else {
#line 360 "MB04WU.f"
	for (i__ = *k; i__ >= 1; --i__) {

/*           Apply F(I) to Q1(I:M,I+1:N) and Q2(I:M,I+1:N) from the left. */

#line 364 "MB04WU.f"
	    i__1 = *m - i__ + 1;
#line 364 "MB04WU.f"
	    dcopy_(&i__1, &q2[i__ + i__ * q2_dim1], &c__1, &dwork[1], &c__1);
#line 365 "MB04WU.f"
	    if (i__ < *n) {
#line 366 "MB04WU.f"
		q1[i__ + i__ * q1_dim1] = 1.;
#line 367 "MB04WU.f"
		i__1 = *m - i__ + 1;
#line 367 "MB04WU.f"
		i__2 = *n - i__;
#line 367 "MB04WU.f"
		dlarf_("Left", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], &c__1, 
			&tau[i__], &q1[i__ + (i__ + 1) * q1_dim1], ldq1, &
			dwork[*m + 1], (ftnlen)4);
#line 369 "MB04WU.f"
		i__1 = *m - i__ + 1;
#line 369 "MB04WU.f"
		i__2 = *n - i__;
#line 369 "MB04WU.f"
		dlarf_("Left", &i__1, &i__2, &q1[i__ + i__ * q1_dim1], &c__1, 
			&tau[i__], &q2[i__ + (i__ + 1) * q2_dim1], ldq2, &
			dwork[*m + 1], (ftnlen)4);
#line 371 "MB04WU.f"
	    }
#line 372 "MB04WU.f"
	    if (i__ < *m) {
#line 372 "MB04WU.f"
		i__1 = *m - i__;
#line 372 "MB04WU.f"
		d__1 = -tau[i__];
#line 372 "MB04WU.f"
		dscal_(&i__1, &d__1, &q1[i__ + 1 + i__ * q1_dim1], &c__1);
#line 372 "MB04WU.f"
	    }
#line 374 "MB04WU.f"
	    q1[i__ + i__ * q1_dim1] = 1. - tau[i__];

/*           Set Q1(1:I-1,I) and Q2(1:M,I) to zero. */

#line 378 "MB04WU.f"
	    i__1 = i__ - 1;
#line 378 "MB04WU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 379 "MB04WU.f"
		q1[j + i__ * q1_dim1] = 0.;
#line 380 "MB04WU.f"
/* L120: */
#line 380 "MB04WU.f"
	    }
#line 381 "MB04WU.f"
	    i__1 = *m;
#line 381 "MB04WU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 382 "MB04WU.f"
		q2[j + i__ * q2_dim1] = 0.;
#line 383 "MB04WU.f"
/* L130: */
#line 383 "MB04WU.f"
	    }

/*           Apply G(I) to Q1(I,I:N) and Q2(I,I:N) from the left. */

#line 387 "MB04WU.f"
	    i__1 = *n - i__ + 1;
#line 387 "MB04WU.f"
	    drot_(&i__1, &q1[i__ + i__ * q1_dim1], ldq1, &q2[i__ + i__ * 
		    q2_dim1], ldq2, &cs[(i__ << 1) - 1], &cs[i__ * 2]);

/*           Apply H(I) to Q1(I:M,I:N) and Q2(I:M,I:N) from the left. */

#line 392 "MB04WU.f"
	    nu = dwork[1];
#line 393 "MB04WU.f"
	    dwork[1] = 1.;
#line 394 "MB04WU.f"
	    i__1 = *m - i__ + 1;
#line 394 "MB04WU.f"
	    i__2 = *n - i__ + 1;
#line 394 "MB04WU.f"
	    dlarf_("Left", &i__1, &i__2, &dwork[1], &c__1, &nu, &q1[i__ + i__ 
		    * q1_dim1], ldq1, &dwork[*m + 1], (ftnlen)4);
#line 396 "MB04WU.f"
	    i__1 = *m - i__ + 1;
#line 396 "MB04WU.f"
	    i__2 = *n - i__ + 1;
#line 396 "MB04WU.f"
	    dlarf_("Left", &i__1, &i__2, &dwork[1], &c__1, &nu, &q2[i__ + i__ 
		    * q2_dim1], ldq2, &dwork[*m + 1], (ftnlen)4);
#line 398 "MB04WU.f"
/* L140: */
#line 398 "MB04WU.f"
	}
#line 399 "MB04WU.f"
    }
/* Computing MAX */
#line 400 "MB04WU.f"
    i__1 = 1, i__2 = *m + *n;
#line 400 "MB04WU.f"
    dwork[1] = (doublereal) max(i__1,i__2);
/* *** Last line of MB04WU *** */
#line 402 "MB04WU.f"
    return 0;
} /* mb04wu_ */

