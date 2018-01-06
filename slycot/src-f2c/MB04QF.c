#line 1 "MB04QF.f"
/* MB04QF.f -- translated by f2c (version 20100827).
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

#line 1 "MB04QF.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b11 = 0.;
static doublereal c_b49 = 1.;

/* Subroutine */ int mb04qf_(char *direct, char *storev, char *storew, 
	integer *n, integer *k, doublereal *v, integer *ldv, doublereal *w, 
	integer *ldw, doublereal *cs, doublereal *tau, doublereal *rs, 
	integer *ldrs, doublereal *t, integer *ldt, doublereal *dwork, ftnlen 
	direct_len, ftnlen storev_len, ftnlen storew_len)
{
    /* System generated locals */
    integer rs_dim1, rs_offset, t_dim1, t_offset, v_dim1, v_offset, w_dim1, 
	    w_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k2;
    static doublereal cm1;
    static integer pr1, pr2, pr3, ps1, ps2, ps3, pt11, pt12, pt13, pt21, pt22,
	     pt23, pt31, pt32, pt33;
    static doublereal vii, wii, taui;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical lcolv, lcolw;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dtrmv_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);


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

/*     To form the triangular block factors R, S and T of a symplectic */
/*     block reflector SH, which is defined as a product of 2k */
/*     concatenated Householder reflectors and k Givens rotators, */

/*         SH = diag( H(1),H(1) ) G(1) diag( F(1),F(1) ) */
/*              diag( H(2),H(2) ) G(2) diag( F(2),F(2) ) */
/*                                .... */
/*              diag( H(k),H(k) ) G(k) diag( F(k),F(k) ). */

/*     The upper triangular blocks of the matrices */

/*                                 [ S1 ]       [ T11 T12 T13 ] */
/*         R  = [ R1 R2 R3 ],  S = [ S2 ],  T = [ T21 T22 T23 ], */
/*                                 [ S3 ]       [ T31 T32 T33 ] */

/*     with R2 unit and S1, R3, T21, T31, T32 strictly upper triangular, */
/*     are stored rowwise in the arrays RS and T, respectively. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DIRECT  CHARACTER*1 */
/*             This is a dummy argument, which is reserved for future */
/*             extensions of this subroutine. Not referenced. */

/*     STOREV  CHARACTER*1 */
/*             Specifies how the vectors which define the concatenated */
/*             Householder F(i) reflectors are stored: */
/*             = 'C':  columnwise; */
/*             = 'R':  rowwise. */

/*     STOREW  CHARACTER*1 */
/*             Specifies how the vectors which define the concatenated */
/*             Householder H(i) reflectors are stored: */
/*             = 'C':  columnwise; */
/*             = 'R':  rowwise. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the Householder reflectors F(i) and H(i). */
/*             N >= 0. */

/*     K       (input) INTEGER */
/*             The number of Givens rotators.  K >= 1. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*                     (LDV,K) if STOREV = 'C', */
/*                     (LDV,N) if STOREV = 'R' */
/*             On entry with STOREV = 'C', the leading N-by-K part of */
/*             this array must contain in its i-th column the vector */
/*             which defines the elementary reflector F(i). */
/*             On entry with STOREV = 'R', the leading K-by-N part of */
/*             this array must contain in its i-th row the vector */
/*             which defines the elementary reflector F(i). */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= MAX(1,N),  if STOREV = 'C'; */
/*             LDV >= K,         if STOREV = 'R'. */

/*     W       (input) DOUBLE PRECISION array, dimension */
/*                     (LDW,K) if STOREW = 'C', */
/*                     (LDW,N) if STOREW = 'R' */
/*             On entry with STOREW = 'C', the leading N-by-K part of */
/*             this array must contain in its i-th column the vector */
/*             which defines the elementary reflector H(i). */
/*             On entry with STOREV = 'R', the leading K-by-N part of */
/*             this array must contain in its i-th row the vector */
/*             which defines the elementary reflector H(i). */

/*     LDW     INTEGER */
/*             The leading dimension of the array W. */
/*             LDW >= MAX(1,N),  if STOREW = 'C'; */
/*             LDW >= K,         if STOREW = 'R'. */

/*     CS      (input) DOUBLE PRECISION array, dimension (2*K) */
/*             On entry, the first 2*K elements of this array must */
/*             contain the cosines and sines of the symplectic Givens */
/*             rotators G(i). */

/*     TAU     (input) DOUBLE PRECISION array, dimension (K) */
/*             On entry, the first K elements of this array must */
/*             contain the scalar factors of the elementary reflectors */
/*             F(i). */

/*     RS      (output) DOUBLE PRECISION array, dimension (K,6*K) */
/*             On exit, the leading K-by-6*K part of this array contains */
/*             the upper triangular matrices defining the factors R and */
/*             S of the symplectic block reflector SH. The (strictly) */
/*             lower portions of this array are not used. */

/*     LDRS    INTEGER */
/*             The leading dimension of the array RS.  LDRS >= K. */

/*     T       (output) DOUBLE PRECISION array, dimension (K,9*K) */
/*             On exit, the leading K-by-9*K part of this array contains */
/*             the upper triangular matrices defining the factor T of the */
/*             symplectic block reflector SH. The (strictly) lower */
/*             portions of this array are not used. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= K. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (3*K) */

/*     REFERENCES */

/*     [1] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires ( 4*K - 2 )*K*N + 19/3*K*K*K + 1/2*K*K */
/*     + 43/6*K - 4 floating point operations. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAEST). */

/*     KEYWORDS */

/*     Elementary matrix operations, orthogonal symplectic matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

/*     Quick return if possible. */

#line 183 "MB04QF.f"
    /* Parameter adjustments */
#line 183 "MB04QF.f"
    v_dim1 = *ldv;
#line 183 "MB04QF.f"
    v_offset = 1 + v_dim1;
#line 183 "MB04QF.f"
    v -= v_offset;
#line 183 "MB04QF.f"
    w_dim1 = *ldw;
#line 183 "MB04QF.f"
    w_offset = 1 + w_dim1;
#line 183 "MB04QF.f"
    w -= w_offset;
#line 183 "MB04QF.f"
    --cs;
#line 183 "MB04QF.f"
    --tau;
#line 183 "MB04QF.f"
    rs_dim1 = *ldrs;
#line 183 "MB04QF.f"
    rs_offset = 1 + rs_dim1;
#line 183 "MB04QF.f"
    rs -= rs_offset;
#line 183 "MB04QF.f"
    t_dim1 = *ldt;
#line 183 "MB04QF.f"
    t_offset = 1 + t_dim1;
#line 183 "MB04QF.f"
    t -= t_offset;
#line 183 "MB04QF.f"
    --dwork;
#line 183 "MB04QF.f"

#line 183 "MB04QF.f"
    /* Function Body */
#line 183 "MB04QF.f"
    if (*n == 0) {
#line 183 "MB04QF.f"
	return 0;
#line 183 "MB04QF.f"
    }

#line 186 "MB04QF.f"
    lcolv = lsame_(storev, "C", (ftnlen)1, (ftnlen)1);
#line 187 "MB04QF.f"
    lcolw = lsame_(storew, "C", (ftnlen)1, (ftnlen)1);

#line 189 "MB04QF.f"
    k2 = *k + *k;
#line 190 "MB04QF.f"
    pr1 = 0;
#line 191 "MB04QF.f"
    pr2 = pr1 + *k;
#line 192 "MB04QF.f"
    pr3 = pr2 + *k;
#line 193 "MB04QF.f"
    ps1 = pr3 + *k;
#line 194 "MB04QF.f"
    ps2 = ps1 + *k;
#line 195 "MB04QF.f"
    ps3 = ps2 + *k;

#line 197 "MB04QF.f"
    pt11 = 0;
#line 198 "MB04QF.f"
    pt12 = pt11 + *k;
#line 199 "MB04QF.f"
    pt13 = pt12 + *k;
#line 200 "MB04QF.f"
    pt21 = pt13 + *k;
#line 201 "MB04QF.f"
    pt22 = pt21 + *k;
#line 202 "MB04QF.f"
    pt23 = pt22 + *k;
#line 203 "MB04QF.f"
    pt31 = pt23 + *k;
#line 204 "MB04QF.f"
    pt32 = pt31 + *k;
#line 205 "MB04QF.f"
    pt33 = pt32 + *k;

#line 207 "MB04QF.f"
    i__1 = *k;
#line 207 "MB04QF.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 208 "MB04QF.f"
	taui = tau[i__];
#line 209 "MB04QF.f"
	vii = v[i__ + i__ * v_dim1];
#line 210 "MB04QF.f"
	v[i__ + i__ * v_dim1] = 1.;
#line 211 "MB04QF.f"
	wii = w[i__ + i__ * w_dim1];
#line 212 "MB04QF.f"
	w[i__ + i__ * w_dim1] = 1.;
#line 213 "MB04QF.f"
	if (wii == 0.) {
#line 214 "MB04QF.f"
	    i__2 = i__;
#line 214 "MB04QF.f"
	    for (j = 1; j <= i__2; ++j) {
#line 215 "MB04QF.f"
		t[j + (pt11 + i__) * t_dim1] = 0.;
#line 216 "MB04QF.f"
/* L10: */
#line 216 "MB04QF.f"
	    }
#line 217 "MB04QF.f"
	    i__2 = i__ - 1;
#line 217 "MB04QF.f"
	    for (j = 1; j <= i__2; ++j) {
#line 218 "MB04QF.f"
		t[j + (pt21 + i__) * t_dim1] = 0.;
#line 219 "MB04QF.f"
/* L20: */
#line 219 "MB04QF.f"
	    }
#line 220 "MB04QF.f"
	    i__2 = i__ - 1;
#line 220 "MB04QF.f"
	    for (j = 1; j <= i__2; ++j) {
#line 221 "MB04QF.f"
		t[j + (pt31 + i__) * t_dim1] = 0.;
#line 222 "MB04QF.f"
/* L30: */
#line 222 "MB04QF.f"
	    }
#line 223 "MB04QF.f"
	    i__2 = i__ - 1;
#line 223 "MB04QF.f"
	    for (j = 1; j <= i__2; ++j) {
#line 224 "MB04QF.f"
		rs[j + (ps1 + i__) * rs_dim1] = 0.;
#line 225 "MB04QF.f"
/* L40: */
#line 225 "MB04QF.f"
	    }
#line 226 "MB04QF.f"
	} else {

/*           Treat first Householder reflection. */

#line 230 "MB04QF.f"
	    if (lcolv && lcolw) {

/*              Compute t1 = -wii * W(i:n,1:i-1)' * W(i:n,i). */

#line 234 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 234 "MB04QF.f"
		i__3 = i__ - 1;
#line 234 "MB04QF.f"
		d__1 = -wii;
#line 234 "MB04QF.f"
		dgemv_("Transpose", &i__2, &i__3, &d__1, &w[i__ + w_dim1], 
			ldw, &w[i__ + i__ * w_dim1], &c__1, &c_b11, &dwork[1],
			 &c__1, (ftnlen)9);

/*              Compute t2 = -wii * V(i:n,1:i-1)' * W(i:n,i). */

#line 239 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 239 "MB04QF.f"
		i__3 = i__ - 1;
#line 239 "MB04QF.f"
		d__1 = -wii;
#line 239 "MB04QF.f"
		dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1], 
			ldv, &w[i__ + i__ * w_dim1], &c__1, &c_b11, &dwork[*k 
			+ 1], &c__1, (ftnlen)9);
#line 241 "MB04QF.f"
	    } else if (lcolv) {

/*              Compute t1 = -wii * W(1:i-1,i:n) * W(i,i:n)'. */

#line 245 "MB04QF.f"
		i__2 = i__ - 1;
#line 245 "MB04QF.f"
		i__3 = *n - i__ + 1;
#line 245 "MB04QF.f"
		d__1 = -wii;
#line 245 "MB04QF.f"
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &w[i__ * w_dim1 + 
			1], ldw, &w[i__ + i__ * w_dim1], ldw, &c_b11, &dwork[
			1], &c__1, (ftnlen)12);

/*              Compute t2 = -wii * V(i:n,1:i-1)' * W(i,i:n)'. */

#line 250 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 250 "MB04QF.f"
		i__3 = i__ - 1;
#line 250 "MB04QF.f"
		d__1 = -wii;
#line 250 "MB04QF.f"
		dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1], 
			ldv, &w[i__ + i__ * w_dim1], ldw, &c_b11, &dwork[*k + 
			1], &c__1, (ftnlen)9);
#line 252 "MB04QF.f"
	    } else if (lcolw) {

/*              Compute t1 = -wii * W(i:n,1:i-1)' * W(i:n,i). */

#line 256 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 256 "MB04QF.f"
		i__3 = i__ - 1;
#line 256 "MB04QF.f"
		d__1 = -wii;
#line 256 "MB04QF.f"
		dgemv_("Transpose", &i__2, &i__3, &d__1, &w[i__ + w_dim1], 
			ldw, &w[i__ + i__ * w_dim1], &c__1, &c_b11, &dwork[1],
			 &c__1, (ftnlen)9);

/*              Compute t2 = -wii * V(1:i-1,i:n) * W(i:n,i). */

#line 261 "MB04QF.f"
		i__2 = i__ - 1;
#line 261 "MB04QF.f"
		i__3 = *n - i__ + 1;
#line 261 "MB04QF.f"
		d__1 = -wii;
#line 261 "MB04QF.f"
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &v[i__ * v_dim1 + 
			1], ldv, &w[i__ + i__ * w_dim1], &c__1, &c_b11, &
			dwork[*k + 1], &c__1, (ftnlen)12);
#line 263 "MB04QF.f"
	    } else {

/*              Compute t1 = -wii * W(1:i-1,i:n) * W(i,i:n)'. */

#line 267 "MB04QF.f"
		i__2 = i__ - 1;
#line 267 "MB04QF.f"
		i__3 = *n - i__ + 1;
#line 267 "MB04QF.f"
		d__1 = -wii;
#line 267 "MB04QF.f"
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &w[i__ * w_dim1 + 
			1], ldw, &w[i__ + i__ * w_dim1], ldw, &c_b11, &dwork[
			1], &c__1, (ftnlen)12);

/*              Compute t2 = -wii * V(1:i-1,i:n) * W(i,i:n)'. */

#line 272 "MB04QF.f"
		i__2 = i__ - 1;
#line 272 "MB04QF.f"
		i__3 = *n - i__ + 1;
#line 272 "MB04QF.f"
		d__1 = -wii;
#line 272 "MB04QF.f"
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &v[i__ * v_dim1 + 
			1], ldv, &w[i__ + i__ * w_dim1], ldw, &c_b11, &dwork[*
			k + 1], &c__1, (ftnlen)12);
#line 274 "MB04QF.f"
	    }

/*           T11(1:i-1,i) := T11(1:i-1,1:i-1)*t1 + T13(1:i-1,1:i-1)*t2 */

#line 278 "MB04QF.f"
	    i__2 = i__ - 1;
#line 278 "MB04QF.f"
	    dcopy_(&i__2, &dwork[1], &c__1, &t[(pt11 + i__) * t_dim1 + 1], &
		    c__1);
#line 279 "MB04QF.f"
	    i__2 = i__ - 1;
#line 279 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt11 + 1) *
		     t_dim1 + 1], ldt, &t[(pt11 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 281 "MB04QF.f"
	    i__2 = i__ - 1;
#line 281 "MB04QF.f"
	    dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt13 + i__) * t_dim1 + 1]
		    , &c__1);
#line 282 "MB04QF.f"
	    i__2 = i__ - 1;
#line 282 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt13 + 1) *
		     t_dim1 + 1], ldt, &t[(pt13 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 284 "MB04QF.f"
	    i__2 = i__ - 1;
#line 284 "MB04QF.f"
	    daxpy_(&i__2, &c_b49, &t[(pt13 + i__) * t_dim1 + 1], &c__1, &t[(
		    pt11 + i__) * t_dim1 + 1], &c__1);
#line 285 "MB04QF.f"
	    t[i__ + (pt11 + i__) * t_dim1] = -wii;

#line 287 "MB04QF.f"
	    if (i__ > 1) {

/*              T21(1:i-1,i) := T21(1:i-1,1:i-1)*t1 + T23(1:i-1,1:i-1)*t2 */

#line 291 "MB04QF.f"
		i__2 = i__ - 2;
#line 291 "MB04QF.f"
		dcopy_(&i__2, &dwork[2], &c__1, &t[(pt21 + i__) * t_dim1 + 1],
			 &c__1);
#line 292 "MB04QF.f"
		i__2 = i__ - 2;
#line 292 "MB04QF.f"
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt21 + 
			2) * t_dim1 + 1], ldt, &t[(pt21 + i__) * t_dim1 + 1], 
			&c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 294 "MB04QF.f"
		t[i__ - 1 + (pt21 + i__) * t_dim1] = 0.;
#line 295 "MB04QF.f"
		i__2 = i__ - 1;
#line 295 "MB04QF.f"
		dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt23 + i__) * t_dim1 
			+ 1], &c__1);
#line 296 "MB04QF.f"
		i__2 = i__ - 1;
#line 296 "MB04QF.f"
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt23 + 
			1) * t_dim1 + 1], ldt, &t[(pt23 + i__) * t_dim1 + 1], 
			&c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 298 "MB04QF.f"
		i__2 = i__ - 1;
#line 298 "MB04QF.f"
		daxpy_(&i__2, &c_b49, &t[(pt23 + i__) * t_dim1 + 1], &c__1, &
			t[(pt21 + i__) * t_dim1 + 1], &c__1);

/*              T31(1:i-1,i) := T31(1:i-1,1:i-1)*t1 + T33(1:i-1,1:i-1)*t2 */

#line 302 "MB04QF.f"
		i__2 = i__ - 2;
#line 302 "MB04QF.f"
		dcopy_(&i__2, &dwork[2], &c__1, &t[(pt31 + i__) * t_dim1 + 1],
			 &c__1);
#line 303 "MB04QF.f"
		i__2 = i__ - 2;
#line 303 "MB04QF.f"
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt31 + 
			2) * t_dim1 + 1], ldt, &t[(pt31 + i__) * t_dim1 + 1], 
			&c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 305 "MB04QF.f"
		t[i__ - 1 + (pt31 + i__) * t_dim1] = 0.;
#line 306 "MB04QF.f"
		i__2 = i__ - 1;
#line 306 "MB04QF.f"
		dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt33 + i__) * t_dim1 
			+ 1], &c__1);
#line 307 "MB04QF.f"
		i__2 = i__ - 1;
#line 307 "MB04QF.f"
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt33 + 
			1) * t_dim1 + 1], ldt, &t[(pt33 + i__) * t_dim1 + 1], 
			&c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 309 "MB04QF.f"
		i__2 = i__ - 1;
#line 309 "MB04QF.f"
		daxpy_(&i__2, &c_b49, &t[(pt33 + i__) * t_dim1 + 1], &c__1, &
			t[(pt31 + i__) * t_dim1 + 1], &c__1);

/*              S1(1:i-1,i) := S1(1:i-1,1:i-1)*t1 + S3(1:i-1,1:i-1)*t2 */

#line 313 "MB04QF.f"
		i__2 = i__ - 2;
#line 313 "MB04QF.f"
		dcopy_(&i__2, &dwork[2], &c__1, &rs[(ps1 + i__) * rs_dim1 + 1]
			, &c__1);
#line 314 "MB04QF.f"
		i__2 = i__ - 2;
#line 314 "MB04QF.f"
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps1 + 
			2) * rs_dim1 + 1], ldrs, &rs[(ps1 + i__) * rs_dim1 + 
			1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 316 "MB04QF.f"
		rs[i__ - 1 + (ps1 + i__) * rs_dim1] = 0.;
#line 317 "MB04QF.f"
		i__2 = i__ - 1;
#line 317 "MB04QF.f"
		dcopy_(&i__2, &dwork[*k + 1], &c__1, &rs[(ps3 + i__) * 
			rs_dim1 + 1], &c__1);
#line 318 "MB04QF.f"
		i__2 = i__ - 1;
#line 318 "MB04QF.f"
		dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps3 + 
			1) * rs_dim1 + 1], ldrs, &rs[(ps3 + i__) * rs_dim1 + 
			1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 320 "MB04QF.f"
		i__2 = i__ - 1;
#line 320 "MB04QF.f"
		daxpy_(&i__2, &c_b49, &rs[(ps3 + i__) * rs_dim1 + 1], &c__1, &
			rs[(ps1 + i__) * rs_dim1 + 1], &c__1);
#line 321 "MB04QF.f"
	    }
#line 322 "MB04QF.f"
	}

/*        Treat Givens rotation. */

#line 326 "MB04QF.f"
	cm1 = cs[(i__ << 1) - 1] - 1.;
#line 327 "MB04QF.f"
	if (lcolw) {
#line 328 "MB04QF.f"
	    dcopy_(&i__, &w[i__ + w_dim1], ldw, &dwork[1], &c__1);
#line 329 "MB04QF.f"
	} else {
#line 330 "MB04QF.f"
	    dcopy_(&i__, &w[i__ * w_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 331 "MB04QF.f"
	}
#line 332 "MB04QF.f"
	if (lcolv) {
#line 333 "MB04QF.f"
	    i__2 = i__ - 1;
#line 333 "MB04QF.f"
	    dcopy_(&i__2, &v[i__ + v_dim1], ldv, &dwork[*k + 1], &c__1);
#line 334 "MB04QF.f"
	} else {
#line 335 "MB04QF.f"
	    i__2 = i__ - 1;
#line 335 "MB04QF.f"
	    dcopy_(&i__2, &v[i__ * v_dim1 + 1], &c__1, &dwork[*k + 1], &c__1);
#line 336 "MB04QF.f"
	}

/*        R1(1:i,i) = T11(1:i,1:i) * dwork(1:i) */
/*                    + [ T13(1:i-1,1:i-1) * dwork(k+1:k+i-1); 0 ] */

#line 341 "MB04QF.f"
	dcopy_(&i__, &dwork[1], &c__1, &rs[(pr1 + i__) * rs_dim1 + 1], &c__1);
#line 342 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__, &t[(pt11 + 1) * 
		t_dim1 + 1], ldt, &rs[(pr1 + i__) * rs_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 344 "MB04QF.f"
	i__2 = i__ - 1;
#line 344 "MB04QF.f"
	dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt13 + i__) * t_dim1 + 1], &
		c__1);
#line 345 "MB04QF.f"
	i__2 = i__ - 1;
#line 345 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt13 + 1) * 
		t_dim1 + 1], ldt, &t[(pt13 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 347 "MB04QF.f"
	i__2 = i__ - 1;
#line 347 "MB04QF.f"
	daxpy_(&i__2, &c_b49, &t[(pt13 + i__) * t_dim1 + 1], &c__1, &rs[(pr1 
		+ i__) * rs_dim1 + 1], &c__1);

/*        R2(1:i-1,i) = T21(1:i-1,2:i) * W(i,2:i) */
/*                      + T23(1:i-1,1:i-1) * V(i,1:i-1) */

#line 352 "MB04QF.f"
	i__2 = i__ - 1;
#line 352 "MB04QF.f"
	dcopy_(&i__2, &dwork[2], &c__1, &rs[(pr2 + i__) * rs_dim1 + 1], &c__1)
		;
#line 353 "MB04QF.f"
	i__2 = i__ - 1;
#line 353 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt21 + 2) * 
		t_dim1 + 1], ldt, &rs[(pr2 + i__) * rs_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 355 "MB04QF.f"
	i__2 = i__ - 1;
#line 355 "MB04QF.f"
	dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt23 + i__) * t_dim1 + 1], &
		c__1);
#line 356 "MB04QF.f"
	i__2 = i__ - 1;
#line 356 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt23 + 1) * 
		t_dim1 + 1], ldt, &t[(pt23 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 358 "MB04QF.f"
	i__2 = i__ - 1;
#line 358 "MB04QF.f"
	daxpy_(&i__2, &c_b49, &t[(pt23 + i__) * t_dim1 + 1], &c__1, &rs[(pr2 
		+ i__) * rs_dim1 + 1], &c__1);

/*        R3(1:i-1,i) = T31(1:i-1,2:i) * dwork(2:i) */
/*                      + T33(1:i-1,1:i-1) * dwork(k+1:k+i-1) */

#line 363 "MB04QF.f"
	i__2 = i__ - 1;
#line 363 "MB04QF.f"
	dcopy_(&i__2, &dwork[2], &c__1, &rs[(pr3 + i__) * rs_dim1 + 1], &c__1)
		;
#line 364 "MB04QF.f"
	i__2 = i__ - 1;
#line 364 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt31 + 2) * 
		t_dim1 + 1], ldt, &rs[(pr3 + i__) * rs_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 366 "MB04QF.f"
	i__2 = i__ - 1;
#line 366 "MB04QF.f"
	dcopy_(&i__2, &dwork[*k + 1], &c__1, &t[(pt33 + i__) * t_dim1 + 1], &
		c__1);
#line 367 "MB04QF.f"
	i__2 = i__ - 1;
#line 367 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt33 + 1) * 
		t_dim1 + 1], ldt, &t[(pt33 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 369 "MB04QF.f"
	i__2 = i__ - 1;
#line 369 "MB04QF.f"
	daxpy_(&i__2, &c_b49, &t[(pt33 + i__) * t_dim1 + 1], &c__1, &rs[(pr3 
		+ i__) * rs_dim1 + 1], &c__1);

/*        S2(1:i-1,i) = S1(1:i-1,2:i) * dwork(2:i) */
/*                      + S3(1:i-1,1:i-1) * dwork(k+1:k+i-1) */

#line 374 "MB04QF.f"
	i__2 = i__ - 1;
#line 374 "MB04QF.f"
	dcopy_(&i__2, &dwork[2], &c__1, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1)
		;
#line 375 "MB04QF.f"
	i__2 = i__ - 1;
#line 375 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps1 + 2) * 
		rs_dim1 + 1], ldrs, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 377 "MB04QF.f"
	i__2 = i__ - 1;
#line 377 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps3 + 1) * 
		rs_dim1 + 1], ldrs, &dwork[*k + 1], &c__1, (ftnlen)5, (ftnlen)
		12, (ftnlen)8);
#line 379 "MB04QF.f"
	i__2 = i__ - 1;
#line 379 "MB04QF.f"
	daxpy_(&i__2, &c_b49, &dwork[*k + 1], &c__1, &rs[(ps2 + i__) * 
		rs_dim1 + 1], &c__1);
#line 380 "MB04QF.f"
	rs[i__ + (ps2 + i__) * rs_dim1] = -cs[i__ * 2];

/*        T12(1:i,i) = [ R1(1:i-1,1:i-1)*S2(1:i-1,i); 0 ] */
/*                     + (c-1) * R1(1:i,i) */

#line 385 "MB04QF.f"
	i__2 = i__ - 1;
#line 385 "MB04QF.f"
	dcopy_(&i__2, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1, &t[(pt12 + i__) *
		 t_dim1 + 1], &c__1);
#line 386 "MB04QF.f"
	i__2 = i__ - 1;
#line 386 "MB04QF.f"
	dscal_(&i__2, &cm1, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1);
#line 387 "MB04QF.f"
	i__2 = i__ - 1;
#line 387 "MB04QF.f"
	dscal_(&i__2, &cs[i__ * 2], &t[(pt12 + i__) * t_dim1 + 1], &c__1);
#line 388 "MB04QF.f"
	i__2 = i__ - 1;
#line 388 "MB04QF.f"
	dcopy_(&i__2, &t[(pt12 + i__) * t_dim1 + 1], &c__1, &t[(pt22 + i__) * 
		t_dim1 + 1], &c__1);
#line 389 "MB04QF.f"
	i__2 = i__ - 1;
#line 389 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(pr1 + 1) * 
		rs_dim1 + 1], ldrs, &t[(pt12 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 391 "MB04QF.f"
	t[i__ + (pt12 + i__) * t_dim1] = 0.;
#line 392 "MB04QF.f"
	daxpy_(&i__, &cm1, &rs[(pr1 + i__) * rs_dim1 + 1], &c__1, &t[(pt12 + 
		i__) * t_dim1 + 1], &c__1);

/*        T22(1:i-1,i) = R2(1:i-1,1:i-1)*S2(1:i-1,i) + (c-1)*R2(1:i-1,i) */

#line 396 "MB04QF.f"
	if (i__ > 1) {
#line 396 "MB04QF.f"
	    i__2 = i__ - 2;
#line 396 "MB04QF.f"
	    dcopy_(&i__2, &t[(pt22 + i__) * t_dim1 + 2], &c__1, &t[(pt32 + 
		    i__) * t_dim1 + 1], &c__1);
#line 396 "MB04QF.f"
	}
#line 398 "MB04QF.f"
	i__2 = i__ - 1;
#line 398 "MB04QF.f"
	dtrmv_("Upper", "No transpose", "Unit diagonal", &i__2, &rs[(pr2 + 1) 
		* rs_dim1 + 1], ldrs, &t[(pt22 + i__) * t_dim1 + 1], &c__1, (
		ftnlen)5, (ftnlen)12, (ftnlen)13);
#line 400 "MB04QF.f"
	i__2 = i__ - 1;
#line 400 "MB04QF.f"
	daxpy_(&i__2, &cm1, &rs[(pr2 + i__) * rs_dim1 + 1], &c__1, &t[(pt22 + 
		i__) * t_dim1 + 1], &c__1);
#line 401 "MB04QF.f"
	t[i__ + (pt22 + i__) * t_dim1] = cm1;

/*        T32(1:i-1,i) = R3(1:i-1,1:i-1)*S2(1:i-1,i) + (c-1)*R3(1:i-1,i) */

#line 405 "MB04QF.f"
	if (i__ > 1) {
#line 406 "MB04QF.f"
	    i__2 = i__ - 2;
#line 406 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-Unit", &i__2, &rs[(pr3 + 2) *
		     rs_dim1 + 1], ldrs, &t[(pt32 + i__) * t_dim1 + 1], &c__1,
		     (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 408 "MB04QF.f"
	    t[i__ - 1 + (pt32 + i__) * t_dim1] = 0.;
#line 409 "MB04QF.f"
	    i__2 = i__ - 1;
#line 409 "MB04QF.f"
	    daxpy_(&i__2, &cm1, &rs[(pr3 + i__) * rs_dim1 + 1], &c__1, &t[(
		    pt32 + i__) * t_dim1 + 1], &c__1);
#line 410 "MB04QF.f"
	}

#line 412 "MB04QF.f"
	if (taui == 0.) {
#line 413 "MB04QF.f"
	    i__2 = i__;
#line 413 "MB04QF.f"
	    for (j = 1; j <= i__2; ++j) {
#line 414 "MB04QF.f"
		t[j + (pt13 + i__) * t_dim1] = 0.;
#line 415 "MB04QF.f"
/* L50: */
#line 415 "MB04QF.f"
	    }
#line 416 "MB04QF.f"
	    i__2 = i__;
#line 416 "MB04QF.f"
	    for (j = 1; j <= i__2; ++j) {
#line 417 "MB04QF.f"
		t[j + (pt23 + i__) * t_dim1] = 0.;
#line 418 "MB04QF.f"
/* L60: */
#line 418 "MB04QF.f"
	    }
#line 419 "MB04QF.f"
	    i__2 = i__;
#line 419 "MB04QF.f"
	    for (j = 1; j <= i__2; ++j) {
#line 420 "MB04QF.f"
		t[j + (pt33 + i__) * t_dim1] = 0.;
#line 421 "MB04QF.f"
/* L70: */
#line 421 "MB04QF.f"
	    }
#line 422 "MB04QF.f"
	    i__2 = i__;
#line 422 "MB04QF.f"
	    for (j = 1; j <= i__2; ++j) {
#line 423 "MB04QF.f"
		rs[j + (ps3 + i__) * rs_dim1] = 0.;
#line 424 "MB04QF.f"
/* L80: */
#line 424 "MB04QF.f"
	    }
#line 425 "MB04QF.f"
	} else {

/*           Treat second Householder reflection. */

#line 429 "MB04QF.f"
	    if (lcolv && lcolw) {

/*              Compute t1 = -tau(i) * W(i:n,1:i)' * V(i:n,i). */

#line 433 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 433 "MB04QF.f"
		d__1 = -taui;
#line 433 "MB04QF.f"
		dgemv_("Transpose", &i__2, &i__, &d__1, &w[i__ + w_dim1], ldw,
			 &v[i__ + i__ * v_dim1], &c__1, &c_b11, &dwork[1], &
			c__1, (ftnlen)9);

/*              Compute t2 = -tau(i) * V(i:n,1:i-1)' * V(i:n,i). */

#line 438 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 438 "MB04QF.f"
		i__3 = i__ - 1;
#line 438 "MB04QF.f"
		d__1 = -taui;
#line 438 "MB04QF.f"
		dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1], 
			ldv, &v[i__ + i__ * v_dim1], &c__1, &c_b11, &dwork[k2 
			+ 1], &c__1, (ftnlen)9);
#line 440 "MB04QF.f"
	    } else if (lcolv) {

/*              Compute t1 = -tau(i) * W(1:i,i:n) * V(i:n,i). */

#line 444 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 444 "MB04QF.f"
		d__1 = -taui;
#line 444 "MB04QF.f"
		dgemv_("No Transpose", &i__, &i__2, &d__1, &w[i__ * w_dim1 + 
			1], ldw, &v[i__ + i__ * v_dim1], &c__1, &c_b11, &
			dwork[1], &c__1, (ftnlen)12);

/*              Compute t2 = -tau(i) * V(i:n,1:i-1)' * V(i:n,i). */

#line 449 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 449 "MB04QF.f"
		i__3 = i__ - 1;
#line 449 "MB04QF.f"
		d__1 = -taui;
#line 449 "MB04QF.f"
		dgemv_("Transpose", &i__2, &i__3, &d__1, &v[i__ + v_dim1], 
			ldv, &v[i__ + i__ * v_dim1], &c__1, &c_b11, &dwork[k2 
			+ 1], &c__1, (ftnlen)9);
#line 451 "MB04QF.f"
	    } else if (lcolw) {

/*              Compute t1 = -tau(i) * W(i:n,1:i)' * V(i,i:n)'. */

#line 455 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 455 "MB04QF.f"
		d__1 = -taui;
#line 455 "MB04QF.f"
		dgemv_("Transpose", &i__2, &i__, &d__1, &w[i__ + w_dim1], ldw,
			 &v[i__ + i__ * v_dim1], ldv, &c_b11, &dwork[1], &
			c__1, (ftnlen)9);

/*              Compute t2 = -tau(i) * V(1:i-1,i:n) * V(i,i:n)'. */

#line 460 "MB04QF.f"
		i__2 = i__ - 1;
#line 460 "MB04QF.f"
		i__3 = *n - i__ + 1;
#line 460 "MB04QF.f"
		d__1 = -taui;
#line 460 "MB04QF.f"
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &v[i__ * v_dim1 + 
			1], ldv, &v[i__ + i__ * v_dim1], ldv, &c_b11, &dwork[
			k2 + 1], &c__1, (ftnlen)12);
#line 462 "MB04QF.f"
	    } else {

/*              Compute t1 = -tau(i) * W(1:i,i:n) * V(i,i:n)'. */

#line 466 "MB04QF.f"
		i__2 = *n - i__ + 1;
#line 466 "MB04QF.f"
		d__1 = -taui;
#line 466 "MB04QF.f"
		dgemv_("No Transpose", &i__, &i__2, &d__1, &w[i__ * w_dim1 + 
			1], ldw, &v[i__ + i__ * v_dim1], ldv, &c_b11, &dwork[
			1], &c__1, (ftnlen)12);

/*              Compute t2 = -tau(i) * V(1:i-1,i:n) * V(i,i:n)'. */

#line 471 "MB04QF.f"
		i__2 = i__ - 1;
#line 471 "MB04QF.f"
		i__3 = *n - i__ + 1;
#line 471 "MB04QF.f"
		d__1 = -taui;
#line 471 "MB04QF.f"
		dgemv_("No Transpose", &i__2, &i__3, &d__1, &v[i__ * v_dim1 + 
			1], ldv, &v[i__ + i__ * v_dim1], ldv, &c_b11, &dwork[
			k2 + 1], &c__1, (ftnlen)12);
#line 473 "MB04QF.f"
	    }

/*           T13(1:i,i) := T11(1:i,1:i)*t1 - tau(i)*T12(1:i,i) */
/*                                         + [T13(1:i-1,1:i-1)*t2;0] */

#line 478 "MB04QF.f"
	    i__2 = i__ - 1;
#line 478 "MB04QF.f"
	    dcopy_(&i__2, &dwork[k2 + 1], &c__1, &t[(pt13 + i__) * t_dim1 + 1]
		    , &c__1);
#line 479 "MB04QF.f"
	    i__2 = i__ - 1;
#line 479 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt13 + 1) *
		     t_dim1 + 1], ldt, &t[(pt13 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 481 "MB04QF.f"
	    t[i__ + (pt13 + i__) * t_dim1] = 0.;
#line 482 "MB04QF.f"
	    dcopy_(&i__, &dwork[1], &c__1, &dwork[*k + 1], &c__1);
#line 483 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__, &t[(pt11 + 1) * 
		    t_dim1 + 1], ldt, &dwork[*k + 1], &c__1, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8);
#line 485 "MB04QF.f"
	    daxpy_(&i__, &c_b49, &dwork[*k + 1], &c__1, &t[(pt13 + i__) * 
		    t_dim1 + 1], &c__1);
#line 486 "MB04QF.f"
	    d__1 = -taui;
#line 486 "MB04QF.f"
	    daxpy_(&i__, &d__1, &t[(pt12 + i__) * t_dim1 + 1], &c__1, &t[(
		    pt13 + i__) * t_dim1 + 1], &c__1);

/*           T23(1:i,i) := T21(1:i,1:i)*t1 - tau(i)*T22(1:i,i) */
/*                                         + [T23(1:i-1,1:i-1)*t2;0] */

#line 491 "MB04QF.f"
	    i__2 = i__ - 1;
#line 491 "MB04QF.f"
	    dcopy_(&i__2, &dwork[k2 + 1], &c__1, &t[(pt23 + i__) * t_dim1 + 1]
		    , &c__1);
#line 492 "MB04QF.f"
	    i__2 = i__ - 1;
#line 492 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt23 + 1) *
		     t_dim1 + 1], ldt, &t[(pt23 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 494 "MB04QF.f"
	    t[i__ + (pt23 + i__) * t_dim1] = 0.;
#line 495 "MB04QF.f"
	    i__2 = i__ - 1;
#line 495 "MB04QF.f"
	    dcopy_(&i__2, &dwork[2], &c__1, &dwork[*k + 1], &c__1);
#line 496 "MB04QF.f"
	    i__2 = i__ - 1;
#line 496 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt21 + 2) *
		     t_dim1 + 1], ldt, &dwork[*k + 1], &c__1, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8);
#line 498 "MB04QF.f"
	    i__2 = i__ - 1;
#line 498 "MB04QF.f"
	    daxpy_(&i__2, &c_b49, &dwork[*k + 1], &c__1, &t[(pt23 + i__) * 
		    t_dim1 + 1], &c__1);
#line 499 "MB04QF.f"
	    d__1 = -taui;
#line 499 "MB04QF.f"
	    daxpy_(&i__, &d__1, &t[(pt22 + i__) * t_dim1 + 1], &c__1, &t[(
		    pt23 + i__) * t_dim1 + 1], &c__1);

/*           T33(1:i,i) := T31(1:i,1:i)*t1 - tau(i)*T32(1:i,i) */
/*                                         + [T33(1:i-1,1:i-1)*t2;0] */

#line 504 "MB04QF.f"
	    i__2 = i__ - 1;
#line 504 "MB04QF.f"
	    dcopy_(&i__2, &dwork[k2 + 1], &c__1, &t[(pt33 + i__) * t_dim1 + 1]
		    , &c__1);
#line 505 "MB04QF.f"
	    i__2 = i__ - 1;
#line 505 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt33 + 1) *
		     t_dim1 + 1], ldt, &t[(pt33 + i__) * t_dim1 + 1], &c__1, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 507 "MB04QF.f"
	    i__2 = i__ - 1;
#line 507 "MB04QF.f"
	    dcopy_(&i__2, &dwork[2], &c__1, &dwork[*k + 1], &c__1);
#line 508 "MB04QF.f"
	    i__2 = i__ - 1;
#line 508 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &t[(pt31 + 2) *
		     t_dim1 + 1], ldt, &dwork[*k + 1], &c__1, (ftnlen)5, (
		    ftnlen)12, (ftnlen)8);
#line 510 "MB04QF.f"
	    i__2 = i__ - 1;
#line 510 "MB04QF.f"
	    daxpy_(&i__2, &c_b49, &dwork[*k + 1], &c__1, &t[(pt33 + i__) * 
		    t_dim1 + 1], &c__1);
#line 511 "MB04QF.f"
	    i__2 = i__ - 1;
#line 511 "MB04QF.f"
	    d__1 = -taui;
#line 511 "MB04QF.f"
	    daxpy_(&i__2, &d__1, &t[(pt32 + i__) * t_dim1 + 1], &c__1, &t[(
		    pt33 + i__) * t_dim1 + 1], &c__1);
#line 512 "MB04QF.f"
	    t[i__ + (pt33 + i__) * t_dim1] = -taui;

/*           S3(1:i,i) := S1(1:i,1:i)*t1 - tau(i)*S2(1:i,i) */
/*                                       + [S3(1:i-1,1:i-1)*t2;0] */

#line 517 "MB04QF.f"
	    i__2 = i__ - 1;
#line 517 "MB04QF.f"
	    dcopy_(&i__2, &dwork[k2 + 1], &c__1, &rs[(ps3 + i__) * rs_dim1 + 
		    1], &c__1);
#line 518 "MB04QF.f"
	    i__2 = i__ - 1;
#line 518 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps3 + 1) *
		     rs_dim1 + 1], ldrs, &rs[(ps3 + i__) * rs_dim1 + 1], &
		    c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 520 "MB04QF.f"
	    i__2 = i__ - 1;
#line 520 "MB04QF.f"
	    dtrmv_("Upper", "No transpose", "Non-unit", &i__2, &rs[(ps1 + 2) *
		     rs_dim1 + 1], ldrs, &dwork[2], &c__1, (ftnlen)5, (ftnlen)
		    12, (ftnlen)8);
#line 522 "MB04QF.f"
	    i__2 = i__ - 1;
#line 522 "MB04QF.f"
	    daxpy_(&i__2, &c_b49, &dwork[2], &c__1, &rs[(ps3 + i__) * rs_dim1 
		    + 1], &c__1);
#line 523 "MB04QF.f"
	    rs[i__ + (ps3 + i__) * rs_dim1] = 0.;
#line 524 "MB04QF.f"
	    d__1 = -taui;
#line 524 "MB04QF.f"
	    daxpy_(&i__, &d__1, &rs[(ps2 + i__) * rs_dim1 + 1], &c__1, &rs[(
		    ps3 + i__) * rs_dim1 + 1], &c__1);
#line 525 "MB04QF.f"
	}
#line 526 "MB04QF.f"
	w[i__ + i__ * w_dim1] = wii;
#line 527 "MB04QF.f"
	v[i__ + i__ * v_dim1] = vii;
#line 528 "MB04QF.f"
/* L90: */
#line 528 "MB04QF.f"
    }

#line 530 "MB04QF.f"
    return 0;
/* *** Last line of MB04QF *** */
} /* mb04qf_ */

