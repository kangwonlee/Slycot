#line 1 "MB04QC.f"
/* MB04QC.f -- translated by f2c (version 20100827).
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

#line 1 "MB04QC.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b22 = 1.;
static doublereal c_b53 = 0.;
static doublereal c_b367 = -1.;

/* Subroutine */ int mb04qc_(char *struct__, char *trana, char *tranb, char *
	tranq, char *direct, char *storev, char *storew, integer *m, integer *
	n, integer *k, doublereal *v, integer *ldv, doublereal *w, integer *
	ldw, doublereal *rs, integer *ldrs, doublereal *t, integer *ldt, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	dwork, ftnlen struct_len, ftnlen trana_len, ftnlen tranb_len, ftnlen 
	tranq_len, ftnlen direct_len, ftnlen storev_len, ftnlen storew_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, rs_dim1, rs_offset, t_dim1, 
	    t_offset, v_dim1, v_offset, w_dim1, w_offset, i__1;

    /* Local variables */
    static integer i__, pr1, pr2, pr3, ps1, ps2, ps3, pt11, pt12, pt13, pt21, 
	    pt22, pt23, pt31, pt32, pt33, pdw1, pdw2, pdw3, pdw4, pdw5, pdw6, 
	    pdw7, pdw8, pdw9;
    static logical la1b1;
    static doublereal fact;
    static logical ltra, ltrb, ltrq;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer itemp;
    static logical lcolv, lcolw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrmm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen);


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

/*     To apply the orthogonal symplectic block reflector */

/*              [  I+V*T*V'  V*R*S*V'  ] */
/*         Q =  [                      ] */
/*              [ -V*R*S*V'  I+V*T*V'  ] */

/*     or its transpose to a real 2m-by-n matrix [ op(A); op(B) ] from */
/*     the left. */
/*     The k-by-k upper triangular blocks of the matrices */

/*                                 [ S1 ]       [ T11 T12 T13 ] */
/*         R  = [ R1 R2 R3 ],  S = [ S2 ],  T = [ T21 T22 T23 ], */
/*                                 [ S3 ]       [ T31 T32 T33 ] */

/*     with R2 unit and S1, R3, T21, T31, T32 strictly upper triangular, */
/*     are stored rowwise in the arrays RS and T, respectively. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     STRUCT  CHARACTER*1 */
/*             Specifies the structure of the first blocks of A and B: */
/*             = 'Z':  the leading K-by-N submatrices of op(A) and op(B) */
/*                     are (implicitly) assumed to be zero; */
/*             = 'N';  no structure to mention. */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op( A ) as follows: */
/*             = 'N':  op( A ) = A; */
/*             = 'T':  op( A ) = A'; */
/*             = 'C':  op( A ) = A'. */

/*     TRANB   CHARACTER*1 */
/*             Specifies the form of op( B ) as follows: */
/*             = 'N':  op( B ) = B; */
/*             = 'T':  op( B ) = B'; */
/*             = 'C':  op( B ) = B'. */

/*     DIRECT  CHARACTER*1 */
/*             This is a dummy argument, which is reserved for future */
/*             extensions of this subroutine. Not referenced. */

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
/*             The number of rows of the matrices op(A) and op(B). */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrices op(A) and op(B). */
/*             N >= 0. */

/*     K       (input) INTEGER */
/*             The order of the triangular matrices defining R, S and T. */
/*             M >= K >= 0. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*                     (LDV,K) if STOREV = 'C', */
/*                     (LDV,M) if STOREV = 'R' */
/*             On entry with STOREV = 'C', the leading M-by-K part of */
/*             this array must contain in its columns the vectors which */
/*             define the elementary reflector used to form parts of Q. */
/*             On entry with STOREV = 'R', the leading K-by-M part of */
/*             this array must contain in its rows the vectors which */
/*             define the elementary reflector used to form parts of Q. */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= MAX(1,M),  if STOREV = 'C'; */
/*             LDV >= MAX(1,K),  if STOREV = 'R'. */

/*     W       (input) DOUBLE PRECISION array, dimension */
/*                     (LDW,K) if STOREW = 'C', */
/*                     (LDW,M) if STOREW = 'R' */
/*             On entry with STOREW = 'C', the leading M-by-K part of */
/*             this array must contain in its columns the vectors which */
/*             define the elementary reflector used to form parts of Q. */
/*             On entry with STOREW = 'R', the leading K-by-M part of */
/*             this array must contain in its rows the vectors which */
/*             define the elementary reflector used to form parts of Q. */

/*     LDW     INTEGER */
/*             The leading dimension of the array W. */
/*             LDW >= MAX(1,M),  if STOREW = 'C'; */
/*             LDW >= MAX(1,K),  if STOREW = 'R'. */

/*     RS      (input) DOUBLE PRECISION array, dimension (K,6*K) */
/*             On entry, the leading K-by-6*K part of this array must */
/*             contain the upper triangular matrices defining the factors */
/*             R and S of the symplectic block reflector Q. The */
/*             (strictly) lower portions of this array are not */
/*             referenced. */

/*     LDRS    INTEGER */
/*             The leading dimension of the array RS.  LDRS >= MAX(1,K). */

/*     T       (input) DOUBLE PRECISION array, dimension (K,9*K) */
/*             On entry, the leading K-by-9*K part of this array must */
/*             contain the upper triangular matrices defining the factor */
/*             T of the symplectic block reflector Q. The (strictly) */
/*             lower portions of this array are not referenced. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= MAX(1,K). */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDA,N) if TRANA = 'N', */
/*                     (LDA,M) if TRANA = 'C' or TRANA = 'T' */
/*             On entry with TRANA = 'N', the leading M-by-N part of this */
/*             array must contain the matrix A. */
/*             On entry with TRANA = 'T' or TRANA = 'C', the leading */
/*             N-by-M part of this array must contain the matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,M),  if TRANA = 'N'; */
/*             LDA >= MAX(1,N),  if TRANA = 'C' or TRANA = 'T'. */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDB,N) if TRANB = 'N', */
/*                     (LDB,M) if TRANB = 'C' or TRANB = 'T' */
/*             On entry with TRANB = 'N', the leading M-by-N part of this */
/*             array must contain the matrix B. */
/*             On entry with TRANB = 'T' or TRANB = 'C', the leading */
/*             N-by-M part of this array must contain the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,M),  if TRANB = 'N'; */
/*             LDB >= MAX(1,N),  if TRANB = 'C' or TRANB = 'T'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK), where */
/*             LDWORK >= 8*N*K,   if STRUCT = 'Z', */
/*             LDWORK >= 9*N*K,   if STRUCT = 'N'. */

/*     REFERENCES */

/*     [1] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 16*( M - K )*N + ( 26*K - 4 )*K*N floating */
/*     point operations if STRUCT = 'Z' and additional ( 12*K + 2 )*K*N */
/*     floating point operations if STRUCT = 'N'. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAESB). */

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

#line 233 "MB04QC.f"
    /* Parameter adjustments */
#line 233 "MB04QC.f"
    v_dim1 = *ldv;
#line 233 "MB04QC.f"
    v_offset = 1 + v_dim1;
#line 233 "MB04QC.f"
    v -= v_offset;
#line 233 "MB04QC.f"
    w_dim1 = *ldw;
#line 233 "MB04QC.f"
    w_offset = 1 + w_dim1;
#line 233 "MB04QC.f"
    w -= w_offset;
#line 233 "MB04QC.f"
    rs_dim1 = *ldrs;
#line 233 "MB04QC.f"
    rs_offset = 1 + rs_dim1;
#line 233 "MB04QC.f"
    rs -= rs_offset;
#line 233 "MB04QC.f"
    t_dim1 = *ldt;
#line 233 "MB04QC.f"
    t_offset = 1 + t_dim1;
#line 233 "MB04QC.f"
    t -= t_offset;
#line 233 "MB04QC.f"
    a_dim1 = *lda;
#line 233 "MB04QC.f"
    a_offset = 1 + a_dim1;
#line 233 "MB04QC.f"
    a -= a_offset;
#line 233 "MB04QC.f"
    b_dim1 = *ldb;
#line 233 "MB04QC.f"
    b_offset = 1 + b_dim1;
#line 233 "MB04QC.f"
    b -= b_offset;
#line 233 "MB04QC.f"
    --dwork;
#line 233 "MB04QC.f"

#line 233 "MB04QC.f"
    /* Function Body */
#line 233 "MB04QC.f"
    if (*m <= 0 || *n <= 0) {
#line 233 "MB04QC.f"
	return 0;
#line 233 "MB04QC.f"
    }
#line 235 "MB04QC.f"
    la1b1 = lsame_(struct__, "N", (ftnlen)1, (ftnlen)1);
#line 236 "MB04QC.f"
    lcolv = lsame_(storev, "C", (ftnlen)1, (ftnlen)1);
#line 237 "MB04QC.f"
    lcolw = lsame_(storew, "C", (ftnlen)1, (ftnlen)1);
#line 238 "MB04QC.f"
    ltra = lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || lsame_(trana, "C", (
	    ftnlen)1, (ftnlen)1);
#line 239 "MB04QC.f"
    ltrb = lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranb, "C", (
	    ftnlen)1, (ftnlen)1);
#line 240 "MB04QC.f"
    ltrq = lsame_(tranq, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranq, "C", (
	    ftnlen)1, (ftnlen)1);

#line 242 "MB04QC.f"
    pr1 = 1;
#line 243 "MB04QC.f"
    pr2 = pr1 + *k;
#line 244 "MB04QC.f"
    pr3 = pr2 + *k;
#line 245 "MB04QC.f"
    ps1 = pr3 + *k;
#line 246 "MB04QC.f"
    ps2 = ps1 + *k;
#line 247 "MB04QC.f"
    ps3 = ps2 + *k;
#line 248 "MB04QC.f"
    pt11 = 1;
#line 249 "MB04QC.f"
    pt12 = pt11 + *k;
#line 250 "MB04QC.f"
    pt13 = pt12 + *k;
#line 251 "MB04QC.f"
    pt21 = pt13 + *k;
#line 252 "MB04QC.f"
    pt22 = pt21 + *k;
#line 253 "MB04QC.f"
    pt23 = pt22 + *k;
#line 254 "MB04QC.f"
    pt31 = pt23 + *k;
#line 255 "MB04QC.f"
    pt32 = pt31 + *k;
#line 256 "MB04QC.f"
    pt33 = pt32 + *k;
#line 257 "MB04QC.f"
    pdw1 = 1;
#line 258 "MB04QC.f"
    pdw2 = pdw1 + *n * *k;
#line 259 "MB04QC.f"
    pdw3 = pdw2 + *n * *k;
#line 260 "MB04QC.f"
    pdw4 = pdw3 + *n * *k;
#line 261 "MB04QC.f"
    pdw5 = pdw4 + *n * *k;
#line 262 "MB04QC.f"
    pdw6 = pdw5 + *n * *k;
#line 263 "MB04QC.f"
    pdw7 = pdw6 + *n * *k;
#line 264 "MB04QC.f"
    pdw8 = pdw7 + *n * *k;
#line 265 "MB04QC.f"
    pdw9 = pdw8 + *n * *k;

/*     Update the matrix A. */

#line 269 "MB04QC.f"
    if (la1b1) {

/*        NZ1) DW7 := A1' */

#line 273 "MB04QC.f"
	if (ltra) {
#line 274 "MB04QC.f"
	    i__1 = *k;
#line 274 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 275 "MB04QC.f"
		dcopy_(n, &a[i__ * a_dim1 + 1], &c__1, &dwork[pdw7 + (i__ - 1)
			 * *n], &c__1);
#line 276 "MB04QC.f"
/* L10: */
#line 276 "MB04QC.f"
	    }
#line 277 "MB04QC.f"
	} else {
#line 278 "MB04QC.f"
	    i__1 = *n;
#line 278 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 279 "MB04QC.f"
		dcopy_(k, &a[i__ * a_dim1 + 1], &c__1, &dwork[pdw7 + i__ - 1],
			 n);
#line 280 "MB04QC.f"
/* L20: */
#line 280 "MB04QC.f"
	    }
#line 281 "MB04QC.f"
	}

/*        NZ2) DW1 := DW7*W1 */

#line 285 "MB04QC.f"
	i__1 = *n * *k;
#line 285 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw1], &c__1);
#line 286 "MB04QC.f"
	if (lcolw) {
#line 287 "MB04QC.f"
	    dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22, &w[
		    w_offset], ldw, &dwork[pdw1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)4);
#line 289 "MB04QC.f"
	} else {
#line 290 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &w[
		    w_offset], ldw, &dwork[pdw1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)4);
#line 292 "MB04QC.f"
	}

/*        NZ3) DW2 := DW7*V1 */

#line 296 "MB04QC.f"
	i__1 = *n * *k;
#line 296 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw2], &c__1);
#line 297 "MB04QC.f"
	if (lcolv) {
#line 298 "MB04QC.f"
	    dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22, &v[
		    v_offset], ldv, &dwork[pdw2], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)4);
#line 300 "MB04QC.f"
	} else {
#line 301 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &v[
		    v_offset], ldv, &dwork[pdw2], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)4);
#line 303 "MB04QC.f"
	}
#line 304 "MB04QC.f"
	fact = 1.;
#line 305 "MB04QC.f"
    } else {
#line 306 "MB04QC.f"
	fact = 0.;
#line 307 "MB04QC.f"
    }

/*     1) DW1 := A2'*W2 */

#line 311 "MB04QC.f"
    if (*m > *k) {
#line 312 "MB04QC.f"
	if (ltra && lcolw) {
#line 313 "MB04QC.f"
	    i__1 = *m - *k;
#line 313 "MB04QC.f"
	    dgemm_("No Transpose", "No Transpose", n, k, &i__1, &c_b22, &a[(*
		    k + 1) * a_dim1 + 1], lda, &w[*k + 1 + w_dim1], ldw, &
		    fact, &dwork[pdw1], n, (ftnlen)12, (ftnlen)12);
#line 316 "MB04QC.f"
	} else if (ltra) {
#line 317 "MB04QC.f"
	    i__1 = *m - *k;
#line 317 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", n, k, &i__1, &c_b22, &a[(*k + 
		    1) * a_dim1 + 1], lda, &w[(*k + 1) * w_dim1 + 1], ldw, &
		    fact, &dwork[pdw1], n, (ftnlen)12, (ftnlen)9);
#line 320 "MB04QC.f"
	} else if (lcolw) {
#line 321 "MB04QC.f"
	    i__1 = *m - *k;
#line 321 "MB04QC.f"
	    dgemm_("Transpose", "No Transpose", n, k, &i__1, &c_b22, &a[*k + 
		    1 + a_dim1], lda, &w[*k + 1 + w_dim1], ldw, &fact, &dwork[
		    pdw1], n, (ftnlen)9, (ftnlen)12);
#line 324 "MB04QC.f"
	} else {
#line 325 "MB04QC.f"
	    i__1 = *m - *k;
#line 325 "MB04QC.f"
	    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, &w[(*k + 1) * w_dim1 + 1], ldw, &fact, &
		    dwork[pdw1], n, (ftnlen)9, (ftnlen)9);
#line 328 "MB04QC.f"
	}
#line 329 "MB04QC.f"
    } else if (! la1b1) {
#line 330 "MB04QC.f"
	dlaset_("All", n, k, &c_b53, &c_b53, &dwork[pdw1], n, (ftnlen)3);
#line 331 "MB04QC.f"
    }

/*     2) DW2 := A2'*V2 */

#line 335 "MB04QC.f"
    if (*m > *k) {
#line 336 "MB04QC.f"
	if (ltra && lcolv) {
#line 337 "MB04QC.f"
	    i__1 = *m - *k;
#line 337 "MB04QC.f"
	    dgemm_("No Transpose", "No Transpose", n, k, &i__1, &c_b22, &a[(*
		    k + 1) * a_dim1 + 1], lda, &v[*k + 1 + v_dim1], ldv, &
		    fact, &dwork[pdw2], n, (ftnlen)12, (ftnlen)12);
#line 340 "MB04QC.f"
	} else if (ltra) {
#line 341 "MB04QC.f"
	    i__1 = *m - *k;
#line 341 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", n, k, &i__1, &c_b22, &a[(*k + 
		    1) * a_dim1 + 1], lda, &v[(*k + 1) * v_dim1 + 1], ldv, &
		    fact, &dwork[pdw2], n, (ftnlen)12, (ftnlen)9);
#line 344 "MB04QC.f"
	} else if (lcolv) {
#line 345 "MB04QC.f"
	    i__1 = *m - *k;
#line 345 "MB04QC.f"
	    dgemm_("Transpose", "No Transpose", n, k, &i__1, &c_b22, &a[*k + 
		    1 + a_dim1], lda, &v[*k + 1 + v_dim1], ldv, &fact, &dwork[
		    pdw2], n, (ftnlen)9, (ftnlen)12);
#line 348 "MB04QC.f"
	} else {
#line 349 "MB04QC.f"
	    i__1 = *m - *k;
#line 349 "MB04QC.f"
	    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, &v[(*k + 1) * v_dim1 + 1], ldv, &fact, &
		    dwork[pdw2], n, (ftnlen)9, (ftnlen)9);
#line 352 "MB04QC.f"
	}
#line 353 "MB04QC.f"
    } else if (! la1b1) {
#line 354 "MB04QC.f"
	dlaset_("All", n, k, &c_b53, &c_b53, &dwork[pdw2], n, (ftnlen)3);
#line 355 "MB04QC.f"
    }

#line 357 "MB04QC.f"
    if (ltrq) {

/*        3) DW3 := DW1*T11 */

#line 361 "MB04QC.f"
	i__1 = *n * *k;
#line 361 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw3], &c__1);
#line 362 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt11 * t_dim1 + 1], ldt, &dwork[pdw3], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        4) DW4 := DW2*T31 */

#line 367 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 367 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
#line 368 "MB04QC.f"
	i__1 = *k - 1;
#line 368 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &t[(pt31 + 1) * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        5) DW3 := DW3 + DW4 */

#line 373 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 373 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3 + *n], &c__1);

#line 375 "MB04QC.f"
	if (la1b1) {

/*           NZ4) DW8 := DW7*T21 */

#line 379 "MB04QC.f"
	    i__1 = *n * (*k - 1);
#line 379 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
#line 380 "MB04QC.f"
	    i__1 = *k - 1;
#line 380 "MB04QC.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &
		    c_b22, &t[(pt21 + 1) * t_dim1 + 1], ldt, &dwork[pdw8], n, 
		    (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ5) DW3 := DW3 + DW8 */

#line 385 "MB04QC.f"
	    i__1 = *n * (*k - 1);
#line 385 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3 + *n], &
		    c__1);
#line 386 "MB04QC.f"
	}

/*        6) DW4 := DW1*T12 */

#line 390 "MB04QC.f"
	i__1 = *n * *k;
#line 390 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw4], &c__1);
#line 391 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt12 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        7) DW5 := DW2*T32 */

#line 396 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 396 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
#line 397 "MB04QC.f"
	i__1 = *k - 1;
#line 397 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &t[(pt32 + 1) * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        8) DW4 := DW4 + DW5 */

#line 402 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 402 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4 + *n], &c__1);

#line 404 "MB04QC.f"
	if (la1b1) {

/*           NZ6) DW8 := DW7*T22 */

#line 408 "MB04QC.f"
	    i__1 = *n * *k;
#line 408 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
#line 409 "MB04QC.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22,
		     &t[pt22 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ7) DW4 := DW4 + DW8 */

#line 414 "MB04QC.f"
	    i__1 = *n * *k;
#line 414 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw4], &c__1);
#line 415 "MB04QC.f"
	}

/*        9) DW5 := DW2*T33 */

#line 419 "MB04QC.f"
	i__1 = *n * *k;
#line 419 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
#line 420 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt33 * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        10) DW6 := DW1*T13 */

#line 425 "MB04QC.f"
	i__1 = *n * *k;
#line 425 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw6], &c__1);
#line 426 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt13 * t_dim1 + 1], ldt, &dwork[pdw6], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        11) DW5 := DW5 + DW6 */

#line 431 "MB04QC.f"
	i__1 = *n * *k;
#line 431 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw6], &c__1, &dwork[pdw5], &c__1);

#line 433 "MB04QC.f"
	if (la1b1) {

/*           NZ8) DW8 := DW7*T23 */

#line 437 "MB04QC.f"
	    i__1 = *n * *k;
#line 437 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
#line 438 "MB04QC.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22,
		     &t[pt23 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ9) DW5 := DW5 + DW8 */

#line 443 "MB04QC.f"
	    i__1 = *n * *k;
#line 443 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw5], &c__1);
#line 444 "MB04QC.f"
	}

/*        12) DW1 := DW1*R1 */

#line 448 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[pr1 * rs_dim1 + 1], ldrs, &dwork[pdw1], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        13) DW2 := DW2*R3 */

#line 453 "MB04QC.f"
	i__1 = *k - 1;
#line 453 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &rs[(pr3 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw2], n, (ftnlen)
		5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        14) DW1 := DW1 + DW2 */

#line 458 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 458 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw2], &c__1, &dwork[pdw1 + *n], &c__1);

#line 460 "MB04QC.f"
	if (la1b1) {

/*           NZ10) DW7 := DW7*R2 */

#line 464 "MB04QC.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &
		    rs[pr2 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)4);

/*           NZ11) DW1 := DW1 + DW7 */

#line 469 "MB04QC.f"
	    i__1 = *n * *k;
#line 469 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw7], &c__1, &dwork[pdw1], &c__1);
#line 470 "MB04QC.f"
	}

/*        Swap Pointers PDW1 <-> PDW2 */

#line 474 "MB04QC.f"
	itemp = pdw2;
#line 475 "MB04QC.f"
	pdw2 = pdw1;
#line 476 "MB04QC.f"
	pdw1 = itemp;
#line 477 "MB04QC.f"
    } else {

/*        3) DW3 := DW1*T11' */

#line 481 "MB04QC.f"
	i__1 = *n * *k;
#line 481 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw3], &c__1);
#line 482 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt11 * t_dim1 + 1], ldt, &dwork[pdw3], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        4) DW4 := DW2*T13' */

#line 487 "MB04QC.f"
	i__1 = *n * *k;
#line 487 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
#line 488 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt13 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        5) DW3 := DW3 + DW4 */

#line 493 "MB04QC.f"
	i__1 = *n * *k;
#line 493 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3], &c__1);

#line 495 "MB04QC.f"
	if (la1b1) {

/*           NZ4) DW8 := DW7*T12' */

#line 499 "MB04QC.f"
	    i__1 = *n * *k;
#line 499 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
#line 500 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    t[pt12 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ5) DW3 := DW3 + DW8 */

#line 505 "MB04QC.f"
	    i__1 = *n * *k;
#line 505 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3], &c__1);
#line 506 "MB04QC.f"
	}

/*        6) DW4 := DW2*T23' */

#line 510 "MB04QC.f"
	i__1 = *n * *k;
#line 510 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
#line 511 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt23 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        7) DW5 := DW1*T21' */

#line 516 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 516 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1 + *n], &c__1, &dwork[pdw5], &c__1);
#line 517 "MB04QC.f"
	i__1 = *k - 1;
#line 517 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		t[(pt21 + 1) * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        8) DW4 := DW4 + DW5 */

#line 522 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 522 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

#line 524 "MB04QC.f"
	if (la1b1) {

/*           NZ6) DW8 := DW7*T22' */

#line 528 "MB04QC.f"
	    i__1 = *n * *k;
#line 528 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
#line 529 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    t[pt22 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ7) DW4 := DW4 + DW8 */

#line 534 "MB04QC.f"
	    i__1 = *n * *k;
#line 534 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw4], &c__1);
#line 535 "MB04QC.f"
	}

/*        9) DW5 := DW2*T33' */

#line 539 "MB04QC.f"
	i__1 = *n * *k;
#line 539 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
#line 540 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt33 * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        10) DW6 := DW1*T31' */

#line 545 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 545 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1 + *n], &c__1, &dwork[pdw6], &c__1);
#line 546 "MB04QC.f"
	i__1 = *k - 1;
#line 546 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		t[(pt31 + 1) * t_dim1 + 1], ldt, &dwork[pdw6], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        11) DW5 := DW5 + DW6 */

#line 551 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 551 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw6], &c__1, &dwork[pdw5], &c__1);

#line 553 "MB04QC.f"
	if (la1b1) {

/*           NZ8) DW8 := DW7*T32' */

#line 557 "MB04QC.f"
	    i__1 = *n * (*k - 1);
#line 557 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw7 + *n], &c__1, &dwork[pdw8], &c__1);
#line 558 "MB04QC.f"
	    i__1 = *k - 1;
#line 558 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &
		    c_b22, &t[(pt32 + 1) * t_dim1 + 1], ldt, &dwork[pdw8], n, 
		    (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ9) DW5 := DW5 + DW8 */

#line 563 "MB04QC.f"
	    i__1 = *n * (*k - 1);
#line 563 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw5], &c__1);
#line 564 "MB04QC.f"
	}

/*        12) DW1 := DW1*S1' */

#line 568 "MB04QC.f"
	i__1 = *k - 1;
#line 568 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		rs[(ps1 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw1 + *n], n, (
		ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        13) DW2 := DW2*S3' */

#line 573 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &rs[
		ps3 * rs_dim1 + 1], ldrs, &dwork[pdw2], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        14) DW2 := DW1 + DW2 */

#line 578 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 578 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw1 + *n], &c__1, &dwork[pdw2], &c__1);

#line 580 "MB04QC.f"
	if (la1b1) {

/*           NZ10) DW7 := DW7*S2' */

#line 584 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    rs[ps2 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ11) DW2 := DW2 + DW7 */

#line 589 "MB04QC.f"
	    i__1 = *n * *k;
#line 589 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw7], &c__1, &dwork[pdw2], &c__1);
#line 590 "MB04QC.f"
	}
#line 591 "MB04QC.f"
    }

#line 593 "MB04QC.f"
    if (la1b1) {

/*        NZ12) DW9 := B1' */

#line 597 "MB04QC.f"
	if (ltrb) {
#line 598 "MB04QC.f"
	    i__1 = *k;
#line 598 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 599 "MB04QC.f"
		dcopy_(n, &b[i__ * b_dim1 + 1], &c__1, &dwork[pdw9 + (i__ - 1)
			 * *n], &c__1);
#line 600 "MB04QC.f"
/* L30: */
#line 600 "MB04QC.f"
	    }
#line 601 "MB04QC.f"
	} else {
#line 602 "MB04QC.f"
	    i__1 = *n;
#line 602 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 603 "MB04QC.f"
		dcopy_(k, &b[i__ * b_dim1 + 1], &c__1, &dwork[pdw9 + i__ - 1],
			 n);
#line 604 "MB04QC.f"
/* L40: */
#line 604 "MB04QC.f"
	    }
#line 605 "MB04QC.f"
	}

/*        NZ13) DW1 := DW9*W1 */

#line 609 "MB04QC.f"
	i__1 = *n * *k;
#line 609 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw1], &c__1);
#line 610 "MB04QC.f"
	if (lcolw) {
#line 611 "MB04QC.f"
	    dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22, &w[
		    w_offset], ldw, &dwork[pdw1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)4);
#line 613 "MB04QC.f"
	} else {
#line 614 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &w[
		    w_offset], ldw, &dwork[pdw1], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)4);
#line 616 "MB04QC.f"
	}

/*        NZ14) DW6 := DW9*V1 */

#line 620 "MB04QC.f"
	i__1 = *n * *k;
#line 620 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw6], &c__1);
#line 621 "MB04QC.f"
	if (lcolv) {
#line 622 "MB04QC.f"
	    dtrmm_("Right", "Lower", "No transpose", "Unit", n, k, &c_b22, &v[
		    v_offset], ldv, &dwork[pdw6], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)12, (ftnlen)4);
#line 624 "MB04QC.f"
	} else {
#line 625 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &v[
		    v_offset], ldv, &dwork[pdw6], n, (ftnlen)5, (ftnlen)5, (
		    ftnlen)9, (ftnlen)4);
#line 627 "MB04QC.f"
	}
#line 628 "MB04QC.f"
    }

/*     15) DW1 := B2'*W2 */

#line 632 "MB04QC.f"
    if (*m > *k) {
#line 633 "MB04QC.f"
	if (ltrb && lcolw) {
#line 634 "MB04QC.f"
	    i__1 = *m - *k;
#line 634 "MB04QC.f"
	    dgemm_("No Transpose", "No Transpose", n, k, &i__1, &c_b22, &b[(*
		    k + 1) * b_dim1 + 1], ldb, &w[*k + 1 + w_dim1], ldw, &
		    fact, &dwork[pdw1], n, (ftnlen)12, (ftnlen)12);
#line 637 "MB04QC.f"
	} else if (ltrb) {

/*           Critical Position */

#line 641 "MB04QC.f"
	    i__1 = *m - *k;
#line 641 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", n, k, &i__1, &c_b22, &b[(*k + 
		    1) * b_dim1 + 1], ldb, &w[(*k + 1) * w_dim1 + 1], ldw, &
		    fact, &dwork[pdw1], n, (ftnlen)12, (ftnlen)9);
#line 644 "MB04QC.f"
	} else if (lcolw) {
#line 645 "MB04QC.f"
	    i__1 = *m - *k;
#line 645 "MB04QC.f"
	    dgemm_("Transpose", "No Transpose", n, k, &i__1, &c_b22, &b[*k + 
		    1 + b_dim1], ldb, &w[*k + 1 + w_dim1], ldw, &fact, &dwork[
		    pdw1], n, (ftnlen)9, (ftnlen)12);
#line 648 "MB04QC.f"
	} else {
#line 649 "MB04QC.f"
	    i__1 = *m - *k;
#line 649 "MB04QC.f"
	    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, &w[(*k + 1) * w_dim1 + 1], ldw, &fact, &
		    dwork[pdw1], n, (ftnlen)9, (ftnlen)9);
#line 652 "MB04QC.f"
	}
#line 653 "MB04QC.f"
    } else if (! la1b1) {
#line 654 "MB04QC.f"
	dlaset_("All", n, k, &c_b53, &c_b53, &dwork[pdw1], n, (ftnlen)3);
#line 655 "MB04QC.f"
    }

/*     16) DW6 := B2'*V2 */

#line 659 "MB04QC.f"
    if (*m > *k) {
#line 660 "MB04QC.f"
	if (ltrb && lcolv) {
#line 661 "MB04QC.f"
	    i__1 = *m - *k;
#line 661 "MB04QC.f"
	    dgemm_("No Transpose", "No Transpose", n, k, &i__1, &c_b22, &b[(*
		    k + 1) * b_dim1 + 1], ldb, &v[*k + 1 + v_dim1], ldv, &
		    fact, &dwork[pdw6], n, (ftnlen)12, (ftnlen)12);
#line 664 "MB04QC.f"
	} else if (ltrb) {
#line 665 "MB04QC.f"
	    i__1 = *m - *k;
#line 665 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", n, k, &i__1, &c_b22, &b[(*k + 
		    1) * b_dim1 + 1], ldb, &v[(*k + 1) * v_dim1 + 1], ldv, &
		    fact, &dwork[pdw6], n, (ftnlen)12, (ftnlen)9);
#line 668 "MB04QC.f"
	} else if (lcolv) {
#line 669 "MB04QC.f"
	    i__1 = *m - *k;
#line 669 "MB04QC.f"
	    dgemm_("Transpose", "No Transpose", n, k, &i__1, &c_b22, &b[*k + 
		    1 + b_dim1], ldb, &v[*k + 1 + v_dim1], ldv, &fact, &dwork[
		    pdw6], n, (ftnlen)9, (ftnlen)12);
#line 672 "MB04QC.f"
	} else {
#line 673 "MB04QC.f"
	    i__1 = *m - *k;
#line 673 "MB04QC.f"
	    dgemm_("Transpose", "Transpose", n, k, &i__1, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, &v[(*k + 1) * v_dim1 + 1], ldv, &fact, &
		    dwork[pdw6], n, (ftnlen)9, (ftnlen)9);
#line 676 "MB04QC.f"
	}
#line 677 "MB04QC.f"
    } else if (! la1b1) {
#line 678 "MB04QC.f"
	dlaset_("All", n, k, &c_b53, &c_b53, &dwork[pdw6], n, (ftnlen)3);
#line 679 "MB04QC.f"
    }

#line 681 "MB04QC.f"
    if (ltrq) {

/*        17) DW7 := DW1*R1 */

#line 685 "MB04QC.f"
	i__1 = *n * *k;
#line 685 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw7], &c__1);
#line 686 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[pr1 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        18) DW8 := DW6*R3 */

#line 691 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 691 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw8], &c__1);
#line 692 "MB04QC.f"
	i__1 = *k - 1;
#line 692 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &rs[(pr3 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)
		5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        19) DW7 := DW7 + DW8 */

#line 697 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 697 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw7 + *n], &c__1);

#line 699 "MB04QC.f"
	if (la1b1) {

/*           NZ15) DW8 := DW9*R2 */

#line 703 "MB04QC.f"
	    i__1 = *n * *k;
#line 703 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
#line 704 "MB04QC.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &
		    rs[pr2 * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)4);

/*           NZ16) DW7 := DW7 + DW8 */

#line 709 "MB04QC.f"
	    i__1 = *n * *k;
#line 709 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw7], &c__1);
#line 710 "MB04QC.f"
	}

/*        20) DW8 := DW7*S1 */

#line 714 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 714 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
#line 715 "MB04QC.f"
	i__1 = *k - 1;
#line 715 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &rs[(ps1 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)
		5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        21) DW3 := DW3 - DW8 */

#line 720 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 720 "MB04QC.f"
	daxpy_(&i__1, &c_b367, &dwork[pdw8], &c__1, &dwork[pdw3 + *n], &c__1);

/*        22) DW8 := DW7*S3 */

#line 724 "MB04QC.f"
	i__1 = *n * *k;
#line 724 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
#line 725 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[ps3 * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        23) DW5 := DW5 - DW8 */

#line 730 "MB04QC.f"
	i__1 = *n * *k;
#line 730 "MB04QC.f"
	daxpy_(&i__1, &c_b367, &dwork[pdw8], &c__1, &dwork[pdw5], &c__1);

/*        24) DW7 := DW7*S2 */

#line 734 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b367, &
		rs[ps2 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 736 "MB04QC.f"
    } else {

/*        17) DW7 := DW6*S3' */

#line 740 "MB04QC.f"
	i__1 = *n * *k;
#line 740 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw7], &c__1);
#line 741 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &rs[
		ps3 * rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        18) DW8 := DW1*S1' */

#line 746 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 746 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1 + *n], &c__1, &dwork[pdw8], &c__1);
#line 747 "MB04QC.f"
	i__1 = *k - 1;
#line 747 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		rs[(ps1 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        19) DW7 := DW7 + DW8 */

#line 752 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 752 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw7], &c__1);

#line 754 "MB04QC.f"
	if (la1b1) {

/*           NZ15) DW8 := DW9*S2' */

#line 758 "MB04QC.f"
	    i__1 = *n * *k;
#line 758 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
#line 759 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    rs[ps2 * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ16) DW7 := DW7 + DW8 */

#line 764 "MB04QC.f"
	    i__1 = *n * *k;
#line 764 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw7], &c__1);
#line 765 "MB04QC.f"
	}

/*        20) DW8 := DW7*R1' */

#line 769 "MB04QC.f"
	i__1 = *n * *k;
#line 769 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw7], &c__1, &dwork[pdw8], &c__1);
#line 770 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &rs[
		pr1 * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        21) DW3 := DW3 + DW8 */

#line 775 "MB04QC.f"
	i__1 = *n * *k;
#line 775 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3], &c__1);

/*        22) DW8 := DW7*R3' */

#line 779 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 779 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw7 + *n], &c__1, &dwork[pdw8], &c__1);
#line 780 "MB04QC.f"
	i__1 = *k - 1;
#line 780 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		rs[(pr3 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw8], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        23) DW5 := DW5 + DW8 */

#line 785 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 785 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw5], &c__1);

/*        24) DW7 := DW7*R2' */

#line 789 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &rs[pr2 * 
		rs_dim1 + 1], ldrs, &dwork[pdw7], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
#line 791 "MB04QC.f"
    }

/*     25) A2 := A2 + W2*DW3' */

#line 795 "MB04QC.f"
    if (*m > *k) {
#line 796 "MB04QC.f"
	if (ltra && lcolw) {
#line 797 "MB04QC.f"
	    i__1 = *m - *k;
#line 797 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", n, &i__1, k, &c_b22, &dwork[
		    pdw3], n, &w[*k + 1 + w_dim1], ldw, &c_b22, &a[(*k + 1) * 
		    a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
#line 800 "MB04QC.f"
	} else if (ltra) {
#line 801 "MB04QC.f"
	    i__1 = *m - *k;
#line 801 "MB04QC.f"
	    dgemm_("No Transpose", "No Transpose", n, &i__1, k, &c_b22, &
		    dwork[pdw3], n, &w[(*k + 1) * w_dim1 + 1], ldw, &c_b22, &
		    a[(*k + 1) * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)12);
#line 804 "MB04QC.f"
	} else if (lcolw) {
#line 805 "MB04QC.f"
	    i__1 = *m - *k;
#line 805 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", &i__1, n, k, &c_b22, &w[*k + 
		    1 + w_dim1], ldw, &dwork[pdw3], n, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 808 "MB04QC.f"
	} else {
#line 809 "MB04QC.f"
	    i__1 = *m - *k;
#line 809 "MB04QC.f"
	    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b22, &w[(*k + 1) 
		    * w_dim1 + 1], ldw, &dwork[pdw3], n, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, (ftnlen)9, (ftnlen)9);
#line 812 "MB04QC.f"
	}
#line 813 "MB04QC.f"
    }

/*     26) A2 := A2 + V2*DW5' */

#line 817 "MB04QC.f"
    if (*m > *k) {
#line 818 "MB04QC.f"
	if (ltra && lcolv) {
#line 819 "MB04QC.f"
	    i__1 = *m - *k;
#line 819 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", n, &i__1, k, &c_b22, &dwork[
		    pdw5], n, &v[*k + 1 + v_dim1], ldv, &c_b22, &a[(*k + 1) * 
		    a_dim1 + 1], lda, (ftnlen)12, (ftnlen)9);
#line 822 "MB04QC.f"
	} else if (ltra) {
#line 823 "MB04QC.f"
	    i__1 = *m - *k;
#line 823 "MB04QC.f"
	    dgemm_("No Transpose", "No Transpose", n, &i__1, k, &c_b22, &
		    dwork[pdw5], n, &v[(*k + 1) * v_dim1 + 1], ldv, &c_b22, &
		    a[(*k + 1) * a_dim1 + 1], lda, (ftnlen)12, (ftnlen)12);
#line 826 "MB04QC.f"
	} else if (lcolv) {
#line 827 "MB04QC.f"
	    i__1 = *m - *k;
#line 827 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", &i__1, n, k, &c_b22, &v[*k + 
		    1 + v_dim1], ldv, &dwork[pdw5], n, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, (ftnlen)12, (ftnlen)9);
#line 830 "MB04QC.f"
	} else {
#line 831 "MB04QC.f"
	    i__1 = *m - *k;
#line 831 "MB04QC.f"
	    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b22, &v[(*k + 1) 
		    * v_dim1 + 1], ldv, &dwork[pdw5], n, &c_b22, &a[*k + 1 + 
		    a_dim1], lda, (ftnlen)9, (ftnlen)9);
#line 834 "MB04QC.f"
	}
#line 835 "MB04QC.f"
    }

/*     27) DW4 := DW4 + DW7 */

#line 839 "MB04QC.f"
    i__1 = *n * *k;
#line 839 "MB04QC.f"
    daxpy_(&i__1, &c_b22, &dwork[pdw7], &c__1, &dwork[pdw4], &c__1);

/*     28) DW3 := DW3*W1' */

#line 843 "MB04QC.f"
    if (lcolw) {
#line 844 "MB04QC.f"
	dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &w[
		w_offset], ldw, &dwork[pdw3], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
#line 846 "MB04QC.f"
    } else {
#line 847 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &w[
		w_offset], ldw, &dwork[pdw3], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)12, (ftnlen)4);
#line 849 "MB04QC.f"
    }

/*     29) DW4 := DW4 + DW3 */

#line 853 "MB04QC.f"
    i__1 = *n * *k;
#line 853 "MB04QC.f"
    daxpy_(&i__1, &c_b22, &dwork[pdw3], &c__1, &dwork[pdw4], &c__1);

/*     30) DW5 := DW5*V1' */

#line 857 "MB04QC.f"
    if (lcolv) {
#line 858 "MB04QC.f"
	dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &v[
		v_offset], ldv, &dwork[pdw5], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
#line 860 "MB04QC.f"
    } else {
#line 861 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &v[
		v_offset], ldv, &dwork[pdw5], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)12, (ftnlen)4);
#line 863 "MB04QC.f"
    }

/*     31) DW4 := DW4 + DW5 */

#line 867 "MB04QC.f"
    i__1 = *n * *k;
#line 867 "MB04QC.f"
    daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

/*     32) A1 := A1 + DW4' */

#line 871 "MB04QC.f"
    if (la1b1) {
#line 872 "MB04QC.f"
	if (ltra) {
#line 873 "MB04QC.f"
	    i__1 = *k;
#line 873 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 874 "MB04QC.f"
		daxpy_(n, &c_b22, &dwork[pdw4 + (i__ - 1) * *n], &c__1, &a[
			i__ * a_dim1 + 1], &c__1);
#line 875 "MB04QC.f"
/* L50: */
#line 875 "MB04QC.f"
	    }
#line 876 "MB04QC.f"
	} else {
#line 877 "MB04QC.f"
	    i__1 = *n;
#line 877 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 878 "MB04QC.f"
		daxpy_(k, &c_b22, &dwork[pdw4 + i__ - 1], n, &a[i__ * a_dim1 
			+ 1], &c__1);
#line 879 "MB04QC.f"
/* L60: */
#line 879 "MB04QC.f"
	    }
#line 880 "MB04QC.f"
	}
#line 881 "MB04QC.f"
    } else {
#line 882 "MB04QC.f"
	if (ltra) {
#line 883 "MB04QC.f"
	    i__1 = *k;
#line 883 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 884 "MB04QC.f"
		dcopy_(n, &dwork[pdw4 + (i__ - 1) * *n], &c__1, &a[i__ * 
			a_dim1 + 1], &c__1);
#line 885 "MB04QC.f"
/* L70: */
#line 885 "MB04QC.f"
	    }
#line 886 "MB04QC.f"
	} else {
#line 887 "MB04QC.f"
	    i__1 = *n;
#line 887 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 888 "MB04QC.f"
		dcopy_(k, &dwork[pdw4 + i__ - 1], n, &a[i__ * a_dim1 + 1], &
			c__1);
#line 889 "MB04QC.f"
/* L80: */
#line 889 "MB04QC.f"
	    }
#line 890 "MB04QC.f"
	}
#line 891 "MB04QC.f"
    }

/*     Update the matrix B. */

#line 895 "MB04QC.f"
    if (ltrq) {

/*        33) DW3 := DW1*T11 */

#line 899 "MB04QC.f"
	i__1 = *n * *k;
#line 899 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw3], &c__1);
#line 900 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt11 * t_dim1 + 1], ldt, &dwork[pdw3], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        34) DW4 := DW6*T31 */

#line 905 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 905 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw4], &c__1);
#line 906 "MB04QC.f"
	i__1 = *k - 1;
#line 906 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &t[(pt31 + 1) * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        35) DW3 := DW3 + DW4 */

#line 911 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 911 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3 + *n], &c__1);

#line 913 "MB04QC.f"
	if (la1b1) {

/*           NZ17) DW8 := DW9*T21 */

#line 917 "MB04QC.f"
	    i__1 = *n * (*k - 1);
#line 917 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
#line 918 "MB04QC.f"
	    i__1 = *k - 1;
#line 918 "MB04QC.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &
		    c_b22, &t[(pt21 + 1) * t_dim1 + 1], ldt, &dwork[pdw8], n, 
		    (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ18) DW3 := DW3 + DW8 */

#line 923 "MB04QC.f"
	    i__1 = *n * (*k - 1);
#line 923 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3 + *n], &
		    c__1);
#line 924 "MB04QC.f"
	}

/*        36) DW4 := DW2*S1 */

#line 928 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 928 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
#line 929 "MB04QC.f"
	i__1 = *k - 1;
#line 929 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &rs[(ps1 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw4], n, (ftnlen)
		5, (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        37) DW3 := DW3 + DW4 */

#line 934 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 934 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3 + *n], &c__1);

/*        38) DW4 := DW1*T12 */

#line 938 "MB04QC.f"
	i__1 = *n * *k;
#line 938 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw4], &c__1);
#line 939 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt12 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        38) DW5 := DW6*T32 */

#line 944 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 944 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw5], &c__1);
#line 945 "MB04QC.f"
	i__1 = *k - 1;
#line 945 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, &i__1, &c_b22,
		 &t[(pt32 + 1) * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5,
		 (ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        40) DW4 := DW4 + DW5 */

#line 950 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 950 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4 + *n], &c__1);

#line 952 "MB04QC.f"
	if (la1b1) {

/*           NZ19) DW8 := DW9*T22 */

#line 956 "MB04QC.f"
	    i__1 = *n * *k;
#line 956 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
#line 957 "MB04QC.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22,
		     &t[pt22 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ20) DW4 := DW4 + DW8 */

#line 962 "MB04QC.f"
	    i__1 = *n * *k;
#line 962 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw4], &c__1);
#line 963 "MB04QC.f"
	}

/*        41) DW5 := DW2*S2 */

#line 967 "MB04QC.f"
	i__1 = *n * *k;
#line 967 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
#line 968 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[ps2 * rs_dim1 + 1], ldrs, &dwork[pdw5], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        42) DW4 := DW4 + DW5 */

#line 973 "MB04QC.f"
	i__1 = *n * *k;
#line 973 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

/*        43) DW6 := DW6*T33 */

#line 977 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt33 * t_dim1 + 1], ldt, &dwork[pdw6], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        44) DW1 := DW1*T13 */

#line 982 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt13 * t_dim1 + 1], ldt, &dwork[pdw1], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)12, (ftnlen)8);

/*        45) DW6 := DW6 + DW1 */

#line 987 "MB04QC.f"
	i__1 = *n * *k;
#line 987 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw1], &c__1, &dwork[pdw6], &c__1);

#line 989 "MB04QC.f"
	if (la1b1) {

/*           NZ19) DW9 := DW9*T23 */

#line 993 "MB04QC.f"
	    dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22,
		     &t[pt23 * t_dim1 + 1], ldt, &dwork[pdw9], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)12, (ftnlen)8);

/*           NZ20) DW6 := DW6 + DW9 */

#line 998 "MB04QC.f"
	    i__1 = *n * *k;
#line 998 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw9], &c__1, &dwork[pdw6], &c__1);
#line 999 "MB04QC.f"
	}

/*        46) DW2 := DW2*S3 */

#line 1003 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Non-Unit", n, k, &c_b22, &
		rs[ps3 * rs_dim1 + 1], ldrs, &dwork[pdw2], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)12, (ftnlen)8);

/*        45) DW6 := DW6 + DW2 */

#line 1008 "MB04QC.f"
	i__1 = *n * *k;
#line 1008 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw2], &c__1, &dwork[pdw6], &c__1);
#line 1009 "MB04QC.f"
    } else {

/*        33) DW3 := DW1*T11' */

#line 1013 "MB04QC.f"
	i__1 = *n * *k;
#line 1013 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1], &c__1, &dwork[pdw3], &c__1);
#line 1014 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt11 * t_dim1 + 1], ldt, &dwork[pdw3], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        34) DW4 := DW6*T13' */

#line 1019 "MB04QC.f"
	i__1 = *n * *k;
#line 1019 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw4], &c__1);
#line 1020 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt13 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        35) DW3 := DW3 + DW4 */

#line 1025 "MB04QC.f"
	i__1 = *n * *k;
#line 1025 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw4], &c__1, &dwork[pdw3], &c__1);

#line 1027 "MB04QC.f"
	if (la1b1) {

/*           NZ17) DW8 := DW9*T12' */

#line 1031 "MB04QC.f"
	    i__1 = *n * *k;
#line 1031 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
#line 1032 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    t[pt12 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ18) DW3 := DW3 + DW8 */

#line 1037 "MB04QC.f"
	    i__1 = *n * *k;
#line 1037 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw3], &c__1);
#line 1038 "MB04QC.f"
	}

/*        36) DW4 := DW2*R1' */

#line 1042 "MB04QC.f"
	i__1 = *n * *k;
#line 1042 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw4], &c__1);
#line 1043 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &rs[
		pr1 * rs_dim1 + 1], ldrs, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        37) DW3 := DW3 - DW4 */

#line 1048 "MB04QC.f"
	i__1 = *n * *k;
#line 1048 "MB04QC.f"
	daxpy_(&i__1, &c_b367, &dwork[pdw4], &c__1, &dwork[pdw3], &c__1);

/*        38) DW4 := DW6*T23' */

#line 1052 "MB04QC.f"
	i__1 = *n * *k;
#line 1052 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw6], &c__1, &dwork[pdw4], &c__1);
#line 1053 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt23 * t_dim1 + 1], ldt, &dwork[pdw4], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        39) DW5 := DW1*T21' */

#line 1058 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 1058 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw1 + *n], &c__1, &dwork[pdw5], &c__1);
#line 1059 "MB04QC.f"
	i__1 = *k - 1;
#line 1059 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		t[(pt21 + 1) * t_dim1 + 1], ldt, &dwork[pdw5], n, (ftnlen)5, (
		ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        40) DW4 := DW4 + DW5 */

#line 1064 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 1064 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

#line 1066 "MB04QC.f"
	if (la1b1) {

/*           NZ19) DW8 := DW9*T22' */

#line 1070 "MB04QC.f"
	    i__1 = *n * *k;
#line 1070 "MB04QC.f"
	    dcopy_(&i__1, &dwork[pdw9], &c__1, &dwork[pdw8], &c__1);
#line 1071 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &
		    t[pt22 * t_dim1 + 1], ldt, &dwork[pdw8], n, (ftnlen)5, (
		    ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ20) DW4 := DW4 + DW8 */

#line 1076 "MB04QC.f"
	    i__1 = *n * *k;
#line 1076 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw8], &c__1, &dwork[pdw4], &c__1);
#line 1077 "MB04QC.f"
	}

/*        41) DW5 := DW2*R2' */

#line 1081 "MB04QC.f"
	i__1 = *n * *k;
#line 1081 "MB04QC.f"
	dcopy_(&i__1, &dwork[pdw2], &c__1, &dwork[pdw5], &c__1);
#line 1082 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Unit", n, k, &c_b22, &rs[pr2 * 
		rs_dim1 + 1], ldrs, &dwork[pdw5], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);

/*        42) DW4 := DW4 - DW5 */

#line 1087 "MB04QC.f"
	i__1 = *n * *k;
#line 1087 "MB04QC.f"
	daxpy_(&i__1, &c_b367, &dwork[pdw5], &c__1, &dwork[pdw4], &c__1);

/*        43) DW6 := DW6*T33' */

#line 1091 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, k, &c_b22, &t[
		pt33 * t_dim1 + 1], ldt, &dwork[pdw6], n, (ftnlen)5, (ftnlen)
		5, (ftnlen)9, (ftnlen)8);

/*        44) DW1 := DW1*T31' */

#line 1096 "MB04QC.f"
	i__1 = *k - 1;
#line 1096 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		t[(pt31 + 1) * t_dim1 + 1], ldt, &dwork[pdw1 + *n], n, (
		ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        45) DW6 := DW6 + DW1 */

#line 1101 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 1101 "MB04QC.f"
	daxpy_(&i__1, &c_b22, &dwork[pdw1 + *n], &c__1, &dwork[pdw6], &c__1);

#line 1103 "MB04QC.f"
	if (la1b1) {

/*           NZ19) DW9 := DW9*T32' */

#line 1107 "MB04QC.f"
	    i__1 = *k - 1;
#line 1107 "MB04QC.f"
	    dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &
		    c_b22, &t[(pt32 + 1) * t_dim1 + 1], ldt, &dwork[pdw9 + *n]
		    , n, (ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*           NZ20) DW6 := DW6 + DW9 */

#line 1112 "MB04QC.f"
	    i__1 = *n * (*k - 1);
#line 1112 "MB04QC.f"
	    daxpy_(&i__1, &c_b22, &dwork[pdw9 + *n], &c__1, &dwork[pdw6], &
		    c__1);
#line 1113 "MB04QC.f"
	}

/*        46) DW2 := DW2*R3' */

#line 1117 "MB04QC.f"
	i__1 = *k - 1;
#line 1117 "MB04QC.f"
	dtrmm_("Right", "Upper", "Transpose", "Non-Unit", n, &i__1, &c_b22, &
		rs[(pr3 + 1) * rs_dim1 + 1], ldrs, &dwork[pdw2 + *n], n, (
		ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);

/*        45) DW6 := DW6 - DW2 */

#line 1122 "MB04QC.f"
	i__1 = *n * (*k - 1);
#line 1122 "MB04QC.f"
	daxpy_(&i__1, &c_b367, &dwork[pdw2 + *n], &c__1, &dwork[pdw6], &c__1);
#line 1123 "MB04QC.f"
    }

/*     46) B2 := B2 + W2*DW3' */

#line 1127 "MB04QC.f"
    if (*m > *k) {
#line 1128 "MB04QC.f"
	if (ltrb && lcolw) {
#line 1129 "MB04QC.f"
	    i__1 = *m - *k;
#line 1129 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", n, &i__1, k, &c_b22, &dwork[
		    pdw3], n, &w[*k + 1 + w_dim1], ldw, &c_b22, &b[(*k + 1) * 
		    b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)9);
#line 1132 "MB04QC.f"
	} else if (ltrb) {
#line 1133 "MB04QC.f"
	    i__1 = *m - *k;
#line 1133 "MB04QC.f"
	    dgemm_("No Transpose", "No Transpose", n, &i__1, k, &c_b22, &
		    dwork[pdw3], n, &w[(*k + 1) * w_dim1 + 1], ldw, &c_b22, &
		    b[(*k + 1) * b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)12);
#line 1136 "MB04QC.f"
	} else if (lcolw) {
#line 1137 "MB04QC.f"
	    i__1 = *m - *k;
#line 1137 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", &i__1, n, k, &c_b22, &w[*k + 
		    1 + w_dim1], ldw, &dwork[pdw3], n, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, (ftnlen)12, (ftnlen)9);
#line 1140 "MB04QC.f"
	} else {
#line 1141 "MB04QC.f"
	    i__1 = *m - *k;
#line 1141 "MB04QC.f"
	    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b22, &w[(*k + 1) 
		    * w_dim1 + 1], ldw, &dwork[pdw3], n, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, (ftnlen)9, (ftnlen)9);
#line 1144 "MB04QC.f"
	}
#line 1145 "MB04QC.f"
    }

/*     47) B2 := B2 + V2*DW6' */

#line 1149 "MB04QC.f"
    if (*m > *k) {
#line 1150 "MB04QC.f"
	if (ltrb && lcolv) {
#line 1151 "MB04QC.f"
	    i__1 = *m - *k;
#line 1151 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", n, &i__1, k, &c_b22, &dwork[
		    pdw6], n, &v[*k + 1 + v_dim1], ldv, &c_b22, &b[(*k + 1) * 
		    b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)9);
#line 1154 "MB04QC.f"
	} else if (ltrb) {
#line 1155 "MB04QC.f"
	    i__1 = *m - *k;
#line 1155 "MB04QC.f"
	    dgemm_("No Transpose", "No Transpose", n, &i__1, k, &c_b22, &
		    dwork[pdw6], n, &v[(*k + 1) * v_dim1 + 1], ldv, &c_b22, &
		    b[(*k + 1) * b_dim1 + 1], ldb, (ftnlen)12, (ftnlen)12);
#line 1158 "MB04QC.f"
	} else if (lcolv) {
#line 1159 "MB04QC.f"
	    i__1 = *m - *k;
#line 1159 "MB04QC.f"
	    dgemm_("No Transpose", "Transpose", &i__1, n, k, &c_b22, &v[*k + 
		    1 + v_dim1], ldv, &dwork[pdw6], n, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, (ftnlen)12, (ftnlen)9);
#line 1162 "MB04QC.f"
	} else {
#line 1163 "MB04QC.f"
	    i__1 = *m - *k;
#line 1163 "MB04QC.f"
	    dgemm_("Transpose", "Transpose", &i__1, n, k, &c_b22, &v[(*k + 1) 
		    * v_dim1 + 1], ldv, &dwork[pdw6], n, &c_b22, &b[*k + 1 + 
		    b_dim1], ldb, (ftnlen)9, (ftnlen)9);
#line 1166 "MB04QC.f"
	}
#line 1167 "MB04QC.f"
    }

/*     48) DW3 := DW3*W1' */

#line 1171 "MB04QC.f"
    if (lcolw) {
#line 1172 "MB04QC.f"
	dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &w[
		w_offset], ldw, &dwork[pdw3], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
#line 1174 "MB04QC.f"
    } else {
#line 1175 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &w[
		w_offset], ldw, &dwork[pdw3], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)12, (ftnlen)4);
#line 1177 "MB04QC.f"
    }

/*     49) DW4 := DW4 + DW3 */

#line 1181 "MB04QC.f"
    i__1 = *n * *k;
#line 1181 "MB04QC.f"
    daxpy_(&i__1, &c_b22, &dwork[pdw3], &c__1, &dwork[pdw4], &c__1);

/*     50) DW6 := DW6*V1' */

#line 1185 "MB04QC.f"
    if (lcolv) {
#line 1186 "MB04QC.f"
	dtrmm_("Right", "Lower", "Transpose", "Unit", n, k, &c_b22, &v[
		v_offset], ldv, &dwork[pdw6], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)9, (ftnlen)4);
#line 1188 "MB04QC.f"
    } else {
#line 1189 "MB04QC.f"
	dtrmm_("Right", "Upper", "No Transpose", "Unit", n, k, &c_b22, &v[
		v_offset], ldv, &dwork[pdw6], n, (ftnlen)5, (ftnlen)5, (
		ftnlen)12, (ftnlen)4);
#line 1191 "MB04QC.f"
    }

/*     51) DW4 := DW4 + DW6 */

#line 1195 "MB04QC.f"
    i__1 = *n * *k;
#line 1195 "MB04QC.f"
    daxpy_(&i__1, &c_b22, &dwork[pdw6], &c__1, &dwork[pdw4], &c__1);

/*     52) B1 := B1 + DW4' */

#line 1199 "MB04QC.f"
    if (la1b1) {
#line 1200 "MB04QC.f"
	if (ltrb) {
#line 1201 "MB04QC.f"
	    i__1 = *k;
#line 1201 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1202 "MB04QC.f"
		daxpy_(n, &c_b22, &dwork[pdw4 + (i__ - 1) * *n], &c__1, &b[
			i__ * b_dim1 + 1], &c__1);
#line 1203 "MB04QC.f"
/* L90: */
#line 1203 "MB04QC.f"
	    }
#line 1204 "MB04QC.f"
	} else {
#line 1205 "MB04QC.f"
	    i__1 = *n;
#line 1205 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1206 "MB04QC.f"
		daxpy_(k, &c_b22, &dwork[pdw4 + i__ - 1], n, &b[i__ * b_dim1 
			+ 1], &c__1);
#line 1207 "MB04QC.f"
/* L100: */
#line 1207 "MB04QC.f"
	    }
#line 1208 "MB04QC.f"
	}
#line 1209 "MB04QC.f"
    } else {
#line 1210 "MB04QC.f"
	if (ltrb) {
#line 1211 "MB04QC.f"
	    i__1 = *k;
#line 1211 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1212 "MB04QC.f"
		dcopy_(n, &dwork[pdw4 + (i__ - 1) * *n], &c__1, &b[i__ * 
			b_dim1 + 1], &c__1);
#line 1213 "MB04QC.f"
/* L110: */
#line 1213 "MB04QC.f"
	    }
#line 1214 "MB04QC.f"
	} else {
#line 1215 "MB04QC.f"
	    i__1 = *n;
#line 1215 "MB04QC.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1216 "MB04QC.f"
		dcopy_(k, &dwork[pdw4 + i__ - 1], n, &b[i__ * b_dim1 + 1], &
			c__1);
#line 1217 "MB04QC.f"
/* L120: */
#line 1217 "MB04QC.f"
	    }
#line 1218 "MB04QC.f"
	}
#line 1219 "MB04QC.f"
    }

#line 1221 "MB04QC.f"
    return 0;
/* *** Last line of MB04QC *** */
} /* mb04qc_ */

