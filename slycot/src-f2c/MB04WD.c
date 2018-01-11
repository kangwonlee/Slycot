#line 1 "MB04WD.f"
/* MB04WD.f -- translated by f2c (version 20100827).
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

#line 1 "MB04WD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;

/* Subroutine */ int mb04wd_(char *tranq1, char *tranq2, integer *m, integer *
	n, integer *k, doublereal *q1, integer *ldq1, doublereal *q2, integer 
	*ldq2, doublereal *cs, doublereal *tau, doublereal *dwork, integer *
	ldwork, integer *info, ftnlen tranq1_len, ftnlen tranq2_len)
{
    /* System generated locals */
    address a__1[2];
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1[2], i__2, i__3, i__4;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, ib, nb, ki, kk, nx, pdt, pdw, ierr, pdrs;
    static logical ltrq1, ltrq2;
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
    extern /* Subroutine */ int mb04wu_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
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

/*     Blocked version of the SLICOT Library routine MB04WU. */

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
/*             value of LDWORK, MAX(M+N,8*N*NB + 15*NB*NB), where NB is */
/*             the optimal block size determined by the function UE01MD. */
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

/*     [1] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSGSB). */

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

#line 185 "MB04WD.f"
    /* Parameter adjustments */
#line 185 "MB04WD.f"
    q1_dim1 = *ldq1;
#line 185 "MB04WD.f"
    q1_offset = 1 + q1_dim1;
#line 185 "MB04WD.f"
    q1 -= q1_offset;
#line 185 "MB04WD.f"
    q2_dim1 = *ldq2;
#line 185 "MB04WD.f"
    q2_offset = 1 + q2_dim1;
#line 185 "MB04WD.f"
    q2 -= q2_offset;
#line 185 "MB04WD.f"
    --cs;
#line 185 "MB04WD.f"
    --tau;
#line 185 "MB04WD.f"
    --dwork;
#line 185 "MB04WD.f"

#line 185 "MB04WD.f"
    /* Function Body */
#line 185 "MB04WD.f"
    *info = 0;
#line 186 "MB04WD.f"
    ltrq1 = lsame_(tranq1, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranq1, "C", (
	    ftnlen)1, (ftnlen)1);
#line 187 "MB04WD.f"
    ltrq2 = lsame_(tranq2, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranq2, "C", (
	    ftnlen)1, (ftnlen)1);
/* Writing concatenation */
#line 188 "MB04WD.f"
    i__1[0] = 1, a__1[0] = tranq1;
#line 188 "MB04WD.f"
    i__1[1] = 1, a__1[1] = tranq2;
#line 188 "MB04WD.f"
    s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 188 "MB04WD.f"
    nb = ue01md_(&c__1, "MB04WD", ch__1, m, n, k, (ftnlen)6, (ftnlen)2);

/*     Check the scalar input parameters. */

#line 192 "MB04WD.f"
    if (! (ltrq1 || lsame_(tranq1, "N", (ftnlen)1, (ftnlen)1))) {
#line 193 "MB04WD.f"
	*info = -1;
#line 194 "MB04WD.f"
    } else if (! (ltrq2 || lsame_(tranq2, "N", (ftnlen)1, (ftnlen)1))) {
#line 195 "MB04WD.f"
	*info = -2;
#line 196 "MB04WD.f"
    } else if (*m < 0) {
#line 197 "MB04WD.f"
	*info = -3;
#line 198 "MB04WD.f"
    } else if (*n < 0 || *n > *m) {
#line 199 "MB04WD.f"
	*info = -4;
#line 200 "MB04WD.f"
    } else if (*k < 0 || *k > *n) {
#line 201 "MB04WD.f"
	*info = -5;
#line 202 "MB04WD.f"
    } else if (ltrq1 && *ldq1 < max(1,*n) || ! ltrq1 && *ldq1 < max(1,*m)) {
#line 204 "MB04WD.f"
	*info = -7;
#line 205 "MB04WD.f"
    } else if (ltrq2 && *ldq2 < max(1,*n) || ! ltrq2 && *ldq2 < max(1,*m)) {
#line 207 "MB04WD.f"
	*info = -9;
#line 208 "MB04WD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 208 "MB04WD.f"
	i__2 = 1, i__3 = *m + *n;
#line 208 "MB04WD.f"
	if (*ldwork < max(i__2,i__3)) {
/* Computing MAX */
#line 209 "MB04WD.f"
	    i__2 = 1, i__3 = *m + *n;
#line 209 "MB04WD.f"
	    dwork[1] = (doublereal) max(i__2,i__3);
#line 210 "MB04WD.f"
	    *info = -13;
#line 211 "MB04WD.f"
	}
#line 211 "MB04WD.f"
    }

/*     Return if there were illegal values. */

#line 215 "MB04WD.f"
    if (*info != 0) {
#line 216 "MB04WD.f"
	i__2 = -(*info);
#line 216 "MB04WD.f"
	xerbla_("MB04WD", &i__2, (ftnlen)6);
#line 217 "MB04WD.f"
	return 0;
#line 218 "MB04WD.f"
    }

/*     Quick return if possible. */

#line 222 "MB04WD.f"
    if (*n == 0) {
#line 223 "MB04WD.f"
	dwork[1] = 1.;
#line 224 "MB04WD.f"
	return 0;
#line 225 "MB04WD.f"
    }

#line 227 "MB04WD.f"
    nbmin = 2;
#line 228 "MB04WD.f"
    nx = 0;
#line 229 "MB04WD.f"
    wrkopt = *m + *n;
#line 230 "MB04WD.f"
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code. */

/* Computing MAX */
/* Writing concatenation */
#line 234 "MB04WD.f"
	i__1[0] = 1, a__1[0] = tranq1;
#line 234 "MB04WD.f"
	i__1[1] = 1, a__1[1] = tranq2;
#line 234 "MB04WD.f"
	s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 234 "MB04WD.f"
	i__2 = 0, i__3 = ue01md_(&c__3, "MB04WD", ch__1, m, n, k, (ftnlen)6, (
		ftnlen)2);
#line 234 "MB04WD.f"
	nx = max(i__2,i__3);
#line 235 "MB04WD.f"
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked code. */

/* Computing MAX */
#line 239 "MB04WD.f"
	    i__2 = wrkopt, i__3 = (*n << 3) * nb + nb * 15 * nb;
#line 239 "MB04WD.f"
	    wrkopt = max(i__2,i__3);
#line 240 "MB04WD.f"
	    if (*ldwork < wrkopt) {

/*              Not enough workspace to use optimal NB:  reduce NB and */
/*              determine the minimum value of NB. */

#line 245 "MB04WD.f"
		nb = (integer) ((sqrt((doublereal) ((*n << 4) * *n + *ldwork *
			 15)) - (doublereal) (*n << 2)) / 15.);
/* Computing MAX */
/* Writing concatenation */
#line 247 "MB04WD.f"
		i__1[0] = 1, a__1[0] = tranq1;
#line 247 "MB04WD.f"
		i__1[1] = 1, a__1[1] = tranq2;
#line 247 "MB04WD.f"
		s_cat(ch__1, a__1, i__1, &c__2, (ftnlen)2);
#line 247 "MB04WD.f"
		i__2 = 2, i__3 = ue01md_(&c__2, "MB04WD", ch__1, m, n, k, (
			ftnlen)6, (ftnlen)2);
#line 247 "MB04WD.f"
		nbmin = max(i__2,i__3);
#line 249 "MB04WD.f"
	    }
#line 250 "MB04WD.f"
	}
#line 251 "MB04WD.f"
    }

#line 253 "MB04WD.f"
    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block. */
/*        The first kk columns are handled by the block method. */

#line 258 "MB04WD.f"
	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
#line 259 "MB04WD.f"
	i__2 = *k, i__3 = ki + nb;
#line 259 "MB04WD.f"
	kk = min(i__2,i__3);
#line 260 "MB04WD.f"
    } else {
#line 261 "MB04WD.f"
	kk = 0;
#line 262 "MB04WD.f"
    }

/*     Use unblocked code for the last or only block. */

#line 266 "MB04WD.f"
    if (kk < *n) {
#line 266 "MB04WD.f"
	i__2 = *m - kk;
#line 266 "MB04WD.f"
	i__3 = *n - kk;
#line 266 "MB04WD.f"
	i__4 = *k - kk;
#line 266 "MB04WD.f"
	mb04wu_(tranq1, tranq2, &i__2, &i__3, &i__4, &q1[kk + 1 + (kk + 1) * 
		q1_dim1], ldq1, &q2[kk + 1 + (kk + 1) * q2_dim1], ldq2, &cs[(
		kk << 1) + 1], &tau[kk + 1], &dwork[1], ldwork, &ierr, (
		ftnlen)1, (ftnlen)1);
#line 266 "MB04WD.f"
    }

/*     Blocked code. */

#line 273 "MB04WD.f"
    if (kk > 0) {
#line 274 "MB04WD.f"
	pdrs = 1;
#line 275 "MB04WD.f"
	pdt = pdrs + nb * 6 * nb;
#line 276 "MB04WD.f"
	pdw = pdt + nb * 9 * nb;
#line 277 "MB04WD.f"
	if (ltrq1 && ltrq2) {
#line 278 "MB04WD.f"
	    i__2 = -nb;
#line 278 "MB04WD.f"
	    for (i__ = ki + 1; i__2 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__2) {
/* Computing MIN */
#line 279 "MB04WD.f"
		i__3 = nb, i__4 = *k - i__ + 1;
#line 279 "MB04WD.f"
		ib = min(i__3,i__4);
#line 280 "MB04WD.f"
		if (i__ + ib <= *n) {

/*                 Form the triangular factors of the symplectic block */
/*                 reflector SH. */

#line 285 "MB04WD.f"
		    i__3 = *m - i__ + 1;
#line 285 "MB04WD.f"
		    mb04qf_("Forward", "Rowwise", "Rowwise", &i__3, &ib, &q1[
			    i__ + i__ * q1_dim1], ldq1, &q2[i__ + i__ * 
			    q2_dim1], ldq2, &cs[(i__ << 1) - 1], &tau[i__], &
			    dwork[pdrs], &nb, &dwork[pdt], &nb, &dwork[pdw], (
			    ftnlen)7, (ftnlen)7, (ftnlen)7);

/*                 Apply SH to Q1(i+ib:n,i:m) and Q2(i+ib:n,i:m) from */
/*                 the right. */

#line 293 "MB04WD.f"
		    i__3 = *m - i__ + 1;
#line 293 "MB04WD.f"
		    i__4 = *n - i__ - ib + 1;
#line 293 "MB04WD.f"
		    mb04qc_("Zero Structure", "Transpose", "Transpose", "No "\
			    "Transpose", "Forward", "Rowwise", "Rowwise", &
			    i__3, &i__4, &ib, &q1[i__ + i__ * q1_dim1], ldq1, 
			    &q2[i__ + i__ * q2_dim1], ldq2, &dwork[pdrs], &nb,
			     &dwork[pdt], &nb, &q2[i__ + ib + i__ * q2_dim1], 
			    ldq2, &q1[i__ + ib + i__ * q1_dim1], ldq1, &dwork[
			    pdw], (ftnlen)14, (ftnlen)9, (ftnlen)9, (ftnlen)
			    12, (ftnlen)7, (ftnlen)7, (ftnlen)7);
#line 300 "MB04WD.f"
		}

/*              Apply SH to columns i:m of the current block. */

#line 304 "MB04WD.f"
		i__3 = *m - i__ + 1;
#line 304 "MB04WD.f"
		mb04wu_("Transpose", "Transpose", &i__3, &ib, &ib, &q1[i__ + 
			i__ * q1_dim1], ldq1, &q2[i__ + i__ * q2_dim1], ldq2, 
			&cs[(i__ << 1) - 1], &tau[i__], &dwork[1], ldwork, &
			ierr, (ftnlen)9, (ftnlen)9);
#line 307 "MB04WD.f"
/* L10: */
#line 307 "MB04WD.f"
	    }

#line 309 "MB04WD.f"
	} else if (ltrq1) {
#line 310 "MB04WD.f"
	    i__2 = -nb;
#line 310 "MB04WD.f"
	    for (i__ = ki + 1; i__2 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__2) {
/* Computing MIN */
#line 311 "MB04WD.f"
		i__3 = nb, i__4 = *k - i__ + 1;
#line 311 "MB04WD.f"
		ib = min(i__3,i__4);
#line 312 "MB04WD.f"
		if (i__ + ib <= *n) {

/*                 Form the triangular factors of the symplectic block */
/*                 reflector SH. */

#line 317 "MB04WD.f"
		    i__3 = *m - i__ + 1;
#line 317 "MB04WD.f"
		    mb04qf_("Forward", "Rowwise", "Columnwise", &i__3, &ib, &
			    q1[i__ + i__ * q1_dim1], ldq1, &q2[i__ + i__ * 
			    q2_dim1], ldq2, &cs[(i__ << 1) - 1], &tau[i__], &
			    dwork[pdrs], &nb, &dwork[pdt], &nb, &dwork[pdw], (
			    ftnlen)7, (ftnlen)7, (ftnlen)10);

/*                 Apply SH to Q1(i+ib:n,i:m) from the right and to */
/*                 Q2(i:m,i+ib:n) from the left. */

#line 325 "MB04WD.f"
		    i__3 = *m - i__ + 1;
#line 325 "MB04WD.f"
		    i__4 = *n - i__ - ib + 1;
#line 325 "MB04WD.f"
		    mb04qc_("Zero Structure", "No Transpose", "Transpose", 
			    "No Transpose", "Forward", "Rowwise", "Columnwise"
			    , &i__3, &i__4, &ib, &q1[i__ + i__ * q1_dim1], 
			    ldq1, &q2[i__ + i__ * q2_dim1], ldq2, &dwork[pdrs]
			    , &nb, &dwork[pdt], &nb, &q2[i__ + (i__ + ib) * 
			    q2_dim1], ldq2, &q1[i__ + ib + i__ * q1_dim1], 
			    ldq1, &dwork[pdw], (ftnlen)14, (ftnlen)12, (
			    ftnlen)9, (ftnlen)12, (ftnlen)7, (ftnlen)7, (
			    ftnlen)10);
#line 332 "MB04WD.f"
		}

/*              Apply SH to columns/rows i:m of the current block. */

#line 336 "MB04WD.f"
		i__3 = *m - i__ + 1;
#line 336 "MB04WD.f"
		mb04wu_("Transpose", "No Transpose", &i__3, &ib, &ib, &q1[i__ 
			+ i__ * q1_dim1], ldq1, &q2[i__ + i__ * q2_dim1], 
			ldq2, &cs[(i__ << 1) - 1], &tau[i__], &dwork[1], 
			ldwork, &ierr, (ftnlen)9, (ftnlen)12);
#line 339 "MB04WD.f"
/* L20: */
#line 339 "MB04WD.f"
	    }

#line 341 "MB04WD.f"
	} else if (ltrq2) {
#line 342 "MB04WD.f"
	    i__2 = -nb;
#line 342 "MB04WD.f"
	    for (i__ = ki + 1; i__2 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__2) {
/* Computing MIN */
#line 343 "MB04WD.f"
		i__3 = nb, i__4 = *k - i__ + 1;
#line 343 "MB04WD.f"
		ib = min(i__3,i__4);
#line 344 "MB04WD.f"
		if (i__ + ib <= *n) {

/*                 Form the triangular factors of the symplectic block */
/*                 reflector SH. */

#line 349 "MB04WD.f"
		    i__3 = *m - i__ + 1;
#line 349 "MB04WD.f"
		    mb04qf_("Forward", "Columnwise", "Rowwise", &i__3, &ib, &
			    q1[i__ + i__ * q1_dim1], ldq1, &q2[i__ + i__ * 
			    q2_dim1], ldq2, &cs[(i__ << 1) - 1], &tau[i__], &
			    dwork[pdrs], &nb, &dwork[pdt], &nb, &dwork[pdw], (
			    ftnlen)7, (ftnlen)10, (ftnlen)7);

/*                 Apply SH to Q1(i:m,i+ib:n) from the left and to */
/*                 Q2(i+ib:n,i:m) from the right. */

#line 357 "MB04WD.f"
		    i__3 = *m - i__ + 1;
#line 357 "MB04WD.f"
		    i__4 = *n - i__ - ib + 1;
#line 357 "MB04WD.f"
		    mb04qc_("Zero Structure", "Transpose", "No Transpose", 
			    "No Transpose", "Forward", "Columnwise", "Rowwise"
			    , &i__3, &i__4, &ib, &q1[i__ + i__ * q1_dim1], 
			    ldq1, &q2[i__ + i__ * q2_dim1], ldq2, &dwork[pdrs]
			    , &nb, &dwork[pdt], &nb, &q2[i__ + ib + i__ * 
			    q2_dim1], ldq2, &q1[i__ + (i__ + ib) * q1_dim1], 
			    ldq1, &dwork[pdw], (ftnlen)14, (ftnlen)9, (ftnlen)
			    12, (ftnlen)12, (ftnlen)7, (ftnlen)10, (ftnlen)7);
#line 364 "MB04WD.f"
		}

/*              Apply SH to columns/rows i:m of the current block. */

#line 368 "MB04WD.f"
		i__3 = *m - i__ + 1;
#line 368 "MB04WD.f"
		mb04wu_("No Transpose", "Transpose", &i__3, &ib, &ib, &q1[i__ 
			+ i__ * q1_dim1], ldq1, &q2[i__ + i__ * q2_dim1], 
			ldq2, &cs[(i__ << 1) - 1], &tau[i__], &dwork[1], 
			ldwork, &ierr, (ftnlen)12, (ftnlen)9);
#line 371 "MB04WD.f"
/* L30: */
#line 371 "MB04WD.f"
	    }

#line 373 "MB04WD.f"
	} else {
#line 374 "MB04WD.f"
	    i__2 = -nb;
#line 374 "MB04WD.f"
	    for (i__ = ki + 1; i__2 < 0 ? i__ >= 1 : i__ <= 1; i__ += i__2) {
/* Computing MIN */
#line 375 "MB04WD.f"
		i__3 = nb, i__4 = *k - i__ + 1;
#line 375 "MB04WD.f"
		ib = min(i__3,i__4);
#line 376 "MB04WD.f"
		if (i__ + ib <= *n) {

/*                 Form the triangular factors of the symplectic block */
/*                 reflector SH. */

#line 381 "MB04WD.f"
		    i__3 = *m - i__ + 1;
#line 381 "MB04WD.f"
		    mb04qf_("Forward", "Columnwise", "Columnwise", &i__3, &ib,
			     &q1[i__ + i__ * q1_dim1], ldq1, &q2[i__ + i__ * 
			    q2_dim1], ldq2, &cs[(i__ << 1) - 1], &tau[i__], &
			    dwork[pdrs], &nb, &dwork[pdt], &nb, &dwork[pdw], (
			    ftnlen)7, (ftnlen)10, (ftnlen)10);

/*                 Apply SH to Q1(i:m,i+ib:n) and Q2(i:m,i+ib:n) from */
/*                 the left. */

#line 389 "MB04WD.f"
		    i__3 = *m - i__ + 1;
#line 389 "MB04WD.f"
		    i__4 = *n - i__ - ib + 1;
#line 389 "MB04WD.f"
		    mb04qc_("Zero Structure", "No Transpose", "No Transpose", 
			    "No Transpose", "Forward", "Columnwise", "Column"\
			    "wise", &i__3, &i__4, &ib, &q1[i__ + i__ * q1_dim1]
			    , ldq1, &q2[i__ + i__ * q2_dim1], ldq2, &dwork[
			    pdrs], &nb, &dwork[pdt], &nb, &q2[i__ + (i__ + ib)
			     * q2_dim1], ldq2, &q1[i__ + (i__ + ib) * q1_dim1]
			    , ldq1, &dwork[pdw], (ftnlen)14, (ftnlen)12, (
			    ftnlen)12, (ftnlen)12, (ftnlen)7, (ftnlen)10, (
			    ftnlen)10);
#line 396 "MB04WD.f"
		}

/*              Apply SH to rows i:m of the current block. */

#line 400 "MB04WD.f"
		i__3 = *m - i__ + 1;
#line 400 "MB04WD.f"
		mb04wu_("No Transpose", "No Transpose", &i__3, &ib, &ib, &q1[
			i__ + i__ * q1_dim1], ldq1, &q2[i__ + i__ * q2_dim1], 
			ldq2, &cs[(i__ << 1) - 1], &tau[i__], &dwork[1], 
			ldwork, &ierr, (ftnlen)12, (ftnlen)12);
#line 403 "MB04WD.f"
/* L40: */
#line 403 "MB04WD.f"
	    }
#line 404 "MB04WD.f"
	}
#line 405 "MB04WD.f"
    }

#line 407 "MB04WD.f"
    dwork[1] = (doublereal) wrkopt;

#line 409 "MB04WD.f"
    return 0;
/* *** Last line of MB04WD *** */
} /* mb04wd_ */

