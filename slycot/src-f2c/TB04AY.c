#line 1 "TB04AY.f"
/* TB04AY.f -- translated by f2c (version 20100827).
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

#line 1 "TB04AY.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b19 = 1.;

/* Subroutine */ int tb04ay_(integer *n, integer *mwork, integer *pwork, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, integer *ncont, 
	integer *indexd, doublereal *dcoeff, integer *lddcoe, doublereal *
	ucoeff, integer *lduco1, integer *lduco2, doublereal *at, integer *n1,
	 doublereal *tau, doublereal *tol1, doublereal *tol2, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, at_dim1, at_offset, b_dim1, b_offset, c_dim1, 
	    c_offset, d_dim1, d_offset, dcoeff_dim1, dcoeff_offset, 
	    ucoeff_dim1, ucoeff_dim2, ucoeff_offset, i__1, i__2, i__3, i__4, 
	    i__5;

    /* Local variables */
    static integer i__, j, k, l, ib, ic, is, iv, iz, ibi;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer maxm;
    static doublereal temp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), tb01ud_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, integer *, doublereal *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     integer *, ftnlen), tb01zd_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer nminl;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer jwork, lwork, nplus, ivmin1, indcon, iwplus, wrkopt;


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

/*     Calculates the (PWORK x MWORK) transfer matrix T(s), in the form */
/*     of polynomial row vectors over monic least common denominator */
/*     polynomials, of a given state-space representation (ssr).  Each */
/*     such row of T(s) is simply a single-output relatively left prime */
/*     polynomial matrix representation (pmr), so can be calculated by */
/*     applying a simplified version of the Orthogonal Structure */
/*     Theorem to a minimal ssr for the corresponding row of the given */
/*     system: such an ssr is obtained by using the Orthogonal Canon- */
/*     ical Form to first separate out a completely controllable one */
/*     for the overall system and then, for each row in turn, applying */
/*     it again to the resulting dual SIMO system.  The Orthogonal */
/*     Structure Theorem produces non-monic denominator and V:I(s) */
/*     polynomials: this is avoided here by first scaling AT (the */
/*     transpose of the controllable part of A, found in this routine) */
/*     by suitable products of its sub-diagonal elements (these are then */
/*     no longer needed, so freeing the entire lower triangle for */
/*     storing the coefficients of V(s) apart from the leading 1's, */
/*     which are treated implicitly).  These polynomials are calculated */
/*     in reverse order (IW = NMINL - 1,...,1), the monic denominator */
/*     D:I(s) found exactly as if it were V:0(s), and finally the */
/*     numerator vector U:I(s) obtained from the Orthogonal Structure */
/*     Theorem relation. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Separate out controllable subsystem (of order NCONT). */

/*     Workspace: MAX(N, 3*MWORK, PWORK). */

#line 78 "TB04AY.f"
    /* Parameter adjustments */
#line 78 "TB04AY.f"
    a_dim1 = *lda;
#line 78 "TB04AY.f"
    a_offset = 1 + a_dim1;
#line 78 "TB04AY.f"
    a -= a_offset;
#line 78 "TB04AY.f"
    b_dim1 = *ldb;
#line 78 "TB04AY.f"
    b_offset = 1 + b_dim1;
#line 78 "TB04AY.f"
    b -= b_offset;
#line 78 "TB04AY.f"
    c_dim1 = *ldc;
#line 78 "TB04AY.f"
    c_offset = 1 + c_dim1;
#line 78 "TB04AY.f"
    c__ -= c_offset;
#line 78 "TB04AY.f"
    d_dim1 = *ldd;
#line 78 "TB04AY.f"
    d_offset = 1 + d_dim1;
#line 78 "TB04AY.f"
    d__ -= d_offset;
#line 78 "TB04AY.f"
    --indexd;
#line 78 "TB04AY.f"
    dcoeff_dim1 = *lddcoe;
#line 78 "TB04AY.f"
    dcoeff_offset = 1 + dcoeff_dim1;
#line 78 "TB04AY.f"
    dcoeff -= dcoeff_offset;
#line 78 "TB04AY.f"
    ucoeff_dim1 = *lduco1;
#line 78 "TB04AY.f"
    ucoeff_dim2 = *lduco2;
#line 78 "TB04AY.f"
    ucoeff_offset = 1 + ucoeff_dim1 * (1 + ucoeff_dim2);
#line 78 "TB04AY.f"
    ucoeff -= ucoeff_offset;
#line 78 "TB04AY.f"
    at_dim1 = *n1;
#line 78 "TB04AY.f"
    at_offset = 1 + at_dim1;
#line 78 "TB04AY.f"
    at -= at_offset;
#line 78 "TB04AY.f"
    --tau;
#line 78 "TB04AY.f"
    --iwork;
#line 78 "TB04AY.f"
    --dwork;
#line 78 "TB04AY.f"

#line 78 "TB04AY.f"
    /* Function Body */
#line 78 "TB04AY.f"
    tb01ud_("No Z", n, mwork, pwork, &a[a_offset], lda, &b[b_offset], ldb, &
	    c__[c_offset], ldc, ncont, &indcon, &iwork[1], &at[at_offset], &
	    c__1, &tau[1], tol2, &iwork[*n + 1], &dwork[1], ldwork, info, (
	    ftnlen)4);
#line 81 "TB04AY.f"
    wrkopt = (integer) dwork[1];

#line 83 "TB04AY.f"
    is = 1;
#line 84 "TB04AY.f"
    ic = is + *ncont;
#line 85 "TB04AY.f"
    iz = ic;
#line 86 "TB04AY.f"
    ib = ic + *ncont;
#line 87 "TB04AY.f"
    lwork = ib + *mwork * *ncont;
#line 88 "TB04AY.f"
    maxm = max(1,*mwork);

/*     Calculate each row of T(s) in turn. */

#line 92 "TB04AY.f"
    i__1 = *pwork;
#line 92 "TB04AY.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Form the dual of I-th NCONT-order MISO subsystem ... */

#line 96 "TB04AY.f"
	dcopy_(ncont, &c__[i__ + c_dim1], ldc, &dwork[ic], &c__1);

#line 98 "TB04AY.f"
	i__2 = *ncont;
#line 98 "TB04AY.f"
	for (j = 1; j <= i__2; ++j) {
#line 99 "TB04AY.f"
	    dcopy_(ncont, &a[j + a_dim1], lda, &at[j * at_dim1 + 1], &c__1);
#line 100 "TB04AY.f"
	    dcopy_(mwork, &b[j + b_dim1], ldb, &dwork[(j - 1) * maxm + ib], &
		    c__1);
#line 101 "TB04AY.f"
/* L10: */
#line 101 "TB04AY.f"
	}

/*        and separate out its controllable part, giving minimal */
/*        state-space realization for row I. */

/*        Workspace: MWORK*NCONT + 2*NCONT + MAX(NCONT,MWORK). */

#line 108 "TB04AY.f"
	i__2 = *ldwork - lwork + 1;
#line 108 "TB04AY.f"
	tb01zd_("No Z", ncont, mwork, &at[at_offset], n1, &dwork[ic], &dwork[
		ib], &maxm, &nminl, &dwork[iz], &c__1, &tau[1], tol1, &dwork[
		lwork], &i__2, info, (ftnlen)4);
/* Computing MAX */
#line 111 "TB04AY.f"
	i__2 = wrkopt, i__3 = (integer) dwork[lwork] + lwork - 1;
#line 111 "TB04AY.f"
	wrkopt = max(i__2,i__3);

/*        Store degree of (monic) denominator, and leading coefficient */
/*        vector of numerator. */

#line 116 "TB04AY.f"
	indexd[i__] = nminl;
#line 117 "TB04AY.f"
	dcoeff[i__ + dcoeff_dim1] = 1.;
#line 118 "TB04AY.f"
	dcopy_(mwork, &d__[i__ + d_dim1], ldd, &ucoeff[i__ + (ucoeff_dim2 + 1)
		 * ucoeff_dim1], lduco1);

#line 120 "TB04AY.f"
	if (nminl == 1) {

/*           Finish off numerator, denominator for simple case NMINL=1. */

#line 124 "TB04AY.f"
	    temp = -at[at_dim1 + 1];
#line 125 "TB04AY.f"
	    dcoeff[i__ + (dcoeff_dim1 << 1)] = temp;
#line 126 "TB04AY.f"
	    dcopy_(mwork, &d__[i__ + d_dim1], ldd, &ucoeff[i__ + ((
		    ucoeff_dim2 << 1) + 1) * ucoeff_dim1], lduco1);
#line 127 "TB04AY.f"
	    dscal_(mwork, &temp, &ucoeff[i__ + ((ucoeff_dim2 << 1) + 1) * 
		    ucoeff_dim1], lduco1);
#line 128 "TB04AY.f"
	    daxpy_(mwork, &dwork[ic], &dwork[ib], &c__1, &ucoeff[i__ + ((
		    ucoeff_dim2 << 1) + 1) * ucoeff_dim1], lduco1);
#line 130 "TB04AY.f"
	} else if (nminl > 1) {

/*           Set up factors for scaling upper triangle of AT ... */

#line 134 "TB04AY.f"
	    i__2 = nminl - 1;
#line 134 "TB04AY.f"
	    i__3 = *n1 + 1;
#line 134 "TB04AY.f"
	    dcopy_(&i__2, &at[at_dim1 + 2], &i__3, &dwork[ic + 1], &c__1);
#line 135 "TB04AY.f"
	    nplus = nminl + 1;

#line 137 "TB04AY.f"
	    i__2 = is + nminl - 1;
#line 137 "TB04AY.f"
	    for (l = is; l <= i__2; ++l) {
#line 138 "TB04AY.f"
		dwork[l] = 1.;
#line 139 "TB04AY.f"
/* L20: */
#line 139 "TB04AY.f"
	    }

/*           and scale it, row by row, starting with row NMINL. */

#line 143 "TB04AY.f"
	    for (jwork = nminl; jwork >= 1; --jwork) {

#line 145 "TB04AY.f"
		i__2 = nminl;
#line 145 "TB04AY.f"
		for (j = jwork; j <= i__2; ++j) {
#line 146 "TB04AY.f"
		    at[jwork + j * at_dim1] = dwork[is + j - 1] * at[jwork + 
			    j * at_dim1];
#line 147 "TB04AY.f"
/* L30: */
#line 147 "TB04AY.f"
		}

/*              Update scale factors for next row. */

#line 151 "TB04AY.f"
		i__2 = nminl - jwork + 1;
#line 151 "TB04AY.f"
		dscal_(&i__2, &dwork[ic + jwork - 1], &dwork[is + jwork - 1], 
			&c__1);
#line 153 "TB04AY.f"
/* L40: */
#line 153 "TB04AY.f"
	    }

/*           Calculate each monic polynomial V:JWORK(s) in turn: */
/*           K-th coefficient stored as AT(IV,K-1). */

#line 158 "TB04AY.f"
	    i__2 = nminl;
#line 158 "TB04AY.f"
	    for (iv = 2; iv <= i__2; ++iv) {
#line 159 "TB04AY.f"
		jwork = nplus - iv;
#line 160 "TB04AY.f"
		iwplus = jwork + 1;
#line 161 "TB04AY.f"
		ivmin1 = iv - 1;

/*              Set up coefficients due to leading 1's of existing */
/*              V:I(s)'s. */

#line 166 "TB04AY.f"
		i__3 = ivmin1;
#line 166 "TB04AY.f"
		for (k = 1; k <= i__3; ++k) {
#line 167 "TB04AY.f"
		    at[iv + k * at_dim1] = -at[iwplus + (jwork + k) * at_dim1]
			    ;
#line 168 "TB04AY.f"
/* L50: */
#line 168 "TB04AY.f"
		}

#line 170 "TB04AY.f"
		if (iv != 2) {

/*                 Then add contribution from s * V:JWORK+1(s) term. */

#line 174 "TB04AY.f"
		    i__3 = iv - 2;
#line 174 "TB04AY.f"
		    daxpy_(&i__3, &c_b19, &at[ivmin1 + at_dim1], n1, &at[iv + 
			    at_dim1], n1);

/*                 Finally, add effect of lower coefficients of existing */
/*                 V:I(s)'s. */

#line 180 "TB04AY.f"
		    i__3 = ivmin1;
#line 180 "TB04AY.f"
		    for (k = 2; k <= i__3; ++k) {
#line 181 "TB04AY.f"
			i__4 = k - 1;
#line 181 "TB04AY.f"
			i__5 = -(*n1 + 1);
#line 181 "TB04AY.f"
			at[iv + k * at_dim1] -= ddot_(&i__4, &at[iwplus + (
				jwork + 1) * at_dim1], n1, &at[iv - k + 1 + 
				at_dim1], &i__5);
#line 184 "TB04AY.f"
/* L60: */
#line 184 "TB04AY.f"
		    }

#line 186 "TB04AY.f"
		}
#line 187 "TB04AY.f"
/* L70: */
#line 187 "TB04AY.f"
	    }

/*           Determine denominator polynomial D(s) as if it were V:0(s). */

#line 191 "TB04AY.f"
	    i__2 = nplus;
#line 191 "TB04AY.f"
	    for (k = 2; k <= i__2; ++k) {
#line 192 "TB04AY.f"
		dcoeff[i__ + k * dcoeff_dim1] = -at[(k - 1) * at_dim1 + 1];
#line 193 "TB04AY.f"
/* L80: */
#line 193 "TB04AY.f"
	    }

#line 195 "TB04AY.f"
	    i__2 = nminl - 1;
#line 195 "TB04AY.f"
	    daxpy_(&i__2, &c_b19, &at[nminl + at_dim1], n1, &dcoeff[i__ + (
		    dcoeff_dim1 << 1)], lddcoe);

#line 198 "TB04AY.f"
	    i__2 = nplus;
#line 198 "TB04AY.f"
	    for (k = 3; k <= i__2; ++k) {
#line 199 "TB04AY.f"
		i__3 = k - 2;
#line 199 "TB04AY.f"
		i__4 = -(*n1 + 1);
#line 199 "TB04AY.f"
		dcoeff[i__ + k * dcoeff_dim1] -= ddot_(&i__3, &at[at_offset], 
			n1, &at[nminl - k + 3 + at_dim1], &i__4);
#line 201 "TB04AY.f"
/* L90: */
#line 201 "TB04AY.f"
	    }

/*           Scale (B' * Z), stored in DWORK(IB). */

#line 205 "TB04AY.f"
	    ibi = ib;

#line 207 "TB04AY.f"
	    i__2 = nminl;
#line 207 "TB04AY.f"
	    for (l = 1; l <= i__2; ++l) {
#line 208 "TB04AY.f"
		dscal_(mwork, &dwork[is + l - 1], &dwork[ibi], &c__1);
#line 209 "TB04AY.f"
		ibi += maxm;
#line 210 "TB04AY.f"
/* L100: */
#line 210 "TB04AY.f"
	    }

/*           Evaluate numerator polynomial vector (V(s) * B) + (D(s) */
/*           * D:I): first set up coefficients due to D:I and leading */
/*           1's of V(s). */

#line 216 "TB04AY.f"
	    ibi = ib;

#line 218 "TB04AY.f"
	    i__2 = nplus;
#line 218 "TB04AY.f"
	    for (k = 2; k <= i__2; ++k) {
#line 219 "TB04AY.f"
		dcopy_(mwork, &dwork[ibi], &c__1, &ucoeff[i__ + (k * 
			ucoeff_dim2 + 1) * ucoeff_dim1], lduco1);
#line 220 "TB04AY.f"
		daxpy_(mwork, &dcoeff[i__ + k * dcoeff_dim1], &d__[i__ + 
			d_dim1], ldd, &ucoeff[i__ + (k * ucoeff_dim2 + 1) * 
			ucoeff_dim1], lduco1);
#line 222 "TB04AY.f"
		ibi += maxm;
#line 223 "TB04AY.f"
/* L110: */
#line 223 "TB04AY.f"
	    }

/*           Add contribution from lower coefficients of V(s). */

#line 227 "TB04AY.f"
	    i__2 = nplus;
#line 227 "TB04AY.f"
	    for (k = 3; k <= i__2; ++k) {

#line 229 "TB04AY.f"
		i__3 = *mwork;
#line 229 "TB04AY.f"
		for (j = 1; j <= i__3; ++j) {
#line 230 "TB04AY.f"
		    i__4 = k - 2;
#line 230 "TB04AY.f"
		    i__5 = -(*n1 + 1);
#line 230 "TB04AY.f"
		    ucoeff[i__ + (j + k * ucoeff_dim2) * ucoeff_dim1] += 
			    ddot_(&i__4, &at[nminl - k + 3 + at_dim1], &i__5, 
			    &dwork[ib + j - 1], &maxm);
#line 233 "TB04AY.f"
/* L120: */
#line 233 "TB04AY.f"
		}

#line 235 "TB04AY.f"
/* L130: */
#line 235 "TB04AY.f"
	    }

#line 237 "TB04AY.f"
	}
#line 238 "TB04AY.f"
/* L140: */
#line 238 "TB04AY.f"
    }

/*     Set optimal workspace dimension. */

#line 242 "TB04AY.f"
    dwork[1] = (doublereal) wrkopt;

#line 244 "TB04AY.f"
    return 0;
/* *** Last line of TB04AY *** */
} /* tb04ay_ */

