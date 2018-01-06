#line 1 "TB03AY.f"
/* TB03AY.f -- translated by f2c (version 20100827).
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

#line 1 "TB03AY.f"
/* Table of constant values */

static doublereal c_b6 = 1.;
static doublereal c_b7 = 0.;
static doublereal c_b10 = -1.;
static integer c__1 = 1;

/* Subroutine */ int tb03ay_(integer *nr, doublereal *a, integer *lda, 
	integer *indblk, integer *nblk, doublereal *vcoeff, integer *ldvco1, 
	integer *ldvco2, doublereal *pcoeff, integer *ldpco1, integer *ldpco2,
	 integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, pcoeff_dim1, pcoeff_dim2, pcoeff_offset, 
	    vcoeff_dim1, vcoeff_dim2, vcoeff_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l, ioff, joff, ncol, nrow;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemm_(char *, char *, integer *, integer *, integer *
	    , doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), daxpy_(
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *), dtrsm_(char *, char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer kplus, lwork, lstop;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer lstart, inplus;


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

/*     To calculate the (PWORK-by-NR) polynomial matrix V(s) one */
/*     (PWORK-by-NBLK(L-1)) block V:L-1(s) at a time, in reverse order */
/*     (L = INDBLK,...,1).  At each stage, the (NBLK(L)-by-NBLK(L)) poly- */
/*     nomial matrix W(s) = V2(s) * A2 is formed, where V2(s) is that */
/*     part of V(s) already computed and A2 is the subdiagonal (incl.) */
/*     part of the L-th column block of A; W(s) is temporarily stored in */
/*     the top left part of P(s), as is subsequently the further matrix */
/*     Wbar(s) = s * V:L(s) - W(s).  Then, except for the final stage */
/*     L = 1 (when the next step is to calculate P(s) itself, not here), */
/*     the top left part of V:L-1(s) is given by Wbar(s) * inv(R), where */
/*     R is the upper triangular part of the L-th superdiagonal block of */
/*     A.  Finally, note that the coefficient matrices W(.,.,K) can only */
/*     be non-zero for K = L + 1,...,INPLUS, with each of these matrices */
/*     having only its first NBLK(L-1) rows non-trivial.  Similarly, */
/*     Wbar(.,.,K) (and so clearly V:L-1(.,.,K) ) can only be non-zero */
/*     for K = L,...,INPLUS, with each of these having only its first */
/*     NBLK(K-1) rows non-trivial except for K = L, which has NBLK(L) */
/*     such rows. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Coprime matrix fraction, elementary polynomial operations, */
/*     polynomial matrix, state-space representation, transfer matrix. */

/*     NOTE: In the interests of speed, this routine does not check the */
/*           inputs for errors. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

#line 74 "TB03AY.f"
    /* Parameter adjustments */
#line 74 "TB03AY.f"
    a_dim1 = *lda;
#line 74 "TB03AY.f"
    a_offset = 1 + a_dim1;
#line 74 "TB03AY.f"
    a -= a_offset;
#line 74 "TB03AY.f"
    --nblk;
#line 74 "TB03AY.f"
    vcoeff_dim1 = *ldvco1;
#line 74 "TB03AY.f"
    vcoeff_dim2 = *ldvco2;
#line 74 "TB03AY.f"
    vcoeff_offset = 1 + vcoeff_dim1 * (1 + vcoeff_dim2);
#line 74 "TB03AY.f"
    vcoeff -= vcoeff_offset;
#line 74 "TB03AY.f"
    pcoeff_dim1 = *ldpco1;
#line 74 "TB03AY.f"
    pcoeff_dim2 = *ldpco2;
#line 74 "TB03AY.f"
    pcoeff_offset = 1 + pcoeff_dim1 * (1 + pcoeff_dim2);
#line 74 "TB03AY.f"
    pcoeff -= pcoeff_offset;
#line 74 "TB03AY.f"

#line 74 "TB03AY.f"
    /* Function Body */
#line 74 "TB03AY.f"
    *info = 0;
#line 75 "TB03AY.f"
    inplus = *indblk + 1;
#line 76 "TB03AY.f"
    joff = *nr;

/*     Calculate each column block V:LWORK-1(s) of V(s) in turn. */

#line 80 "TB03AY.f"
    i__1 = *indblk;
#line 80 "TB03AY.f"
    for (l = 1; l <= i__1; ++l) {
#line 81 "TB03AY.f"
	lwork = inplus - l;

/*        Determine number of columns of V:LWORK(s) & its position in V. */

#line 85 "TB03AY.f"
	ncol = nblk[lwork];
#line 86 "TB03AY.f"
	joff -= ncol;

/*        Find limits for V2(s) * A2 calculation: skips zero rows */
/*        in V(s). */

#line 91 "TB03AY.f"
	lstart = joff + 1;
#line 92 "TB03AY.f"
	lstop = joff;

/*        Calculate W(s) and store (temporarily) in top left part */
/*        of P(s). */

#line 97 "TB03AY.f"
	i__2 = inplus;
#line 97 "TB03AY.f"
	for (k = lwork + 1; k <= i__2; ++k) {
#line 98 "TB03AY.f"
	    nrow = nblk[k - 1];
#line 99 "TB03AY.f"
	    lstop += nrow;
#line 100 "TB03AY.f"
	    i__3 = lstop - lstart + 1;
#line 100 "TB03AY.f"
	    dgemm_("No transpose", "No transpose", &nrow, &ncol, &i__3, &c_b6,
		     &vcoeff[(lstart + k * vcoeff_dim2) * vcoeff_dim1 + 1], 
		    ldvco1, &a[lstart + (joff + 1) * a_dim1], lda, &c_b7, &
		    pcoeff[(k * pcoeff_dim2 + 1) * pcoeff_dim1 + 1], ldpco1, (
		    ftnlen)12, (ftnlen)12);
#line 104 "TB03AY.f"
/* L10: */
#line 104 "TB03AY.f"
	}

/*        Replace W(s) by Wbar(s) = s * V:L(s) - W(s). */

#line 108 "TB03AY.f"
	nrow = ncol;

#line 110 "TB03AY.f"
	i__2 = *indblk;
#line 110 "TB03AY.f"
	for (k = lwork; k <= i__2; ++k) {
#line 111 "TB03AY.f"
	    kplus = k + 1;

#line 113 "TB03AY.f"
	    i__3 = ncol;
#line 113 "TB03AY.f"
	    for (j = 1; j <= i__3; ++j) {
#line 114 "TB03AY.f"
		dscal_(&nrow, &c_b10, &pcoeff[(j + k * pcoeff_dim2) * 
			pcoeff_dim1 + 1], &c__1);
#line 115 "TB03AY.f"
		daxpy_(&nrow, &c_b6, &vcoeff[(joff + j + kplus * vcoeff_dim2) 
			* vcoeff_dim1 + 1], &c__1, &pcoeff[(j + k * 
			pcoeff_dim2) * pcoeff_dim1 + 1], &c__1);
#line 117 "TB03AY.f"
/* L20: */
#line 117 "TB03AY.f"
	    }

#line 119 "TB03AY.f"
	    nrow = nblk[k];
#line 120 "TB03AY.f"
/* L30: */
#line 120 "TB03AY.f"
	}

#line 122 "TB03AY.f"
	i__2 = ncol;
#line 122 "TB03AY.f"
	for (j = 1; j <= i__2; ++j) {
#line 123 "TB03AY.f"
	    dscal_(&nrow, &c_b10, &pcoeff[(j + inplus * pcoeff_dim2) * 
		    pcoeff_dim1 + 1], &c__1);
#line 124 "TB03AY.f"
/* L40: */
#line 124 "TB03AY.f"
	}

#line 126 "TB03AY.f"
	if (lwork != 1) {

/*           If not final stage, use the upper triangular R (from A) */
/*           to calculate V:L-1(s), finally storing this new block. */

#line 131 "TB03AY.f"
	    ioff = joff - nblk[lwork - 1];

#line 133 "TB03AY.f"
	    i__2 = ncol;
#line 133 "TB03AY.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 134 "TB03AY.f"
		if (a[ioff + i__ + (joff + i__) * a_dim1] == 0.) {

/*                 Error return. */

#line 138 "TB03AY.f"
		    *info = i__;
#line 139 "TB03AY.f"
		    return 0;
#line 140 "TB03AY.f"
		}
#line 141 "TB03AY.f"
/* L50: */
#line 141 "TB03AY.f"
	    }

#line 143 "TB03AY.f"
	    nrow = nblk[lwork];

#line 145 "TB03AY.f"
	    i__2 = inplus;
#line 145 "TB03AY.f"
	    for (k = lwork; k <= i__2; ++k) {
#line 146 "TB03AY.f"
		dlacpy_("Full", &nrow, &ncol, &pcoeff[(k * pcoeff_dim2 + 1) * 
			pcoeff_dim1 + 1], ldpco1, &vcoeff[(ioff + 1 + k * 
			vcoeff_dim2) * vcoeff_dim1 + 1], ldvco1, (ftnlen)4);
#line 148 "TB03AY.f"
		dtrsm_("Right", "Upper", "No Transpose", "Non-unit", &nrow, &
			ncol, &c_b6, &a[ioff + 1 + (joff + 1) * a_dim1], lda, 
			&vcoeff[(ioff + 1 + k * vcoeff_dim2) * vcoeff_dim1 + 
			1], ldvco1, (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)
			8);
#line 151 "TB03AY.f"
		nrow = nblk[k];
#line 152 "TB03AY.f"
/* L60: */
#line 152 "TB03AY.f"
	    }

#line 154 "TB03AY.f"
	}
#line 155 "TB03AY.f"
/* L70: */
#line 155 "TB03AY.f"
    }

#line 157 "TB03AY.f"
    return 0;
/* *** Last line of TB03AY *** */
} /* tb03ay_ */

