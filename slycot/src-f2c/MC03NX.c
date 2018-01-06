#line 1 "MC03NX.f"
/* MC03NX.f -- translated by f2c (version 20100827).
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

#line 1 "MC03NX.f"
/* Table of constant values */

static doublereal c_b3 = 0.;
static doublereal c_b4 = 1.;
static doublereal c_b13 = -1.;
static integer c__1 = 1;

/* Subroutine */ int mc03nx_(integer *mp, integer *np, integer *dp, 
	doublereal *p, integer *ldp1, integer *ldp2, doublereal *a, integer *
	lda, doublereal *e, integer *lde)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, p_dim1, p_dim2, p_offset, 
	    i__1;

    /* Local variables */
    static integer j, k, h1, hb, he, hi;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), dlaset_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, ftnlen);


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

/*     Given an MP-by-NP polynomial matrix of degree dp */
/*                                    dp-1            dp */
/*     P(s) = P(0) + ... + P(dp-1) * s     + P(dp) * s            (1) */

/*     the routine composes the related pencil s*E-A where */

/*         | I              |           | O          -P(dp) | */
/*         |   .            |           | I .           .   | */
/*     A = |     .          |  and  E = |   . .         .   |.    (2) */
/*         |       .        |           |     . O       .   | */
/*         |         I      |           |       I  O -P(2)  | */
/*         |           P(0) |           |          I -P(1)  | */

/*     ================================================================== */
/*     REMARK: This routine is intended to be called only from the SLICOT */
/*             routine MC03ND. */
/*     ================================================================== */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     MP      (input) INTEGER */
/*             The number of rows of the polynomial matrix P(s). */
/*             MP >= 0. */

/*     NP      (input) INTEGER */
/*             The number of columns of the polynomial matrix P(s). */
/*             NP >= 0. */

/*     DP      (input) INTEGER */
/*             The degree of the polynomial matrix P(s).  DP >= 1. */

/*     P       (input) DOUBLE PRECISION array, dimension (LDP1,LDP2,DP+1) */
/*             The leading MP-by-NP-by-(DP+1) part of this array must */
/*             contain the coefficients of the polynomial matrix P(s) */
/*             in (1) in increasing powers of s. */

/*     LDP1    INTEGER */
/*             The leading dimension of array P.  LDP1 >= MAX(1,MP). */

/*     LDP2    INTEGER */
/*             The second dimension of array P.   LDP2 >= MAX(1,NP). */

/*     A       (output) DOUBLE PRECISION array, dimension */
/*             (LDA,(DP-1)*MP+NP) */
/*             The leading DP*MP-by-((DP-1)*MP+NP) part of this array */
/*             contains the matrix A as described in (2). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,DP*MP). */

/*     E       (output) DOUBLE PRECISION array, dimension */
/*             (LDE,(DP-1)*MP+NP) */
/*             The leading DP*MP-by-((DP-1)*MP+NP) part of this array */
/*             contains the matrix E as described in (2). */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,DP*MP). */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC03BX by G.J.H.H. van den Hurk. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, input output description, */
/*     polynomial matrix, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

#line 116 "MC03NX.f"
    /* Parameter adjustments */
#line 116 "MC03NX.f"
    p_dim1 = *ldp1;
#line 116 "MC03NX.f"
    p_dim2 = *ldp2;
#line 116 "MC03NX.f"
    p_offset = 1 + p_dim1 * (1 + p_dim2);
#line 116 "MC03NX.f"
    p -= p_offset;
#line 116 "MC03NX.f"
    a_dim1 = *lda;
#line 116 "MC03NX.f"
    a_offset = 1 + a_dim1;
#line 116 "MC03NX.f"
    a -= a_offset;
#line 116 "MC03NX.f"
    e_dim1 = *lde;
#line 116 "MC03NX.f"
    e_offset = 1 + e_dim1;
#line 116 "MC03NX.f"
    e -= e_offset;
#line 116 "MC03NX.f"

#line 116 "MC03NX.f"
    /* Function Body */
#line 116 "MC03NX.f"
    if (*mp <= 0 || *np <= 0) {
#line 116 "MC03NX.f"
	return 0;
#line 116 "MC03NX.f"
    }

/*     Initialisation of matrices A and E. */

#line 121 "MC03NX.f"
    h1 = *dp * *mp;
#line 122 "MC03NX.f"
    hb = h1 - *mp;
#line 123 "MC03NX.f"
    he = hb + *np;
#line 124 "MC03NX.f"
    dlaset_("Full", &h1, &he, &c_b3, &c_b4, &a[a_offset], lda, (ftnlen)4);
#line 125 "MC03NX.f"
    dlaset_("Full", mp, &hb, &c_b3, &c_b3, &e[e_offset], lde, (ftnlen)4);
#line 126 "MC03NX.f"
    dlacpy_("Full", &hb, &hb, &a[a_offset], lda, &e[*mp + 1 + e_dim1], lde, (
	    ftnlen)4);

/*     Insert the matrices P(0), P(1), ..., P(dp) at the right places */
/*     in the matrices A and E. */

#line 131 "MC03NX.f"
    ++hb;
#line 132 "MC03NX.f"
    dlacpy_("Full", mp, np, &p[(p_dim2 + 1) * p_dim1 + 1], ldp1, &a[hb + hb * 
	    a_dim1], lda, (ftnlen)4);
#line 133 "MC03NX.f"
    hi = 1;

#line 135 "MC03NX.f"
    for (k = *dp + 1; k >= 2; --k) {
#line 136 "MC03NX.f"
	dlacpy_("Full", mp, np, &p[(k * p_dim2 + 1) * p_dim1 + 1], ldp1, &e[
		hi + hb * e_dim1], lde, (ftnlen)4);
#line 137 "MC03NX.f"
	hi += *mp;
#line 138 "MC03NX.f"
/* L20: */
#line 138 "MC03NX.f"
    }

#line 140 "MC03NX.f"
    i__1 = he;
#line 140 "MC03NX.f"
    for (j = hb; j <= i__1; ++j) {
#line 141 "MC03NX.f"
	dscal_(&h1, &c_b13, &e[j * e_dim1 + 1], &c__1);
#line 142 "MC03NX.f"
/* L40: */
#line 142 "MC03NX.f"
    }

#line 144 "MC03NX.f"
    return 0;
/* *** Last line of MC03NX *** */
} /* mc03nx_ */

