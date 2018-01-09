#line 1 "UD01CD.f"
/* UD01CD.f -- translated by f2c (version 20100827).
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

#line 1 "UD01CD.f"
/* Table of constant values */

static doublereal c_b5 = 0.;
static integer c__3 = 3;
static integer c__1 = 1;
static integer c__5 = 5;

/* Subroutine */ int ud01cd_(integer *mp, integer *np, integer *dp, integer *
	nin, doublereal *p, integer *ldp1, integer *ldp2, integer *info)
{
    /* System generated locals */
    integer p_dim1, p_dim2, p_offset, i__1;

    /* Builtin functions */
    integer s_rsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_rsle(void);

    /* Local variables */
    static integer d__, i__, j, k;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 0, 1, 0, 0 };
    static cilist io___6 = { 0, 0, 0, 0, 0 };
    static cilist io___7 = { 0, 0, 0, 0, 0 };



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

/*     To read the elements of a sparse matrix polynomial */
/*                                                    dp-1           dp */
/*        P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s    + P(dp) * s  . */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     MP      (input) INTEGER */
/*             The number of rows of the matrix polynomial P(s). */
/*             MP >= 1. */

/*     NP      (input) INTEGER */
/*             The number of columns of the matrix polynomial P(s). */
/*             NP >= 1. */

/*     DP      (input) INTEGER */
/*             The degree of the matrix polynomial P(s).  DP >= 0. */

/*     NIN     (input) INTEGER */
/*             The input channel from which the elements of P(s) are */
/*             read.  NIN >= 0. */

/*     P       (output) DOUBLE PRECISION array, dimension */
/*             (LDP1,LDP2,DP+1) */
/*             The leading MP-by-NP-by-(DP+1) part of this array contains */
/*             the coefficients of the matrix polynomial P(s). */
/*             Specifically, P(i,j,k) contains the coefficient of */
/*             s**(k-1) of the polynomial which is the (i,j)-th element */
/*             of P(s), where i = 1,2,...,MP, j = 1,2,...,NP and */
/*             k = 1,2,...,DP+1. */
/*             The not assigned elements are set to zero. */

/*     LDP1    INTEGER */
/*             The leading dimension of array P.  LDP1 >= MP. */

/*     LDP2    INTEGER */
/*             The second dimension of array P.  LDP2 >= NP. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1 : if a row index i is read with i < 1 or i > MP or */
/*                   a column index j is read with j < 1 or j > NP or */
/*                   a coefficient degree d is read with d < 0 or */
/*                   d > DP + 1. This is a warning. */

/*     METHOD */

/*     First, the elements P(i,j,k) with 1 <= i <= MP, 1 <= j <= NP and */
/*     1 <= k <= DP + 1 are set to zero. Next the nonzero (polynomial) */
/*     elements are read from the input file NIN. Each nonzero element is */
/*     given by the values i, j, d, P(i,j,k), k = 1, ..., d+1, where d is */
/*     the degree and P(i,j,k) is the coefficient of s**(k-1) in the */
/*     (i,j)-th element of P(s), i.e., let */
/*                                                              d */
/*         P   (s) = P   (0) + P   (1) * s + . . . + P   (d) * s */
/*          i,j       i,j       i,j                   i,j */

/*     be the nonzero (i,j)-th element of the matrix polynomial P(s). */

/*     Then P(i,j,k) corresponds to coefficient P   (k-1), k = 1,...,d+1. */
/*                                               i,j */
/*     For each nonzero element, the values i, j, and d are read as one */
/*     record of the file NIN, and the values P(i,j,k), k = 1,...,d+1, */
/*     are read as the following record. */
/*     The routine terminates after the last line has been read. */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */
/*     Based on routine RDSPOM by A.J. Geurts, Eindhoven University of */
/*     Technology, Holland. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

#line 128 "UD01CD.f"
    /* Parameter adjustments */
#line 128 "UD01CD.f"
    p_dim1 = *ldp1;
#line 128 "UD01CD.f"
    p_dim2 = *ldp2;
#line 128 "UD01CD.f"
    p_offset = 1 + p_dim1 * (1 + p_dim2);
#line 128 "UD01CD.f"
    p -= p_offset;
#line 128 "UD01CD.f"

#line 128 "UD01CD.f"
    /* Function Body */
#line 128 "UD01CD.f"
    *info = 0;

/*     Check the input scalar arguments. */

#line 132 "UD01CD.f"
    if (*mp < 1) {
#line 133 "UD01CD.f"
	*info = -1;
#line 134 "UD01CD.f"
    } else if (*np < 1) {
#line 135 "UD01CD.f"
	*info = -2;
#line 136 "UD01CD.f"
    } else if (*dp < 0) {
#line 137 "UD01CD.f"
	*info = -3;
#line 138 "UD01CD.f"
    } else if (*nin < 0) {
#line 139 "UD01CD.f"
	*info = -4;
#line 140 "UD01CD.f"
    } else if (*ldp1 < *mp) {
#line 141 "UD01CD.f"
	*info = -6;
#line 142 "UD01CD.f"
    } else if (*ldp2 < *np) {
#line 143 "UD01CD.f"
	*info = -7;
#line 144 "UD01CD.f"
    }

#line 146 "UD01CD.f"
    if (*info != 0) {

/*        Error return. */

#line 150 "UD01CD.f"
	i__1 = -(*info);
#line 150 "UD01CD.f"
	xerbla_("UD01CD", &i__1, (ftnlen)6);
#line 151 "UD01CD.f"
	return 0;
#line 152 "UD01CD.f"
    }

#line 154 "UD01CD.f"
    i__1 = *dp + 1;
#line 154 "UD01CD.f"
    for (k = 1; k <= i__1; ++k) {
#line 155 "UD01CD.f"
	dlaset_("Full", mp, np, &c_b5, &c_b5, &p[(k * p_dim2 + 1) * p_dim1 + 
		1], ldp1, (ftnlen)4);
#line 156 "UD01CD.f"
/* L10: */
#line 156 "UD01CD.f"
    }

/*     Read (i, j, d, P(i,j,k), k=1,...,d+1) of the nonzero elements one */
/*     by one. */

#line 161 "UD01CD.f"
L20:
#line 161 "UD01CD.f"
    io___2.ciunit = *nin;
#line 161 "UD01CD.f"
    i__1 = s_rsle(&io___2);
#line 161 "UD01CD.f"
    if (i__1 != 0) {
#line 161 "UD01CD.f"
	goto L30;
#line 161 "UD01CD.f"
    }
#line 161 "UD01CD.f"
    i__1 = do_lio(&c__3, &c__1, (char *)&i__, (ftnlen)sizeof(integer));
#line 161 "UD01CD.f"
    if (i__1 != 0) {
#line 161 "UD01CD.f"
	goto L30;
#line 161 "UD01CD.f"
    }
#line 161 "UD01CD.f"
    i__1 = do_lio(&c__3, &c__1, (char *)&j, (ftnlen)sizeof(integer));
#line 161 "UD01CD.f"
    if (i__1 != 0) {
#line 161 "UD01CD.f"
	goto L30;
#line 161 "UD01CD.f"
    }
#line 161 "UD01CD.f"
    i__1 = do_lio(&c__3, &c__1, (char *)&d__, (ftnlen)sizeof(integer));
#line 161 "UD01CD.f"
    if (i__1 != 0) {
#line 161 "UD01CD.f"
	goto L30;
#line 161 "UD01CD.f"
    }
#line 161 "UD01CD.f"
    i__1 = e_rsle();
#line 161 "UD01CD.f"
    if (i__1 != 0) {
#line 161 "UD01CD.f"
	goto L30;
#line 161 "UD01CD.f"
    }
#line 162 "UD01CD.f"
    if (i__ < 1 || i__ > *mp || j < 1 || j > *np || d__ < 0 || d__ > *dp + 1) 
	    {
#line 164 "UD01CD.f"
	*info = 1;
#line 165 "UD01CD.f"
	io___6.ciunit = *nin;
#line 165 "UD01CD.f"
	s_rsle(&io___6);
#line 165 "UD01CD.f"
	e_rsle();
#line 166 "UD01CD.f"
    } else {
#line 167 "UD01CD.f"
	io___7.ciunit = *nin;
#line 167 "UD01CD.f"
	s_rsle(&io___7);
#line 167 "UD01CD.f"
	i__1 = d__ + 1;
#line 167 "UD01CD.f"
	for (k = 1; k <= i__1; ++k) {
#line 167 "UD01CD.f"
	    do_lio(&c__5, &c__1, (char *)&p[i__ + (j + k * p_dim2) * p_dim1], 
		    (ftnlen)sizeof(doublereal));
#line 167 "UD01CD.f"
	}
#line 167 "UD01CD.f"
	e_rsle();
#line 168 "UD01CD.f"
    }
#line 169 "UD01CD.f"
    goto L20;

#line 171 "UD01CD.f"
L30:
#line 172 "UD01CD.f"
    return 0;
/* *** Last line of UD01CD *** */
} /* ud01cd_ */

