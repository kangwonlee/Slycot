#line 1 "SB04MW.f"
/* SB04MW.f -- translated by f2c (version 20100827).
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

#line 1 "SB04MW.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int sb04mw_(integer *m, doublereal *d__, integer *ipr, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k;
    static doublereal d1, d2;
    static integer i1, m1, m2, mpi, iprm;
    static doublereal mult;
    static integer iprm1;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


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

/*     To solve a linear algebraic system of order M whose coefficient */
/*     matrix is in upper Hessenberg form, stored compactly, row-wise. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the system.  M >= 0. */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*             (M*(M+1)/2+2*M) */
/*             On entry, the first M*(M+1)/2 + M elements of this array */
/*             must contain an upper Hessenberg matrix, stored compactly, */
/*             row-wise, and the next M elements must contain the right */
/*             hand side of the linear system, as set by SLICOT Library */
/*             routine SB04MY. */
/*             On exit, the content of this array is updated, the last M */
/*             elements containing the solution with components */
/*             interchanged (see IPR). */

/*     IPR     (output) INTEGER array, dimension (2*M) */
/*             The leading M elements contain information about the */
/*             row interchanges performed for solving the system. */
/*             Specifically, the i-th component of the solution is */
/*             specified by IPR(i). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if a singular matrix was encountered. */

/*     METHOD */

/*     Gaussian elimination with partial pivoting is used. The rows of */
/*     the matrix are not actually permuted, only their indices are */
/*     interchanged in array IPR. */

/*     REFERENCES */

/*     [1] Golub, G.H., Nash, S. and Van Loan, C.F. */
/*         A Hessenberg-Schur method for the problem AX + XB = C. */
/*         IEEE Trans. Auto. Contr., AC-24, pp. 909-913, 1979. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Sep. 1997. */
/*     Supersedes Release 2.0 routine SB04AW by G. Golub, S. Nash, and */
/*     C. Van Loan, Stanford University, California, United States of */
/*     America, January 1982. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Hessenberg form, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 106 "SB04MW.f"
    /* Parameter adjustments */
#line 106 "SB04MW.f"
    --ipr;
#line 106 "SB04MW.f"
    --d__;
#line 106 "SB04MW.f"

#line 106 "SB04MW.f"
    /* Function Body */
#line 106 "SB04MW.f"
    *info = 0;
#line 107 "SB04MW.f"
    m1 = *m * (*m + 3) / 2;
#line 108 "SB04MW.f"
    m2 = *m + *m;
#line 109 "SB04MW.f"
    mpi = *m;
#line 110 "SB04MW.f"
    iprm = m1;
#line 111 "SB04MW.f"
    m1 = *m;
#line 112 "SB04MW.f"
    i1 = 1;

#line 114 "SB04MW.f"
    i__1 = *m;
#line 114 "SB04MW.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 115 "SB04MW.f"
	++mpi;
#line 116 "SB04MW.f"
	++iprm;
#line 117 "SB04MW.f"
	ipr[mpi] = i1;
#line 118 "SB04MW.f"
	ipr[i__] = iprm;
#line 119 "SB04MW.f"
	i1 += m1;
#line 120 "SB04MW.f"
	if (i__ > 1) {
#line 120 "SB04MW.f"
	    --m1;
#line 120 "SB04MW.f"
	}
#line 121 "SB04MW.f"
/* L20: */
#line 121 "SB04MW.f"
    }

#line 123 "SB04MW.f"
    m1 = *m - 1;
#line 124 "SB04MW.f"
    mpi = *m;

/*     Reduce to upper triangular form. */

#line 128 "SB04MW.f"
    i__1 = m1;
#line 128 "SB04MW.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 129 "SB04MW.f"
	i1 = i__ + 1;
#line 130 "SB04MW.f"
	++mpi;
#line 131 "SB04MW.f"
	iprm = ipr[mpi];
#line 132 "SB04MW.f"
	iprm1 = ipr[mpi + 1];
#line 133 "SB04MW.f"
	d1 = d__[iprm];
#line 134 "SB04MW.f"
	d2 = d__[iprm1];
#line 135 "SB04MW.f"
	if (abs(d1) <= abs(d2)) {

/*           Permute the row indices. */

#line 139 "SB04MW.f"
	    k = iprm;
#line 140 "SB04MW.f"
	    ipr[mpi] = iprm1;
#line 141 "SB04MW.f"
	    iprm = iprm1;
#line 142 "SB04MW.f"
	    iprm1 = k;
#line 143 "SB04MW.f"
	    k = ipr[i__];
#line 144 "SB04MW.f"
	    ipr[i__] = ipr[i1];
#line 145 "SB04MW.f"
	    ipr[i1] = k;
#line 146 "SB04MW.f"
	    d1 = d2;
#line 147 "SB04MW.f"
	}

/*        Check singularity. */

#line 151 "SB04MW.f"
	if (d1 == 0.) {
#line 152 "SB04MW.f"
	    *info = 1;
#line 153 "SB04MW.f"
	    return 0;
#line 154 "SB04MW.f"
	}

#line 156 "SB04MW.f"
	mult = -d__[iprm1] / d1;
#line 157 "SB04MW.f"
	++iprm1;
#line 158 "SB04MW.f"
	ipr[mpi + 1] = iprm1;

/*        Annihilate the subdiagonal elements of the matrix. */

#line 162 "SB04MW.f"
	d__[ipr[i1]] += mult * d__[ipr[i__]];
#line 163 "SB04MW.f"
	i__2 = *m - i__;
#line 163 "SB04MW.f"
	daxpy_(&i__2, &mult, &d__[iprm + 1], &c__1, &d__[iprm1], &c__1);
#line 164 "SB04MW.f"
/* L40: */
#line 164 "SB04MW.f"
    }

/*     Check singularity. */

#line 168 "SB04MW.f"
    if (d__[ipr[m2]] == 0.) {
#line 169 "SB04MW.f"
	*info = 1;
#line 170 "SB04MW.f"
	return 0;
#line 171 "SB04MW.f"
    }

/*     Back substitution. */

#line 175 "SB04MW.f"
    d__[ipr[*m]] /= d__[ipr[m2]];
#line 176 "SB04MW.f"
    mpi = m2;

#line 178 "SB04MW.f"
    for (i__ = m1; i__ >= 1; --i__) {
#line 179 "SB04MW.f"
	--mpi;
#line 180 "SB04MW.f"
	iprm = ipr[mpi];
#line 181 "SB04MW.f"
	iprm1 = iprm;
#line 182 "SB04MW.f"
	mult = 0.;

#line 184 "SB04MW.f"
	i__1 = *m;
#line 184 "SB04MW.f"
	for (i1 = i__ + 1; i1 <= i__1; ++i1) {
#line 185 "SB04MW.f"
	    ++iprm1;
#line 186 "SB04MW.f"
	    mult += d__[ipr[i1]] * d__[iprm1];
#line 187 "SB04MW.f"
/* L60: */
#line 187 "SB04MW.f"
	}

#line 189 "SB04MW.f"
	d__[ipr[i__]] = (d__[ipr[i__]] - mult) / d__[iprm];
#line 190 "SB04MW.f"
/* L80: */
#line 190 "SB04MW.f"
    }

#line 192 "SB04MW.f"
    return 0;
/* *** Last line of SB04MW *** */
} /* sb04mw_ */

