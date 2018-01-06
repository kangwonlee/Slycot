#line 1 "SB04MR.f"
/* SB04MR.f -- translated by f2c (version 20100827).
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

#line 1 "SB04MR.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int sb04mr_(integer *m, doublereal *d__, integer *ipr, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal d1, d2, d3;
    static integer i1, i2, m1, mpi, mpi1, mpi2;
    static doublereal dmax__;
    static integer iprm, iprm1;
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
/*     matrix has zeros below the second subdiagonal. The matrix is */
/*     stored compactly, row-wise. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the system.  M >= 0. */
/*             Note that parameter M should have twice the value in the */
/*             original problem (see SLICOT Library routine SB04MU). */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*             (M*(M+1)/2+3*M) */
/*             On entry, the first M*(M+1)/2 + 2*M elements of this array */
/*             must contain the coefficient matrix, stored compactly, */
/*             row-wise, and the next M elements must contain the right */
/*             hand side of the linear system, as set by SLICOT Library */
/*             routine SB04MU. */
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
/*     Supersedes Release 2.0 routine SB04AR by G. Golub, S. Nash, and */
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

#line 110 "SB04MR.f"
    /* Parameter adjustments */
#line 110 "SB04MR.f"
    --ipr;
#line 110 "SB04MR.f"
    --d__;
#line 110 "SB04MR.f"

#line 110 "SB04MR.f"
    /* Function Body */
#line 110 "SB04MR.f"
    *info = 0;
#line 111 "SB04MR.f"
    i2 = *m * (*m + 5) / 2;
#line 112 "SB04MR.f"
    mpi = *m;
#line 113 "SB04MR.f"
    iprm = i2;
#line 114 "SB04MR.f"
    m1 = *m;
#line 115 "SB04MR.f"
    i1 = 1;

#line 117 "SB04MR.f"
    i__1 = *m;
#line 117 "SB04MR.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 118 "SB04MR.f"
	++mpi;
#line 119 "SB04MR.f"
	++iprm;
#line 120 "SB04MR.f"
	ipr[mpi] = i1;
#line 121 "SB04MR.f"
	ipr[i__] = iprm;
#line 122 "SB04MR.f"
	i1 += m1;
#line 123 "SB04MR.f"
	if (i__ >= 3) {
#line 123 "SB04MR.f"
	    --m1;
#line 123 "SB04MR.f"
	}
#line 124 "SB04MR.f"
/* L20: */
#line 124 "SB04MR.f"
    }

#line 126 "SB04MR.f"
    m1 = *m - 1;
#line 127 "SB04MR.f"
    mpi1 = *m + 1;

/*     Reduce to upper triangular form. */

#line 131 "SB04MR.f"
    i__1 = m1;
#line 131 "SB04MR.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 132 "SB04MR.f"
	mpi = mpi1;
#line 133 "SB04MR.f"
	++mpi1;
#line 134 "SB04MR.f"
	iprm = ipr[mpi];
#line 135 "SB04MR.f"
	d1 = d__[iprm];
#line 136 "SB04MR.f"
	i1 = 2;
#line 137 "SB04MR.f"
	if (i__ == m1) {
#line 137 "SB04MR.f"
	    i1 = 1;
#line 137 "SB04MR.f"
	}
#line 138 "SB04MR.f"
	mpi2 = mpi + i1;
#line 139 "SB04MR.f"
	l = 0;
#line 140 "SB04MR.f"
	dmax__ = abs(d1);

#line 142 "SB04MR.f"
	i__2 = mpi2;
#line 142 "SB04MR.f"
	for (j = mpi1; j <= i__2; ++j) {
#line 143 "SB04MR.f"
	    d2 = d__[ipr[j]];
#line 144 "SB04MR.f"
	    d3 = abs(d2);
#line 145 "SB04MR.f"
	    if (d3 > dmax__) {
#line 146 "SB04MR.f"
		dmax__ = d3;
#line 147 "SB04MR.f"
		d1 = d2;
#line 148 "SB04MR.f"
		l = j - mpi;
#line 149 "SB04MR.f"
	    }
#line 150 "SB04MR.f"
/* L40: */
#line 150 "SB04MR.f"
	}

/*        Check singularity. */

#line 154 "SB04MR.f"
	if (dmax__ == 0.) {
#line 155 "SB04MR.f"
	    *info = 1;
#line 156 "SB04MR.f"
	    return 0;
#line 157 "SB04MR.f"
	}

#line 159 "SB04MR.f"
	if (l > 0) {

/*           Permute the row indices. */

#line 163 "SB04MR.f"
	    k = iprm;
#line 164 "SB04MR.f"
	    j = mpi + l;
#line 165 "SB04MR.f"
	    iprm = ipr[j];
#line 166 "SB04MR.f"
	    ipr[j] = k;
#line 167 "SB04MR.f"
	    ipr[mpi] = iprm;
#line 168 "SB04MR.f"
	    k = ipr[i__];
#line 169 "SB04MR.f"
	    i2 = i__ + l;
#line 170 "SB04MR.f"
	    ipr[i__] = ipr[i2];
#line 171 "SB04MR.f"
	    ipr[i2] = k;
#line 172 "SB04MR.f"
	}
#line 173 "SB04MR.f"
	++iprm;

/*        Annihilate the subdiagonal elements of the matrix. */

#line 177 "SB04MR.f"
	i2 = i__;
#line 178 "SB04MR.f"
	d3 = d__[ipr[i__]];

#line 180 "SB04MR.f"
	i__2 = mpi2;
#line 180 "SB04MR.f"
	for (j = mpi1; j <= i__2; ++j) {
#line 181 "SB04MR.f"
	    ++i2;
#line 182 "SB04MR.f"
	    iprm1 = ipr[j];
#line 183 "SB04MR.f"
	    dmax__ = -d__[iprm1] / d1;
#line 184 "SB04MR.f"
	    d__[ipr[i2]] += dmax__ * d3;
#line 185 "SB04MR.f"
	    i__3 = *m - i__;
#line 185 "SB04MR.f"
	    daxpy_(&i__3, &dmax__, &d__[iprm], &c__1, &d__[iprm1 + 1], &c__1);
#line 186 "SB04MR.f"
/* L60: */
#line 186 "SB04MR.f"
	}

#line 188 "SB04MR.f"
	++ipr[mpi1];
#line 189 "SB04MR.f"
	if (i__ != m1) {
#line 189 "SB04MR.f"
	    ++ipr[mpi2];
#line 189 "SB04MR.f"
	}
#line 190 "SB04MR.f"
/* L80: */
#line 190 "SB04MR.f"
    }

#line 192 "SB04MR.f"
    mpi = *m + *m;
#line 193 "SB04MR.f"
    iprm = ipr[mpi];

/*     Check singularity. */

#line 197 "SB04MR.f"
    if (d__[iprm] == 0.) {
#line 198 "SB04MR.f"
	*info = 1;
#line 199 "SB04MR.f"
	return 0;
#line 200 "SB04MR.f"
    }

/*     Back substitution. */

#line 204 "SB04MR.f"
    d__[ipr[*m]] /= d__[iprm];

#line 206 "SB04MR.f"
    for (i__ = m1; i__ >= 1; --i__) {
#line 207 "SB04MR.f"
	--mpi;
#line 208 "SB04MR.f"
	iprm = ipr[mpi];
#line 209 "SB04MR.f"
	iprm1 = iprm;
#line 210 "SB04MR.f"
	dmax__ = 0.;

#line 212 "SB04MR.f"
	i__1 = *m;
#line 212 "SB04MR.f"
	for (k = i__ + 1; k <= i__1; ++k) {
#line 213 "SB04MR.f"
	    ++iprm1;
#line 214 "SB04MR.f"
	    dmax__ += d__[ipr[k]] * d__[iprm1];
#line 215 "SB04MR.f"
/* L100: */
#line 215 "SB04MR.f"
	}

#line 217 "SB04MR.f"
	d__[ipr[i__]] = (d__[ipr[i__]] - dmax__) / d__[iprm];
#line 218 "SB04MR.f"
/* L120: */
#line 218 "SB04MR.f"
    }

#line 220 "SB04MR.f"
    return 0;
/* *** Last line of SB04MR *** */
} /* sb04mr_ */

