#line 1 "SB04QR.f"
/* SB04QR.f -- translated by f2c (version 20100827).
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

#line 1 "SB04QR.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int sb04qr_(integer *m, doublereal *d__, integer *ipr, 
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
/*     matrix has zeros below the third subdiagonal and zero elements on */
/*     the third subdiagonal with even column indices. The matrix is */
/*     stored compactly, row-wise. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the system.  M >= 0, M even. */
/*             Note that parameter M should have twice the value in the */
/*             original problem (see SLICOT Library routine SB04QU). */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*             (M*M/2+4*M) */
/*             On entry, the first M*M/2 + 3*M elements of this array */
/*             must contain the coefficient matrix, stored compactly, */
/*             row-wise, and the next M elements must contain the right */
/*             hand side of the linear system, as set by SLICOT Library */
/*             routine SB04QU. */
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

/*     [2] Sima, V. */
/*         Algorithms for Linear-quadratic Optimization. */
/*         Marcel Dekker, Inc., New York, 1996. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     D. Sima, University of Bucharest, May 2000. */

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

#line 112 "SB04QR.f"
    /* Parameter adjustments */
#line 112 "SB04QR.f"
    --ipr;
#line 112 "SB04QR.f"
    --d__;
#line 112 "SB04QR.f"

#line 112 "SB04QR.f"
    /* Function Body */
#line 112 "SB04QR.f"
    *info = 0;
#line 113 "SB04QR.f"
    i2 = *m * *m / 2 + *m * 3;
#line 114 "SB04QR.f"
    mpi = *m;
#line 115 "SB04QR.f"
    iprm = i2;
#line 116 "SB04QR.f"
    m1 = *m;
#line 117 "SB04QR.f"
    i1 = 1;

#line 119 "SB04QR.f"
    i__1 = *m;
#line 119 "SB04QR.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 120 "SB04QR.f"
	++mpi;
#line 121 "SB04QR.f"
	++iprm;
#line 122 "SB04QR.f"
	ipr[mpi] = i1;
#line 123 "SB04QR.f"
	ipr[i__] = iprm;
#line 124 "SB04QR.f"
	i1 += m1;
#line 125 "SB04QR.f"
	if (i__ >= 4 && i__ % 2 == 0) {
#line 125 "SB04QR.f"
	    m1 += -2;
#line 125 "SB04QR.f"
	}
#line 126 "SB04QR.f"
/* L20: */
#line 126 "SB04QR.f"
    }

#line 128 "SB04QR.f"
    m1 = *m - 1;
#line 129 "SB04QR.f"
    mpi1 = *m + 1;

/*     Reduce to upper triangular form. */

#line 133 "SB04QR.f"
    i__1 = m1;
#line 133 "SB04QR.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 134 "SB04QR.f"
	mpi = mpi1;
#line 135 "SB04QR.f"
	++mpi1;
#line 136 "SB04QR.f"
	iprm = ipr[mpi];
#line 137 "SB04QR.f"
	d1 = d__[iprm];
#line 138 "SB04QR.f"
	i1 = 3;
#line 139 "SB04QR.f"
	if (i__ % 2 == 0) {
#line 139 "SB04QR.f"
	    i1 = 2;
#line 139 "SB04QR.f"
	}
#line 140 "SB04QR.f"
	if (i__ == m1) {
#line 140 "SB04QR.f"
	    i1 = 1;
#line 140 "SB04QR.f"
	}
#line 141 "SB04QR.f"
	mpi2 = mpi + i1;
#line 142 "SB04QR.f"
	l = 0;
#line 143 "SB04QR.f"
	dmax__ = abs(d1);

#line 145 "SB04QR.f"
	i__2 = mpi2;
#line 145 "SB04QR.f"
	for (j = mpi1; j <= i__2; ++j) {
#line 146 "SB04QR.f"
	    d2 = d__[ipr[j]];
#line 147 "SB04QR.f"
	    d3 = abs(d2);
#line 148 "SB04QR.f"
	    if (d3 > dmax__) {
#line 149 "SB04QR.f"
		dmax__ = d3;
#line 150 "SB04QR.f"
		d1 = d2;
#line 151 "SB04QR.f"
		l = j - mpi;
#line 152 "SB04QR.f"
	    }
#line 153 "SB04QR.f"
/* L40: */
#line 153 "SB04QR.f"
	}

/*        Check singularity. */

#line 157 "SB04QR.f"
	if (dmax__ == 0.) {
#line 158 "SB04QR.f"
	    *info = 1;
#line 159 "SB04QR.f"
	    return 0;
#line 160 "SB04QR.f"
	}

#line 162 "SB04QR.f"
	if (l > 0) {

/*           Permute the row indices. */

#line 166 "SB04QR.f"
	    k = iprm;
#line 167 "SB04QR.f"
	    j = mpi + l;
#line 168 "SB04QR.f"
	    iprm = ipr[j];
#line 169 "SB04QR.f"
	    ipr[j] = k;
#line 170 "SB04QR.f"
	    ipr[mpi] = iprm;
#line 171 "SB04QR.f"
	    k = ipr[i__];
#line 172 "SB04QR.f"
	    i2 = i__ + l;
#line 173 "SB04QR.f"
	    ipr[i__] = ipr[i2];
#line 174 "SB04QR.f"
	    ipr[i2] = k;
#line 175 "SB04QR.f"
	}
#line 176 "SB04QR.f"
	++iprm;

/*        Annihilate the subdiagonal elements of the matrix. */

#line 180 "SB04QR.f"
	i2 = i__;
#line 181 "SB04QR.f"
	d3 = d__[ipr[i__]];

#line 183 "SB04QR.f"
	i__2 = mpi2;
#line 183 "SB04QR.f"
	for (j = mpi1; j <= i__2; ++j) {
#line 184 "SB04QR.f"
	    ++i2;
#line 185 "SB04QR.f"
	    iprm1 = ipr[j];
#line 186 "SB04QR.f"
	    dmax__ = -d__[iprm1] / d1;
#line 187 "SB04QR.f"
	    d__[ipr[i2]] += dmax__ * d3;
#line 188 "SB04QR.f"
	    i__3 = *m - i__;
#line 188 "SB04QR.f"
	    daxpy_(&i__3, &dmax__, &d__[iprm], &c__1, &d__[iprm1 + 1], &c__1);
#line 189 "SB04QR.f"
	    ++ipr[j];
#line 190 "SB04QR.f"
/* L60: */
#line 190 "SB04QR.f"
	}

#line 192 "SB04QR.f"
/* L80: */
#line 192 "SB04QR.f"
    }

#line 194 "SB04QR.f"
    mpi = *m + *m;
#line 195 "SB04QR.f"
    iprm = ipr[mpi];

/*     Check singularity. */

#line 199 "SB04QR.f"
    if (d__[iprm] == 0.) {
#line 200 "SB04QR.f"
	*info = 1;
#line 201 "SB04QR.f"
	return 0;
#line 202 "SB04QR.f"
    }

/*     Back substitution. */

#line 206 "SB04QR.f"
    d__[ipr[*m]] /= d__[iprm];

#line 208 "SB04QR.f"
    for (i__ = m1; i__ >= 1; --i__) {
#line 209 "SB04QR.f"
	--mpi;
#line 210 "SB04QR.f"
	iprm = ipr[mpi];
#line 211 "SB04QR.f"
	iprm1 = iprm;
#line 212 "SB04QR.f"
	dmax__ = 0.;

#line 214 "SB04QR.f"
	i__1 = *m;
#line 214 "SB04QR.f"
	for (k = i__ + 1; k <= i__1; ++k) {
#line 215 "SB04QR.f"
	    ++iprm1;
#line 216 "SB04QR.f"
	    dmax__ += d__[ipr[k]] * d__[iprm1];
#line 217 "SB04QR.f"
/* L100: */
#line 217 "SB04QR.f"
	}

#line 219 "SB04QR.f"
	d__[ipr[i__]] = (d__[ipr[i__]] - dmax__) / d__[iprm];
#line 220 "SB04QR.f"
/* L120: */
#line 220 "SB04QR.f"
    }

#line 222 "SB04QR.f"
    return 0;
/* *** Last line of SB04QR *** */
} /* sb04qr_ */

