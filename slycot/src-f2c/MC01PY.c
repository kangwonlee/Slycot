#line 1 "MC01PY.f"
/* MC01PY.f -- translated by f2c (version 20100827).
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

#line 1 "MC01PY.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mc01py_(integer *k, doublereal *rez, doublereal *imz, 
	doublereal *p, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__;
    static doublereal u, v;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *), xerbla_(char *,
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

/*     To compute the coefficients of a real polynomial P(x) from its */
/*     zeros. The coefficients are stored in decreasing order of the */
/*     powers of x. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     K       (input) INTEGER */
/*             The number of zeros (and hence the degree) of P(x). */
/*             K >= 0. */

/*     REZ     (input) DOUBLE PRECISION array, dimension (K) */
/*     IMZ     (input) DOUBLE PRECISION array, dimension (K) */
/*             The real and imaginary parts of the i-th zero of P(x) */
/*             must be stored in REZ(i) and IMZ(i), respectively, where */
/*             i = 1, 2, ..., K. The zeros may be supplied in any order, */
/*             except that complex conjugate zeros must appear */
/*             consecutively. */

/*     P       (output) DOUBLE PRECISION array, dimension (K+1) */
/*             This array contains the coefficients of P(x) in decreasing */
/*             powers of x. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (K) */
/*             If K = 0, this array is not referenced. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, (REZ(i),IMZ(i)) is a complex zero but */
/*                   (REZ(i-1),IMZ(i-1)) is not its conjugate. */

/*     METHOD */

/*     The routine computes the coefficients of the real K-th degree */
/*     polynomial P(x) as */

/*        P(x) = (x - r(1)) * (x - r(2)) * ... * (x - r(K)) */

/*     where r(i) = (REZ(i),IMZ(i)). */

/*     Note that REZ(i) = REZ(j) and IMZ(i) = -IMZ(j) if r(i) and r(j) */
/*     form a complex conjugate pair (where i <> j), and that IMZ(i) = 0 */
/*     if r(i) is real. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2002. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, polynomial operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 108 "MC01PY.f"
    /* Parameter adjustments */
#line 108 "MC01PY.f"
    --dwork;
#line 108 "MC01PY.f"
    --p;
#line 108 "MC01PY.f"
    --imz;
#line 108 "MC01PY.f"
    --rez;
#line 108 "MC01PY.f"

#line 108 "MC01PY.f"
    /* Function Body */
#line 108 "MC01PY.f"
    if (*k < 0) {
#line 109 "MC01PY.f"
	*info = -1;

/*        Error return. */

#line 113 "MC01PY.f"
	i__1 = -(*info);
#line 113 "MC01PY.f"
	xerbla_("MC01PY", &i__1, (ftnlen)6);
#line 114 "MC01PY.f"
	return 0;
#line 115 "MC01PY.f"
    }

/*     Quick return if possible. */

#line 119 "MC01PY.f"
    *info = 0;
#line 120 "MC01PY.f"
    p[1] = 1.;
#line 121 "MC01PY.f"
    if (*k == 0) {
#line 121 "MC01PY.f"
	return 0;
#line 121 "MC01PY.f"
    }

#line 124 "MC01PY.f"
    i__ = 1;
/*     WHILE ( I <= K ) DO */
#line 126 "MC01PY.f"
L20:
#line 126 "MC01PY.f"
    if (i__ <= *k) {
#line 127 "MC01PY.f"
	u = rez[i__];
#line 128 "MC01PY.f"
	v = imz[i__];
#line 129 "MC01PY.f"
	dwork[i__] = 0.;

#line 131 "MC01PY.f"
	if (v == 0.) {
#line 132 "MC01PY.f"
	    d__1 = -u;
#line 132 "MC01PY.f"
	    daxpy_(&i__, &d__1, &p[1], &c__1, &dwork[1], &c__1);

#line 134 "MC01PY.f"
	} else {
#line 135 "MC01PY.f"
	    if (i__ == *k) {
#line 136 "MC01PY.f"
		*info = *k;
#line 137 "MC01PY.f"
		return 0;
#line 138 "MC01PY.f"
	    } else if (u != rez[i__ + 1] || v != -imz[i__ + 1]) {
#line 139 "MC01PY.f"
		*info = i__ + 1;
#line 140 "MC01PY.f"
		return 0;
#line 141 "MC01PY.f"
	    }

#line 143 "MC01PY.f"
	    dwork[i__ + 1] = 0.;
#line 144 "MC01PY.f"
	    d__1 = -(u + u);
#line 144 "MC01PY.f"
	    daxpy_(&i__, &d__1, &p[1], &c__1, &dwork[1], &c__1);
/* Computing 2nd power */
#line 145 "MC01PY.f"
	    d__2 = u;
/* Computing 2nd power */
#line 145 "MC01PY.f"
	    d__3 = v;
#line 145 "MC01PY.f"
	    d__1 = d__2 * d__2 + d__3 * d__3;
#line 145 "MC01PY.f"
	    daxpy_(&i__, &d__1, &p[1], &c__1, &dwork[2], &c__1);
#line 146 "MC01PY.f"
	    ++i__;
#line 147 "MC01PY.f"
	}

#line 149 "MC01PY.f"
	dcopy_(&i__, &dwork[1], &c__1, &p[2], &c__1);
#line 150 "MC01PY.f"
	++i__;
#line 151 "MC01PY.f"
	goto L20;
#line 152 "MC01PY.f"
    }
/*     END WHILE 20 */

#line 155 "MC01PY.f"
    return 0;
/* *** Last line of MC01PY *** */
} /* mc01py_ */

