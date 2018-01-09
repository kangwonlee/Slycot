#line 1 "MC01OD.f"
/* MC01OD.f -- translated by f2c (version 20100827).
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

#line 1 "MC01OD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mc01od_(integer *k, doublereal *rez, doublereal *imz, 
	doublereal *rep, doublereal *imp, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static doublereal u, v;
    static integer k2;
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

/*     To compute the coefficients of a complex polynomial P(x) from its */
/*     zeros. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     K       (input) INTEGER */
/*             The number of zeros (and hence the degree) of P(x). */
/*             K >= 0. */

/*     REZ     (input) DOUBLE PRECISION array, dimension (K) */
/*     IMZ     (input) DOUBLE PRECISION array, dimension (K) */
/*             The real and imaginary parts of the i-th zero of P(x) */
/*             must be stored in REZ(i) and IMZ(i), respectively, where */
/*             i = 1, 2, ..., K. The zeros may be supplied in any order. */

/*     REP     (output) DOUBLE PRECISION array, dimension (K+1) */
/*     IMP     (output) DOUBLE PRECISION array, dimension (K+1) */
/*             These arrays contain the real and imaginary parts, */
/*             respectively, of the coefficients of P(x) in increasing */
/*             powers of x. If K = 0, then REP(1) is set to one and */
/*             IMP(1) is set to zero. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (2*K+2) */
/*             If K = 0, this array is not referenced. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes the coefficients of the complex K-th degree */
/*     polynomial P(x) as */

/*        P(x) = (x - r(1)) * (x - r(2)) * ... * (x - r(K)) */

/*     where r(i) = (REZ(i),IMZ(i)), using real arithmetic. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01CD by Alan Brown and */
/*     A.J. Geurts. */

/*     REVISIONS */

/*     V. Sima, May 2002. */

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

#line 104 "MC01OD.f"
    /* Parameter adjustments */
#line 104 "MC01OD.f"
    --dwork;
#line 104 "MC01OD.f"
    --imp;
#line 104 "MC01OD.f"
    --rep;
#line 104 "MC01OD.f"
    --imz;
#line 104 "MC01OD.f"
    --rez;
#line 104 "MC01OD.f"

#line 104 "MC01OD.f"
    /* Function Body */
#line 104 "MC01OD.f"
    if (*k < 0) {
#line 105 "MC01OD.f"
	*info = -1;

/*        Error return. */

#line 109 "MC01OD.f"
	i__1 = -(*info);
#line 109 "MC01OD.f"
	xerbla_("MC01OD", &i__1, (ftnlen)6);
#line 110 "MC01OD.f"
	return 0;
#line 111 "MC01OD.f"
    }

/*     Quick return if possible. */

#line 115 "MC01OD.f"
    *info = 0;
#line 116 "MC01OD.f"
    rep[1] = 1.;
#line 117 "MC01OD.f"
    imp[1] = 0.;
#line 118 "MC01OD.f"
    if (*k == 0) {
#line 118 "MC01OD.f"
	return 0;
#line 118 "MC01OD.f"
    }

#line 121 "MC01OD.f"
    k2 = *k + 2;

#line 123 "MC01OD.f"
    i__1 = *k;
#line 123 "MC01OD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 124 "MC01OD.f"
	u = rez[i__];
#line 125 "MC01OD.f"
	v = imz[i__];
#line 126 "MC01OD.f"
	dwork[1] = 0.;
#line 127 "MC01OD.f"
	dwork[k2] = 0.;
#line 128 "MC01OD.f"
	dcopy_(&i__, &rep[1], &c__1, &dwork[2], &c__1);
#line 129 "MC01OD.f"
	dcopy_(&i__, &imp[1], &c__1, &dwork[k2 + 1], &c__1);

#line 131 "MC01OD.f"
	if (u != 0.) {
#line 132 "MC01OD.f"
	    d__1 = -u;
#line 132 "MC01OD.f"
	    daxpy_(&i__, &d__1, &rep[1], &c__1, &dwork[1], &c__1);
#line 133 "MC01OD.f"
	    d__1 = -u;
#line 133 "MC01OD.f"
	    daxpy_(&i__, &d__1, &imp[1], &c__1, &dwork[k2], &c__1);
#line 134 "MC01OD.f"
	}

#line 136 "MC01OD.f"
	if (v != 0.) {
#line 137 "MC01OD.f"
	    daxpy_(&i__, &v, &imp[1], &c__1, &dwork[1], &c__1);
#line 138 "MC01OD.f"
	    d__1 = -v;
#line 138 "MC01OD.f"
	    daxpy_(&i__, &d__1, &rep[1], &c__1, &dwork[k2], &c__1);
#line 139 "MC01OD.f"
	}

#line 141 "MC01OD.f"
	i__2 = i__ + 1;
#line 141 "MC01OD.f"
	dcopy_(&i__2, &dwork[1], &c__1, &rep[1], &c__1);
#line 142 "MC01OD.f"
	i__2 = i__ + 1;
#line 142 "MC01OD.f"
	dcopy_(&i__2, &dwork[k2], &c__1, &imp[1], &c__1);
#line 143 "MC01OD.f"
/* L20: */
#line 143 "MC01OD.f"
    }

#line 145 "MC01OD.f"
    return 0;
/* *** Last line of MC01OD *** */
} /* mc01od_ */

