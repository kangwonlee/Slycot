#line 1 "MD03BF.f"
/* MD03BF.f -- translated by f2c (version 20100827).
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

#line 1 "MD03BF.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int md03bf_(integer *iflag, integer *m, integer *n, integer *
	ipar, integer *lipar, doublereal *dpar1, integer *ldpar1, doublereal *
	dpar2, integer *ldpar2, doublereal *x, integer *nfevl, doublereal *e, 
	doublereal *j, integer *ldj, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* Initialized data */

    static doublereal y[15] = { .14,.18,.22,.25,.29,.32,.35,.39,.37,.58,.73,
	    .96,1.34,2.1,4.39 };

    /* System generated locals */
    integer j_dim1, j_offset;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__;
    static doublereal err, tmp1, tmp2, tmp3, tmp4;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 1, 0, "(' Norm of current error = ', D15.6)", 
	    0 };



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

/*     This is the FCN routine for solving a standard nonlinear least */
/*     squares problem using SLICOT Library routine MD03BD. See the */
/*     parameter FCN in the routine MD03BD for the description of */
/*     parameters. */

/*     The example programmed in this routine is adapted from that */
/*     accompanying the MINPACK routine LMDER. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. NOUT is the unit number for printing intermediate results .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. DATA Statements .. */
#line 54 "MD03BF.f"
    /* Parameter adjustments */
#line 54 "MD03BF.f"
    --ipar;
#line 54 "MD03BF.f"
    --dpar1;
#line 54 "MD03BF.f"
    --dpar2;
#line 54 "MD03BF.f"
    --x;
#line 54 "MD03BF.f"
    --e;
#line 54 "MD03BF.f"
    j_dim1 = *ldj;
#line 54 "MD03BF.f"
    j_offset = 1 + j_dim1;
#line 54 "MD03BF.f"
    j -= j_offset;
#line 54 "MD03BF.f"
    --dwork;
#line 54 "MD03BF.f"

#line 54 "MD03BF.f"
    /* Function Body */

/*     .. Executable Statements .. */

#line 62 "MD03BF.f"
    *info = 0;
#line 63 "MD03BF.f"
    if (*iflag == 1) {

/*        Compute the error function values. */

#line 67 "MD03BF.f"
	for (i__ = 1; i__ <= 15; ++i__) {
#line 68 "MD03BF.f"
	    tmp1 = (doublereal) i__;
#line 69 "MD03BF.f"
	    tmp2 = (doublereal) (16 - i__);
#line 70 "MD03BF.f"
	    if (i__ > 8) {
#line 71 "MD03BF.f"
		tmp3 = tmp2;
#line 72 "MD03BF.f"
	    } else {
#line 73 "MD03BF.f"
		tmp3 = tmp1;
#line 74 "MD03BF.f"
	    }
#line 75 "MD03BF.f"
	    e[i__] = y[i__ - 1] - (x[1] + tmp1 / (x[2] * tmp2 + x[3] * tmp3));
#line 76 "MD03BF.f"
/* L10: */
#line 76 "MD03BF.f"
	}

#line 78 "MD03BF.f"
    } else if (*iflag == 2) {

/*        Compute the Jacobian. */

#line 82 "MD03BF.f"
	for (i__ = 1; i__ <= 15; ++i__) {
#line 83 "MD03BF.f"
	    tmp1 = (doublereal) i__;
#line 84 "MD03BF.f"
	    tmp2 = (doublereal) (16 - i__);
#line 85 "MD03BF.f"
	    if (i__ > 8) {
#line 86 "MD03BF.f"
		tmp3 = tmp2;
#line 87 "MD03BF.f"
	    } else {
#line 88 "MD03BF.f"
		tmp3 = tmp1;
#line 89 "MD03BF.f"
	    }
/* Computing 2nd power */
#line 90 "MD03BF.f"
	    d__1 = x[2] * tmp2 + x[3] * tmp3;
#line 90 "MD03BF.f"
	    tmp4 = d__1 * d__1;
#line 91 "MD03BF.f"
	    j[i__ + j_dim1] = -1.;
#line 92 "MD03BF.f"
	    j[i__ + (j_dim1 << 1)] = tmp1 * tmp2 / tmp4;
#line 93 "MD03BF.f"
	    j[i__ + j_dim1 * 3] = tmp1 * tmp3 / tmp4;
#line 94 "MD03BF.f"
/* L30: */
#line 94 "MD03BF.f"
	}

#line 96 "MD03BF.f"
	*nfevl = 0;

#line 98 "MD03BF.f"
    } else if (*iflag == 3) {

/*        Set the parameter LDJ, the length of the array J, and the sizes */
/*        of the workspace for FCN (IFLAG = 1 or 2), MD03BA and MD03BB. */

#line 103 "MD03BF.f"
	*ldj = *m;
#line 104 "MD03BF.f"
	ipar[1] = *m * *n;
#line 105 "MD03BF.f"
	ipar[2] = 0;
#line 106 "MD03BF.f"
	ipar[3] = 0;
#line 107 "MD03BF.f"
	ipar[4] = (*n << 2) + 1;
#line 108 "MD03BF.f"
	ipar[5] = *n << 2;

#line 110 "MD03BF.f"
    } else if (*iflag == 0) {

/*        Special call for printing intermediate results. */

#line 114 "MD03BF.f"
	err = dnrm2_(m, &e[1], &c__1);
#line 115 "MD03BF.f"
	s_wsfe(&io___8);
#line 115 "MD03BF.f"
	do_fio(&c__1, (char *)&err, (ftnlen)sizeof(doublereal));
#line 115 "MD03BF.f"
	e_wsfe();

#line 117 "MD03BF.f"
    }

#line 119 "MD03BF.f"
    return 0;

/* *** Last line of MD03BF *** */
} /* md03bf_ */

