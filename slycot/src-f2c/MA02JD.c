#line 1 "MA02JD.f"
/* MA02JD.f -- translated by f2c (version 20100827).
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

#line 1 "MA02JD.f"
/* Table of constant values */

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static doublereal c_b27 = -1.;

doublereal ma02jd_(logical *ltran1, logical *ltran2, integer *n, doublereal *
	q1, integer *ldq1, doublereal *q2, integer *ldq2, doublereal *res, 
	integer *ldres)
{
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, res_dim1, res_offset, 
	    i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal temp;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal dummy[1];
    extern doublereal dlapy2_(doublereal *, doublereal *), dlange_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    ftnlen);


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

/*     To compute || Q^T Q - I ||_F for a matrix of the form */

/*                       [  op( Q1 )  op( Q2 ) ] */
/*                  Q =  [                     ], */
/*                       [ -op( Q2 )  op( Q1 ) ] */

/*     where Q1 and Q2 are N-by-N matrices. This residual can be used to */
/*     test wether Q is numerically an orthogonal symplectic matrix. */

/*     FUNCTION VALUE */

/*     MA02JD  DOUBLE PRECISION */
/*             The computed residual. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LTRAN1  LOGICAL */
/*             Specifies the form of op( Q1 ) as follows: */
/*             = .FALSE.:  op( Q1 ) = Q1; */
/*             = .TRUE. :  op( Q1 ) = Q1'. */

/*     LTRAN2  LOGICAL */
/*             Specifies the form of op( Q2 ) as follows: */
/*             = .FALSE.:  op( Q2 ) = Q2; */
/*             = .TRUE. :  op( Q2 ) = Q2'. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices Q1 and Q2.  N >= 0. */

/*     Q1      (input) DOUBLE PRECISION array, dimension (LDQ1,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix op( Q1 ). */

/*     LDQ1    INTEGER */
/*             The leading dimension of the array Q1.  LDQ1 >= MAX(1,N). */

/*     Q2      (input) DOUBLE PRECISION array, dimension (LDQ2,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix op( Q2 ). */

/*     LDQ2    INTEGER */
/*             The leading dimension of the array Q2.  LDQ2 >= MAX(1,N). */

/*     Workspace */

/*     RES     DOUBLE PRECISION array, dimension (LDRES,N) */

/*     LDRES   INTEGER */
/*             The leading dimension of the array RES.  LDRES >= MAX(1,N). */

/*     METHOD */

/*     The routine computes the residual by simple elementary operations. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAORS). */

/*     KEYWORDS */

/*     Elementary operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Subroutines .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

#line 120 "MA02JD.f"
    /* Parameter adjustments */
#line 120 "MA02JD.f"
    q1_dim1 = *ldq1;
#line 120 "MA02JD.f"
    q1_offset = 1 + q1_dim1;
#line 120 "MA02JD.f"
    q1 -= q1_offset;
#line 120 "MA02JD.f"
    q2_dim1 = *ldq2;
#line 120 "MA02JD.f"
    q2_offset = 1 + q2_dim1;
#line 120 "MA02JD.f"
    q2 -= q2_offset;
#line 120 "MA02JD.f"
    res_dim1 = *ldres;
#line 120 "MA02JD.f"
    res_offset = 1 + res_dim1;
#line 120 "MA02JD.f"
    res -= res_offset;
#line 120 "MA02JD.f"

#line 120 "MA02JD.f"
    /* Function Body */
#line 120 "MA02JD.f"
    if (*ltran1) {
#line 121 "MA02JD.f"
	dgemm_("No Transpose", "Transpose", n, n, n, &c_b4, &q1[q1_offset], 
		ldq1, &q1[q1_offset], ldq1, &c_b5, &res[res_offset], ldres, (
		ftnlen)12, (ftnlen)9);
#line 123 "MA02JD.f"
    } else {
#line 124 "MA02JD.f"
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b4, &q1[q1_offset], 
		ldq1, &q1[q1_offset], ldq1, &c_b5, &res[res_offset], ldres, (
		ftnlen)9, (ftnlen)12);
#line 126 "MA02JD.f"
    }
#line 127 "MA02JD.f"
    if (*ltran2) {
#line 128 "MA02JD.f"
	dgemm_("No Transpose", "Transpose", n, n, n, &c_b4, &q2[q2_offset], 
		ldq2, &q2[q2_offset], ldq2, &c_b4, &res[res_offset], ldres, (
		ftnlen)12, (ftnlen)9);
#line 130 "MA02JD.f"
    } else {
#line 131 "MA02JD.f"
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b4, &q2[q2_offset], 
		ldq2, &q2[q2_offset], ldq2, &c_b4, &res[res_offset], ldres, (
		ftnlen)9, (ftnlen)12);
#line 133 "MA02JD.f"
    }
#line 134 "MA02JD.f"
    i__1 = *n;
#line 134 "MA02JD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 135 "MA02JD.f"
	res[i__ + i__ * res_dim1] += -1.;
#line 136 "MA02JD.f"
/* L10: */
#line 136 "MA02JD.f"
    }
#line 137 "MA02JD.f"
    temp = dlange_("Frobenius", n, n, &res[res_offset], ldres, dummy, (ftnlen)
	    9);
#line 138 "MA02JD.f"
    if (*ltran1 && *ltran2) {
#line 139 "MA02JD.f"
	dgemm_("No Transpose", "Transpose", n, n, n, &c_b4, &q2[q2_offset], 
		ldq2, &q1[q1_offset], ldq1, &c_b5, &res[res_offset], ldres, (
		ftnlen)12, (ftnlen)9);
#line 141 "MA02JD.f"
	dgemm_("No Transpose", "Transpose", n, n, n, &c_b4, &q1[q1_offset], 
		ldq1, &q2[q2_offset], ldq2, &c_b27, &res[res_offset], ldres, (
		ftnlen)12, (ftnlen)9);
#line 143 "MA02JD.f"
    } else if (*ltran1) {
#line 144 "MA02JD.f"
	dgemm_("Transpose", "Transpose", n, n, n, &c_b4, &q2[q2_offset], ldq2,
		 &q1[q1_offset], ldq1, &c_b5, &res[res_offset], ldres, (
		ftnlen)9, (ftnlen)9);
#line 146 "MA02JD.f"
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b4, &q1[q1_offset],
		 ldq1, &q2[q2_offset], ldq2, &c_b27, &res[res_offset], ldres, 
		(ftnlen)12, (ftnlen)12);
#line 148 "MA02JD.f"
    } else if (*ltran2) {
#line 149 "MA02JD.f"
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b4, &q2[q2_offset],
		 ldq2, &q1[q1_offset], ldq1, &c_b5, &res[res_offset], ldres, (
		ftnlen)12, (ftnlen)12);
#line 151 "MA02JD.f"
	dgemm_("Transpose", "Transpose", n, n, n, &c_b4, &q1[q1_offset], ldq1,
		 &q2[q2_offset], ldq2, &c_b27, &res[res_offset], ldres, (
		ftnlen)9, (ftnlen)9);
#line 153 "MA02JD.f"
    } else {
#line 154 "MA02JD.f"
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b4, &q2[q2_offset], 
		ldq2, &q1[q1_offset], ldq1, &c_b5, &res[res_offset], ldres, (
		ftnlen)9, (ftnlen)12);
#line 156 "MA02JD.f"
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b4, &q1[q1_offset], 
		ldq1, &q2[q2_offset], ldq2, &c_b27, &res[res_offset], ldres, (
		ftnlen)9, (ftnlen)12);
#line 158 "MA02JD.f"
    }
#line 159 "MA02JD.f"
    d__1 = dlange_("Frobenius", n, n, &res[res_offset], ldres, dummy, (ftnlen)
	    9);
#line 159 "MA02JD.f"
    temp = dlapy2_(&temp, &d__1);
#line 161 "MA02JD.f"
    ret_val = sqrt(2.) * temp;
#line 162 "MA02JD.f"
    return ret_val;
/* *** Last line of MA02JD *** */
} /* ma02jd_ */
