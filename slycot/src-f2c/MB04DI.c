#line 1 "MB04DI.f"
/* MB04DI.f -- translated by f2c (version 20100827).
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

#line 1 "MB04DI.f"
/* Table of constant values */

static doublereal c_b14 = -1.;

/* Subroutine */ int mb04di_(char *job, char *sgn, integer *n, integer *ilo, 
	doublereal *scale, integer *m, doublereal *v1, integer *ldv1, 
	doublereal *v2, integer *ldv2, integer *info, ftnlen job_len, ftnlen 
	sgn_len)
{
    /* System generated locals */
    integer v1_dim1, v1_offset, v2_dim1, v2_offset, i__1;

    /* Local variables */
    static integer i__, k;
    static logical lsgn, sysw;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *), dswap_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    static logical lperm;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To apply the inverse of a balancing transformation, computed by */
/*     the SLICOT Library routines MB04DD or MB04DS, to a 2*N-by-M matrix */

/*               [   V1   ] */
/*               [        ], */
/*               [ sgn*V2 ] */

/*     where sgn is either +1 or -1. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the type of inverse transformation required: */
/*             = 'N':  do nothing, return immediately; */
/*             = 'P':  do inverse transformation for permutation only; */
/*             = 'S':  do inverse transformation for scaling only; */
/*             = 'B':  do inverse transformations for both permutation */
/*                     and scaling. */
/*             JOB must be the same as the argument JOB supplied to */
/*             MB04DD or MB04DS. */

/*     SGN     CHARACTER*1 */
/*             Specifies the sign to use for V2: */
/*             = 'P':  sgn = +1; */
/*             = 'N':  sgn = -1. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of rows of the matrices V1 and V2. N >= 0. */

/*     ILO     (input) INTEGER */
/*             The integer ILO determined by MB04DD or MB04DS. */
/*             1 <= ILO <= N+1. */

/*     SCALE   (input) DOUBLE PRECISION array, dimension (N) */
/*             Details of the permutation and scaling factors, as */
/*             returned by MB04DD or MB04DS. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrices V1 and V2.  M >= 0. */

/*     V1      (input/output) DOUBLE PRECISION array, dimension (LDV1,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix V1. */
/*             On exit, the leading N-by-M part of this array is */
/*             overwritten by the updated matrix V1 of the transformed */
/*             matrix. */

/*     LDV1    INTEGER */
/*             The leading dimension of the array V1. LDV1 >= max(1,N). */

/*     V2      (input/output) DOUBLE PRECISION array, dimension (LDV2,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the matrix V2. */
/*             On exit, the leading N-by-M part of this array is */
/*             overwritten by the updated matrix V2 of the transformed */
/*             matrix. */

/*     LDV2    INTEGER */
/*             The leading dimension of the array V2. LDV2 >= max(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     REFERENCES */

/*     [1] Benner, P. */
/*         Symplectic balancing of Hamiltonian matrices. */
/*         SIAM J. Sci. Comput., 22 (5), pp. 1885-1904, 2000. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DHABAK). */

/*     KEYWORDS */

/*     Balancing, Hamiltonian matrix, skew-Hamiltonian matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 139 "MB04DI.f"
    /* Parameter adjustments */
#line 139 "MB04DI.f"
    --scale;
#line 139 "MB04DI.f"
    v1_dim1 = *ldv1;
#line 139 "MB04DI.f"
    v1_offset = 1 + v1_dim1;
#line 139 "MB04DI.f"
    v1 -= v1_offset;
#line 139 "MB04DI.f"
    v2_dim1 = *ldv2;
#line 139 "MB04DI.f"
    v2_offset = 1 + v2_dim1;
#line 139 "MB04DI.f"
    v2 -= v2_offset;
#line 139 "MB04DI.f"

#line 139 "MB04DI.f"
    /* Function Body */
#line 139 "MB04DI.f"
    *info = 0;
#line 140 "MB04DI.f"
    lperm = lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
#line 141 "MB04DI.f"
    lscal = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
#line 142 "MB04DI.f"
    lsgn = lsame_(sgn, "N", (ftnlen)1, (ftnlen)1);
#line 143 "MB04DI.f"
    if (! lperm && ! lscal && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 145 "MB04DI.f"
	*info = -1;
#line 146 "MB04DI.f"
    } else if (! lsgn && ! lsame_(sgn, "P", (ftnlen)1, (ftnlen)1)) {
#line 147 "MB04DI.f"
	*info = -2;
#line 148 "MB04DI.f"
    } else if (*n < 0) {
#line 149 "MB04DI.f"
	*info = -3;
#line 150 "MB04DI.f"
    } else if (*ilo < 1 || *ilo > *n + 1) {
#line 151 "MB04DI.f"
	*info = -4;
#line 152 "MB04DI.f"
    } else if (*m < 0) {
#line 153 "MB04DI.f"
	*info = -6;
#line 154 "MB04DI.f"
    } else if (*ldv1 < max(1,*n)) {
#line 155 "MB04DI.f"
	*info = -8;
#line 156 "MB04DI.f"
    } else if (*ldv2 < max(1,*n)) {
#line 157 "MB04DI.f"
	*info = -10;
#line 158 "MB04DI.f"
    }

/*     Return if there were illegal values. */

#line 162 "MB04DI.f"
    if (*info != 0) {
#line 163 "MB04DI.f"
	i__1 = -(*info);
#line 163 "MB04DI.f"
	xerbla_("MB04DI", &i__1, (ftnlen)6);
#line 164 "MB04DI.f"
	return 0;
#line 165 "MB04DI.f"
    }

/*     Quick return if possible. */

#line 169 "MB04DI.f"
    if (*n == 0 || *m == 0 || lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 169 "MB04DI.f"
	return 0;
#line 169 "MB04DI.f"
    }

/*     Inverse scaling. */

#line 174 "MB04DI.f"
    if (lscal) {
#line 175 "MB04DI.f"
	i__1 = *n;
#line 175 "MB04DI.f"
	for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 176 "MB04DI.f"
	    drscl_(m, &scale[i__], &v1[i__ + v1_dim1], ldv1);
#line 177 "MB04DI.f"
/* L20: */
#line 177 "MB04DI.f"
	}
#line 178 "MB04DI.f"
	i__1 = *n;
#line 178 "MB04DI.f"
	for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 179 "MB04DI.f"
	    drscl_(m, &scale[i__], &v2[i__ + v2_dim1], ldv2);
#line 180 "MB04DI.f"
/* L30: */
#line 180 "MB04DI.f"
	}
#line 181 "MB04DI.f"
    }

/*     Inverse permutation. */

#line 185 "MB04DI.f"
    if (lperm) {
#line 186 "MB04DI.f"
	for (i__ = *ilo - 1; i__ >= 1; --i__) {
#line 187 "MB04DI.f"
	    k = (integer) scale[i__];
#line 188 "MB04DI.f"
	    sysw = k > *n;
#line 189 "MB04DI.f"
	    if (sysw) {
#line 189 "MB04DI.f"
		k -= *n;
#line 189 "MB04DI.f"
	    }

#line 192 "MB04DI.f"
	    if (k != i__) {

/*              Exchange rows k <-> i. */

#line 196 "MB04DI.f"
		dswap_(m, &v1[i__ + v1_dim1], ldv1, &v1[k + v1_dim1], ldv1);
#line 197 "MB04DI.f"
		dswap_(m, &v2[i__ + v2_dim1], ldv2, &v2[k + v2_dim1], ldv2);
#line 198 "MB04DI.f"
	    }

#line 200 "MB04DI.f"
	    if (sysw) {

/*              Exchange V1(k,:) <-> V2(k,:). */

#line 204 "MB04DI.f"
		dswap_(m, &v1[k + v1_dim1], ldv1, &v2[k + v2_dim1], ldv2);
#line 205 "MB04DI.f"
		if (lsgn) {
#line 206 "MB04DI.f"
		    dscal_(m, &c_b14, &v2[k + v2_dim1], ldv2);
#line 207 "MB04DI.f"
		} else {
#line 208 "MB04DI.f"
		    dscal_(m, &c_b14, &v1[k + v1_dim1], ldv1);
#line 209 "MB04DI.f"
		}
#line 210 "MB04DI.f"
	    }
#line 211 "MB04DI.f"
/* L40: */
#line 211 "MB04DI.f"
	}
#line 212 "MB04DI.f"
    }

#line 214 "MB04DI.f"
    return 0;
/* *** Last line of MB04DI *** */
} /* mb04di_ */

