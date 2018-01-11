#line 1 "TC01OD.f"
/* TC01OD.f -- translated by f2c (version 20100827).
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

#line 1 "TC01OD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int tc01od_(char *leri, integer *m, integer *p, integer *
	indlim, doublereal *pcoeff, integer *ldpco1, integer *ldpco2, 
	doublereal *qcoeff, integer *ldqco1, integer *ldqco2, integer *info, 
	ftnlen leri_len)
{
    /* System generated locals */
    integer pcoeff_dim1, pcoeff_dim2, pcoeff_offset, qcoeff_dim1, qcoeff_dim2,
	     qcoeff_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer j, k, porm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lleri;
    static integer mplim;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer minmp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), xerbla_(char *, integer *, ftnlen);


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

/*     To find the dual right (left) polynomial matrix representation of */
/*     a given left (right) polynomial matrix representation, where the */
/*     right and left polynomial matrix representations are of the form */
/*     Q(s)*inv(P(s)) and inv(P(s))*Q(s) respectively. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LERI    CHARACTER*1 */
/*             Indicates whether a left or right matrix fraction is input */
/*             as follows: */
/*             = 'L':  A left matrix fraction is input; */
/*             = 'R':  A right matrix fraction is input. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     INDLIM  (input) INTEGER */
/*             The highest value of K for which PCOEFF(.,.,K) and */
/*             QCOEFF(.,.,K) are to be transposed. */
/*             K = kpcoef + 1, where kpcoef is the maximum degree of the */
/*             polynomials in P(s).  INDLIM >= 1. */

/*     PCOEFF  (input/output) DOUBLE PRECISION array, dimension */
/*             (LDPCO1,LDPCO2,INDLIM) */
/*             If LERI = 'L' then porm = P, otherwise porm = M. */
/*             On entry, the leading porm-by-porm-by-INDLIM part of this */
/*             array must contain the coefficients of the denominator */
/*             matrix P(s). */
/*             PCOEFF(I,J,K) is the coefficient in s**(INDLIM-K) of */
/*             polynomial (I,J) of P(s), where K = 1,2,...,INDLIM. */
/*             On exit, the leading porm-by-porm-by-INDLIM part of this */
/*             array contains the coefficients of the denominator matrix */
/*             P'(s) of the dual system. */

/*     LDPCO1  INTEGER */
/*             The leading dimension of array PCOEFF. */
/*             LDPCO1 >= MAX(1,P) if LERI = 'L', */
/*             LDPCO1 >= MAX(1,M) if LERI = 'R'. */

/*     LDPCO2  INTEGER */
/*             The second dimension of array PCOEFF. */
/*             LDPCO2 >= MAX(1,P) if LERI = 'L', */
/*             LDPCO2 >= MAX(1,M) if LERI = 'R'. */

/*     QCOEFF  (input/output) DOUBLE PRECISION array, dimension */
/*             (LDQCO1,LDQCO2,INDLIM) */
/*             On entry, the leading P-by-M-by-INDLIM part of this array */
/*             must contain the coefficients of the numerator matrix */
/*             Q(s). */
/*             QCOEFF(I,J,K) is the coefficient in s**(INDLIM-K) of */
/*             polynomial (I,J) of Q(s), where K = 1,2,...,INDLIM. */
/*             On exit, the leading M-by-P-by-INDLIM part of the array */
/*             contains the coefficients of the numerator matrix Q'(s) */
/*             of the dual system. */

/*     LDQCO1  INTEGER */
/*             The leading dimension of array QCOEFF. */
/*             LDQCO1 >= MAX(1,M,P). */

/*     LDQCO2  INTEGER */
/*             The second dimension of array QCOEFF. */
/*             LDQCO2 >= MAX(1,M,P). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If the given M-input/P-output left (right) polynomial matrix */
/*     representation has numerator matrix Q(s) and denominator matrix */
/*     P(s), its dual P-input/M-output right (left) polynomial matrix */
/*     representation simply has numerator matrix Q'(s) and denominator */
/*     matrix P'(s). */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TC01CD by T.W.C.Williams, Kingston */
/*     Polytechnic, United Kingdom, March 1982. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Coprime matrix fraction, elementary polynomial operations, */
/*     polynomial matrix, state-space representation, transfer matrix. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 152 "TC01OD.f"
    /* Parameter adjustments */
#line 152 "TC01OD.f"
    pcoeff_dim1 = *ldpco1;
#line 152 "TC01OD.f"
    pcoeff_dim2 = *ldpco2;
#line 152 "TC01OD.f"
    pcoeff_offset = 1 + pcoeff_dim1 * (1 + pcoeff_dim2);
#line 152 "TC01OD.f"
    pcoeff -= pcoeff_offset;
#line 152 "TC01OD.f"
    qcoeff_dim1 = *ldqco1;
#line 152 "TC01OD.f"
    qcoeff_dim2 = *ldqco2;
#line 152 "TC01OD.f"
    qcoeff_offset = 1 + qcoeff_dim1 * (1 + qcoeff_dim2);
#line 152 "TC01OD.f"
    qcoeff -= qcoeff_offset;
#line 152 "TC01OD.f"

#line 152 "TC01OD.f"
    /* Function Body */
#line 152 "TC01OD.f"
    *info = 0;
#line 153 "TC01OD.f"
    lleri = lsame_(leri, "L", (ftnlen)1, (ftnlen)1);
#line 154 "TC01OD.f"
    mplim = max(*m,*p);
#line 155 "TC01OD.f"
    minmp = min(*m,*p);

/*     Test the input scalar arguments. */

#line 159 "TC01OD.f"
    if (! lleri && ! lsame_(leri, "R", (ftnlen)1, (ftnlen)1)) {
#line 160 "TC01OD.f"
	*info = -1;
#line 161 "TC01OD.f"
    } else if (*m < 0) {
#line 162 "TC01OD.f"
	*info = -2;
#line 163 "TC01OD.f"
    } else if (*p < 0) {
#line 164 "TC01OD.f"
	*info = -3;
#line 165 "TC01OD.f"
    } else if (*indlim < 1) {
#line 166 "TC01OD.f"
	*info = -4;
#line 167 "TC01OD.f"
    } else if (lleri && *ldpco1 < max(1,*p) || ! lleri && *ldpco1 < max(1,*m))
	     {
#line 169 "TC01OD.f"
	*info = -6;
#line 170 "TC01OD.f"
    } else if (lleri && *ldpco2 < max(1,*p) || ! lleri && *ldpco2 < max(1,*m))
	     {
#line 172 "TC01OD.f"
	*info = -7;
#line 173 "TC01OD.f"
    } else if (*ldqco1 < max(1,mplim)) {
#line 174 "TC01OD.f"
	*info = -9;
#line 175 "TC01OD.f"
    } else if (*ldqco2 < max(1,mplim)) {
#line 176 "TC01OD.f"
	*info = -10;
#line 177 "TC01OD.f"
    }

#line 179 "TC01OD.f"
    if (*info != 0) {

/*        Error return. */

#line 183 "TC01OD.f"
	i__1 = -(*info);
#line 183 "TC01OD.f"
	xerbla_("TC01OD", &i__1, (ftnlen)6);
#line 184 "TC01OD.f"
	return 0;
#line 185 "TC01OD.f"
    }

/*     Quick return if possible. */

#line 189 "TC01OD.f"
    if (*m == 0 || *p == 0) {
#line 189 "TC01OD.f"
	return 0;
#line 189 "TC01OD.f"
    }

#line 192 "TC01OD.f"
    if (mplim != 1) {

/*        Non-scalar system: transpose numerator matrix Q(s). */

#line 196 "TC01OD.f"
	i__1 = *indlim;
#line 196 "TC01OD.f"
	for (k = 1; k <= i__1; ++k) {

#line 198 "TC01OD.f"
	    i__2 = mplim;
#line 198 "TC01OD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 199 "TC01OD.f"
		if (j < minmp) {
#line 200 "TC01OD.f"
		    i__3 = minmp - j;
#line 200 "TC01OD.f"
		    dswap_(&i__3, &qcoeff[j + 1 + (j + k * qcoeff_dim2) * 
			    qcoeff_dim1], &c__1, &qcoeff[j + (j + 1 + k * 
			    qcoeff_dim2) * qcoeff_dim1], ldqco1);
#line 202 "TC01OD.f"
		} else if (j > *p) {
#line 203 "TC01OD.f"
		    dcopy_(p, &qcoeff[(j + k * qcoeff_dim2) * qcoeff_dim1 + 1]
			    , &c__1, &qcoeff[j + (k * qcoeff_dim2 + 1) * 
			    qcoeff_dim1], ldqco1);
#line 205 "TC01OD.f"
		} else if (j > *m) {
#line 206 "TC01OD.f"
		    dcopy_(m, &qcoeff[j + (k * qcoeff_dim2 + 1) * qcoeff_dim1]
			    , ldqco1, &qcoeff[(j + k * qcoeff_dim2) * 
			    qcoeff_dim1 + 1], &c__1);
#line 208 "TC01OD.f"
		}
#line 209 "TC01OD.f"
/* L10: */
#line 209 "TC01OD.f"
	    }

#line 211 "TC01OD.f"
/* L20: */
#line 211 "TC01OD.f"
	}

/*        Find dimension of denominator matrix P(s): M (P) for */
/*        right (left) polynomial matrix representation. */

#line 216 "TC01OD.f"
	porm = *m;
#line 217 "TC01OD.f"
	if (lleri) {
#line 217 "TC01OD.f"
	    porm = *p;
#line 217 "TC01OD.f"
	}
#line 218 "TC01OD.f"
	if (porm != 1) {

/*           Non-scalar P(s): transpose it. */

#line 222 "TC01OD.f"
	    i__1 = *indlim;
#line 222 "TC01OD.f"
	    for (k = 1; k <= i__1; ++k) {

#line 224 "TC01OD.f"
		i__2 = porm - 1;
#line 224 "TC01OD.f"
		for (j = 1; j <= i__2; ++j) {
#line 225 "TC01OD.f"
		    i__3 = porm - j;
#line 225 "TC01OD.f"
		    dswap_(&i__3, &pcoeff[j + 1 + (j + k * pcoeff_dim2) * 
			    pcoeff_dim1], &c__1, &pcoeff[j + (j + 1 + k * 
			    pcoeff_dim2) * pcoeff_dim1], ldpco1);
#line 227 "TC01OD.f"
/* L30: */
#line 227 "TC01OD.f"
		}

#line 229 "TC01OD.f"
/* L40: */
#line 229 "TC01OD.f"
	    }

#line 231 "TC01OD.f"
	}
#line 232 "TC01OD.f"
    }

#line 234 "TC01OD.f"
    return 0;
/* *** Last line of TC01OD *** */
} /* tc01od_ */

