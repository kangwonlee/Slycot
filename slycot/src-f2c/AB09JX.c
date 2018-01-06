#line 1 "AB09JX.f"
/* AB09JX.f -- translated by f2c (version 20100827).
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

#line 1 "AB09JX.f"
/* Subroutine */ int ab09jx_(char *dico, char *stdom, char *evtype, integer *
	n, doublereal *alpha, doublereal *er, doublereal *ei, doublereal *ed, 
	doublereal *tolinf, integer *info, ftnlen dico_len, ftnlen stdom_len, 
	ftnlen evtype_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__;
    static logical stab;
    static doublereal rpev, scale, absev;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical recevp, stdevp;


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

/*     To check stability/antistability of finite eigenvalues with */
/*     respect to a given stability domain. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the stability domain as follows: */
/*             = 'C':  for a continuous-time system; */
/*             = 'D':  for a discrete-time system. */

/*     STDOM   CHARACTER*1 */
/*             Specifies whether the domain of interest is of stability */
/*             type (left part of complex plane or inside of a circle) */
/*             or of instability type (right part of complex plane or */
/*             outside of a circle) as follows: */
/*             = 'S':  stability type domain; */
/*             = 'U':  instability type domain. */

/*     EVTYPE  CHARACTER*1 */
/*             Specifies whether the eigenvalues arise from a standard */
/*             or a generalized eigenvalue problem as follows: */
/*             = 'S':  standard eigenvalue problem; */
/*             = 'G':  generalized eigenvalue problem; */
/*             = 'R':  reciprocal generalized eigenvalue problem. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of vectors ER, EI and ED.  N >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             Specifies the boundary of the domain of interest for the */
/*             eigenvalues. For a continuous-time system */
/*             (DICO = 'C'), ALPHA is the boundary value for the real */
/*             parts of eigenvalues, while for a discrete-time system */
/*             (DICO = 'D'), ALPHA >= 0 represents the boundary value for */
/*             the moduli of eigenvalues. */

/*     ER, EI, (input) DOUBLE PRECISION arrays, dimension (N) */
/*     ED      If EVTYPE = 'S', ER(j) + EI(j)*i, j = 1,...,N, are */
/*             the eigenvalues of a real matrix. */
/*             ED is not referenced and is implicitly considered as */
/*             a vector having all elements equal to one. */
/*             If EVTYPE = 'G' or EVTYPE = 'R', (ER(j) + EI(j)*i)/ED(j), */
/*             j = 1,...,N, are the generalized eigenvalues of a pair of */
/*             real matrices. If ED(j) is zero, then the j-th generalized */
/*             eigenvalue is infinite. */
/*             Complex conjugate pairs of eigenvalues must appear */
/*             consecutively. */

/*     Tolerances */

/*     TOLINF  DOUBLE PRECISION */
/*             If EVTYPE = 'G' or 'R', TOLINF contains the tolerance for */
/*             detecting infinite generalized eigenvalues. */
/*             0 <= TOLINF < 1. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             =  0:  successful exit, i.e., all eigenvalues lie within */
/*                    the domain of interest defined by DICO, STDOM */
/*                    and ALPHA; */
/*             <  0:  if INFO = -i, the i-th argument had an illegal */
/*                    value; */
/*             =  1:  some eigenvalues lie outside the domain of interest */
/*                    defined by DICO, STDOM and ALPHA. */
/*     METHOD */

/*     The domain of interest for an eigenvalue lambda is defined by the */
/*     parameters ALPHA, DICO and STDOM as follows: */
/*        - for a continuous-time system (DICO = 'C'): */
/*               Real(lambda) < ALPHA if STDOM = 'S'; */
/*               Real(lambda) > ALPHA if STDOM = 'U'; */
/*        - for a discrete-time system (DICO = 'D'): */
/*               Abs(lambda) < ALPHA if STDOM = 'S'; */
/*               Abs(lambda) > ALPHA if STDOM = 'U'. */
/*     If EVTYPE = 'R', the same conditions apply for 1/lambda. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001. */

/*     KEYWORDS */

/*     Stability. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 142 "AB09JX.f"
    /* Parameter adjustments */
#line 142 "AB09JX.f"
    --ed;
#line 142 "AB09JX.f"
    --ei;
#line 142 "AB09JX.f"
    --er;
#line 142 "AB09JX.f"

#line 142 "AB09JX.f"
    /* Function Body */
#line 142 "AB09JX.f"
    *info = 0;
#line 143 "AB09JX.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 144 "AB09JX.f"
    stab = lsame_(stdom, "S", (ftnlen)1, (ftnlen)1);
#line 145 "AB09JX.f"
    stdevp = lsame_(evtype, "S", (ftnlen)1, (ftnlen)1);
#line 146 "AB09JX.f"
    recevp = lsame_(evtype, "R", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input arguments. */

#line 150 "AB09JX.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 151 "AB09JX.f"
	*info = -1;
#line 152 "AB09JX.f"
    } else if (! (stab || lsame_(stdom, "U", (ftnlen)1, (ftnlen)1))) {
#line 153 "AB09JX.f"
	*info = -2;
#line 154 "AB09JX.f"
    } else if (! (stdevp || lsame_(evtype, "G", (ftnlen)1, (ftnlen)1) || 
	    recevp)) {
#line 156 "AB09JX.f"
	*info = -3;
#line 157 "AB09JX.f"
    } else if (*n < 0) {
#line 158 "AB09JX.f"
	*info = -4;
#line 159 "AB09JX.f"
    } else if (discr && *alpha < 0.) {
#line 160 "AB09JX.f"
	*info = -5;
#line 161 "AB09JX.f"
    } else if (*tolinf < 0. || *tolinf >= 1.) {
#line 162 "AB09JX.f"
	*info = -9;
#line 163 "AB09JX.f"
    }

#line 165 "AB09JX.f"
    if (*info != 0) {

/*        Error return. */

#line 169 "AB09JX.f"
	i__1 = -(*info);
#line 169 "AB09JX.f"
	xerbla_("AB09JX", &i__1, (ftnlen)6);
#line 170 "AB09JX.f"
	return 0;
#line 171 "AB09JX.f"
    }

/*     Quick return if possible. */

#line 175 "AB09JX.f"
    if (*n == 0) {
#line 175 "AB09JX.f"
	return 0;
#line 175 "AB09JX.f"
    }

#line 178 "AB09JX.f"
    if (stab) {

/*        Check the stability of finite eigenvalues. */

#line 182 "AB09JX.f"
	scale = 1.;
#line 183 "AB09JX.f"
	if (discr) {
#line 184 "AB09JX.f"
	    i__1 = *n;
#line 184 "AB09JX.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 185 "AB09JX.f"
		absev = dlapy2_(&er[i__], &ei[i__]);
#line 186 "AB09JX.f"
		if (recevp) {
#line 187 "AB09JX.f"
		    scale = absev;
#line 188 "AB09JX.f"
		    absev = (d__1 = ed[i__], abs(d__1));
#line 189 "AB09JX.f"
		} else if (! stdevp) {
#line 190 "AB09JX.f"
		    scale = ed[i__];
#line 191 "AB09JX.f"
		}
#line 192 "AB09JX.f"
		if (abs(scale) > *tolinf && absev >= *alpha * scale) {
#line 194 "AB09JX.f"
		    *info = 1;
#line 195 "AB09JX.f"
		    return 0;
#line 196 "AB09JX.f"
		}
#line 197 "AB09JX.f"
/* L10: */
#line 197 "AB09JX.f"
	    }
#line 198 "AB09JX.f"
	} else {
#line 199 "AB09JX.f"
	    i__1 = *n;
#line 199 "AB09JX.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 200 "AB09JX.f"
		rpev = er[i__];
#line 201 "AB09JX.f"
		if (recevp) {
#line 202 "AB09JX.f"
		    scale = rpev;
#line 203 "AB09JX.f"
		    rpev = ed[i__];
#line 204 "AB09JX.f"
		} else if (! stdevp) {
#line 205 "AB09JX.f"
		    scale = ed[i__];
#line 206 "AB09JX.f"
		}
#line 207 "AB09JX.f"
		if (abs(scale) > *tolinf && rpev >= *alpha * scale) {
#line 209 "AB09JX.f"
		    *info = 1;
#line 210 "AB09JX.f"
		    return 0;
#line 211 "AB09JX.f"
		}
#line 212 "AB09JX.f"
/* L20: */
#line 212 "AB09JX.f"
	    }
#line 213 "AB09JX.f"
	}
#line 214 "AB09JX.f"
    } else {

/*        Check the anti-stability of finite eigenvalues. */

#line 218 "AB09JX.f"
	if (discr) {
#line 219 "AB09JX.f"
	    i__1 = *n;
#line 219 "AB09JX.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 220 "AB09JX.f"
		absev = dlapy2_(&er[i__], &ei[i__]);
#line 221 "AB09JX.f"
		if (recevp) {
#line 222 "AB09JX.f"
		    scale = absev;
#line 223 "AB09JX.f"
		    absev = (d__1 = ed[i__], abs(d__1));
#line 224 "AB09JX.f"
		} else if (! stdevp) {
#line 225 "AB09JX.f"
		    scale = ed[i__];
#line 226 "AB09JX.f"
		}
#line 227 "AB09JX.f"
		if (abs(scale) > *tolinf && absev <= *alpha * scale) {
#line 229 "AB09JX.f"
		    *info = 1;
#line 230 "AB09JX.f"
		    return 0;
#line 231 "AB09JX.f"
		}
#line 232 "AB09JX.f"
/* L30: */
#line 232 "AB09JX.f"
	    }
#line 233 "AB09JX.f"
	} else {
#line 234 "AB09JX.f"
	    i__1 = *n;
#line 234 "AB09JX.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 235 "AB09JX.f"
		rpev = er[i__];
#line 236 "AB09JX.f"
		if (recevp) {
#line 237 "AB09JX.f"
		    scale = rpev;
#line 238 "AB09JX.f"
		    rpev = ed[i__];
#line 239 "AB09JX.f"
		} else if (! stdevp) {
#line 240 "AB09JX.f"
		    scale = ed[i__];
#line 241 "AB09JX.f"
		}
#line 242 "AB09JX.f"
		if (abs(scale) > *tolinf && rpev <= *alpha * scale) {
#line 244 "AB09JX.f"
		    *info = 1;
#line 245 "AB09JX.f"
		    return 0;
#line 246 "AB09JX.f"
		}
#line 247 "AB09JX.f"
/* L40: */
#line 247 "AB09JX.f"
	    }
#line 248 "AB09JX.f"
	}
#line 249 "AB09JX.f"
    }

#line 251 "AB09JX.f"
    return 0;
/* *** Last line of AB09JX *** */
} /* ab09jx_ */

