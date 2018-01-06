#line 1 "SB09MD.f"
/* SB09MD.f -- translated by f2c (version 20100827).
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

#line 1 "SB09MD.f"
/* Subroutine */ int sb09md_(integer *n, integer *nc, integer *nb, doublereal 
	*h1, integer *ldh1, doublereal *h2, integer *ldh2, doublereal *ss, 
	integer *ldss, doublereal *se, integer *ldse, doublereal *pre, 
	integer *ldpre, doublereal *tol, integer *info)
{
    /* System generated locals */
    integer h1_dim1, h1_offset, h2_dim1, h2_offset, pre_dim1, pre_offset, 
	    se_dim1, se_offset, ss_dim1, ss_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k;
    static doublereal var, sse, sss, vare, epso, toler;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noflow;


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

/*     To compare two multivariable sequences M1(k) and M2(k) for */
/*     k = 1,2,...,N, and evaluate their closeness. Each of the */
/*     parameters M1(k) and M2(k) is an NC by NB matrix. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of parameters.  N >= 0. */

/*     NC      (input) INTEGER */
/*             The number of rows in M1(k) and M2(k).  NC >= 0. */

/*     NB      (input) INTEGER */
/*             The number of columns in M1(k) and M2(k).  NB >= 0. */

/*     H1      (input) DOUBLE PRECISION array, dimension (LDH1,N*NB) */
/*             The leading NC-by-N*NB part of this array must contain */
/*             the multivariable sequence M1(k), where k = 1,2,...,N. */
/*             Each parameter M1(k) is an NC-by-NB matrix, whose */
/*             (i,j)-th element must be stored in H1(i,(k-1)*NB+j) for */
/*             i = 1,2,...,NC and j = 1,2,...,NB. */

/*     LDH1    INTEGER */
/*             The leading dimension of array H1.  LDH1 >= MAX(1,NC). */

/*     H2      (input) DOUBLE PRECISION array, dimension (LDH2,N*NB) */
/*             The leading NC-by-N*NB part of this array must contain */
/*             the multivariable sequence M2(k), where k = 1,2,...,N. */
/*             Each parameter M2(k) is an NC-by-NB matrix, whose */
/*             (i,j)-th element must be stored in H2(i,(k-1)*NB+j) for */
/*             i = 1,2,...,NC and j = 1,2,...,NB. */

/*     LDH2    INTEGER */
/*             The leading dimension of array H2.  LDH2 >= MAX(1,NC). */

/*     SS      (output) DOUBLE PRECISION array, dimension (LDSS,NB) */
/*             The leading NC-by-NB part of this array contains the */
/*             matrix SS. */

/*     LDSS    INTEGER */
/*             The leading dimension of array SS.  LDSS >= MAX(1,NC). */

/*     SE      (output) DOUBLE PRECISION array, dimension (LDSE,NB) */
/*             The leading NC-by-NB part of this array contains the */
/*             quadratic error matrix SE. */

/*     LDSE    INTEGER */
/*             The leading dimension of array SE.  LDSE >= MAX(1,NC). */

/*     PRE     (output) DOUBLE PRECISION array, dimension (LDPRE,NB) */
/*             The leading NC-by-NB part of this array contains the */
/*             percentage relative error matrix PRE. */

/*     LDPRE   INTEGER */
/*             The leading dimension of array PRE.  LDPRE >= MAX(1,NC). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in the computation of the error */
/*             matrices SE and PRE. If the user sets TOL to be less than */
/*             EPS then the tolerance is taken as EPS, where EPS is the */
/*             machine precision (see LAPACK Library routine DLAMCH). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The (i,j)-th element of the matrix SS is defined by: */
/*                        N          2 */
/*               SS    = SUM  M1  (k) .                            (1) */
/*                 ij    k=1    ij */

/*     The (i,j)-th element of the quadratic error matrix SE is defined */
/*     by: */
/*                        N                      2 */
/*               SE    = SUM  (M1  (k) - M2  (k)) .                (2) */
/*                 ij    k=1     ij        ij */

/*     The (i,j)-th element of the percentage relative error matrix PRE */
/*     is defined by: */

/*               PRE   = 100 x SQRT( SE  / SS  ).                  (3) */
/*                  ij                 ij    ij */

/*     The following precautions are taken by the routine to guard */
/*     against underflow and overflow: */

/*     (i) if ABS( M1  (k) ) > 1/TOL or ABS( M1  (k) - M2  (k) ) > 1/TOL, */
/*                   ij                        ij        ij */

/*         then SE   and SS   are set to 1/TOL and PRE   is set to 1; and */
/*                ij       ij                         ij */

/*     (ii) if ABS( SS  ) <= TOL, then PRE   is set to 100. */
/*                    ij                  ij */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately */
/*        2xNBxNCx(N+1) multiplications/divisions, */
/*        4xNBxNCxN     additions/subtractions and */
/*          NBxNC       square roots. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB09AD by S. Van Huffel, Katholieke */
/*     University Leuven, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Closeness multivariable sequences, elementary matrix operations, */
/*     real signals, system response. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 174 "SB09MD.f"
    /* Parameter adjustments */
#line 174 "SB09MD.f"
    h1_dim1 = *ldh1;
#line 174 "SB09MD.f"
    h1_offset = 1 + h1_dim1;
#line 174 "SB09MD.f"
    h1 -= h1_offset;
#line 174 "SB09MD.f"
    h2_dim1 = *ldh2;
#line 174 "SB09MD.f"
    h2_offset = 1 + h2_dim1;
#line 174 "SB09MD.f"
    h2 -= h2_offset;
#line 174 "SB09MD.f"
    ss_dim1 = *ldss;
#line 174 "SB09MD.f"
    ss_offset = 1 + ss_dim1;
#line 174 "SB09MD.f"
    ss -= ss_offset;
#line 174 "SB09MD.f"
    se_dim1 = *ldse;
#line 174 "SB09MD.f"
    se_offset = 1 + se_dim1;
#line 174 "SB09MD.f"
    se -= se_offset;
#line 174 "SB09MD.f"
    pre_dim1 = *ldpre;
#line 174 "SB09MD.f"
    pre_offset = 1 + pre_dim1;
#line 174 "SB09MD.f"
    pre -= pre_offset;
#line 174 "SB09MD.f"

#line 174 "SB09MD.f"
    /* Function Body */
#line 174 "SB09MD.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 178 "SB09MD.f"
    if (*n < 0) {
#line 179 "SB09MD.f"
	*info = -1;
#line 180 "SB09MD.f"
    } else if (*nc < 0) {
#line 181 "SB09MD.f"
	*info = -2;
#line 182 "SB09MD.f"
    } else if (*nb < 0) {
#line 183 "SB09MD.f"
	*info = -3;
#line 184 "SB09MD.f"
    } else if (*ldh1 < max(1,*nc)) {
#line 185 "SB09MD.f"
	*info = -5;
#line 186 "SB09MD.f"
    } else if (*ldh2 < max(1,*nc)) {
#line 187 "SB09MD.f"
	*info = -7;
#line 188 "SB09MD.f"
    } else if (*ldss < max(1,*nc)) {
#line 189 "SB09MD.f"
	*info = -9;
#line 190 "SB09MD.f"
    } else if (*ldse < max(1,*nc)) {
#line 191 "SB09MD.f"
	*info = -11;
#line 192 "SB09MD.f"
    } else if (*ldpre < max(1,*nc)) {
#line 193 "SB09MD.f"
	*info = -13;
#line 194 "SB09MD.f"
    }

#line 196 "SB09MD.f"
    if (*info != 0) {

/*        Error return. */

#line 200 "SB09MD.f"
	i__1 = -(*info);
#line 200 "SB09MD.f"
	xerbla_("SB09MD", &i__1, (ftnlen)6);
#line 201 "SB09MD.f"
	return 0;
#line 202 "SB09MD.f"
    }

/*     Quick return if possible. */

#line 206 "SB09MD.f"
    if (*n == 0 || *nc == 0 || *nb == 0) {
#line 206 "SB09MD.f"
	return 0;
#line 206 "SB09MD.f"
    }

/* Computing MAX */
#line 209 "SB09MD.f"
    d__1 = *tol, d__2 = dlamch_("Epsilon", (ftnlen)7);
#line 209 "SB09MD.f"
    toler = max(d__1,d__2);
#line 210 "SB09MD.f"
    epso = 1. / toler;

#line 212 "SB09MD.f"
    i__1 = *nb;
#line 212 "SB09MD.f"
    for (j = 1; j <= i__1; ++j) {

#line 214 "SB09MD.f"
	i__2 = *nc;
#line 214 "SB09MD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 215 "SB09MD.f"
	    sse = 0.;
#line 216 "SB09MD.f"
	    sss = 0.;
#line 217 "SB09MD.f"
	    noflow = TRUE_;
#line 218 "SB09MD.f"
	    k = 0;

/*           WHILE ( ( NOFLOW .AND. ( K .LT. N*NB ) ) DO */
#line 221 "SB09MD.f"
L20:
#line 221 "SB09MD.f"
	    if (noflow && k < *n * *nb) {
#line 222 "SB09MD.f"
		var = h1[i__ + (k + j) * h1_dim1];
#line 223 "SB09MD.f"
		vare = h2[i__ + (k + j) * h2_dim1] - var;
#line 224 "SB09MD.f"
		if (abs(var) > epso || abs(vare) > epso) {
#line 226 "SB09MD.f"
		    se[i__ + j * se_dim1] = epso;
#line 227 "SB09MD.f"
		    ss[i__ + j * ss_dim1] = epso;
#line 228 "SB09MD.f"
		    pre[i__ + j * pre_dim1] = 1.;
#line 229 "SB09MD.f"
		    noflow = FALSE_;
#line 230 "SB09MD.f"
		} else {
#line 231 "SB09MD.f"
		    if (abs(vare) > toler) {
#line 231 "SB09MD.f"
			sse += vare * vare;
#line 231 "SB09MD.f"
		    }
#line 232 "SB09MD.f"
		    if (abs(var) > toler) {
#line 232 "SB09MD.f"
			sss += var * var;
#line 232 "SB09MD.f"
		    }
#line 233 "SB09MD.f"
		    k += *nb;
#line 234 "SB09MD.f"
		}
#line 235 "SB09MD.f"
		goto L20;
#line 236 "SB09MD.f"
	    }
/*           END WHILE 20 */

#line 239 "SB09MD.f"
	    if (noflow) {
#line 240 "SB09MD.f"
		se[i__ + j * se_dim1] = sse;
#line 241 "SB09MD.f"
		ss[i__ + j * ss_dim1] = sss;
#line 242 "SB09MD.f"
		pre[i__ + j * pre_dim1] = 100.;
#line 243 "SB09MD.f"
		if (sss > toler) {
#line 243 "SB09MD.f"
		    pre[i__ + j * pre_dim1] = sqrt(sse / sss) * 100.;
#line 243 "SB09MD.f"
		}
#line 244 "SB09MD.f"
	    }
#line 245 "SB09MD.f"
/* L40: */
#line 245 "SB09MD.f"
	}

#line 247 "SB09MD.f"
/* L60: */
#line 247 "SB09MD.f"
    }

#line 249 "SB09MD.f"
    return 0;
/* *** Last line of SB09MD *** */
} /* sb09md_ */

