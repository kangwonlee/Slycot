#line 1 "MC01TD.f"
/* MC01TD.f -- translated by f2c (version 20100827).
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

#line 1 "MC01TD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b18 = 1.;

/* Subroutine */ int mc01td_(char *dico, integer *dp, doublereal *p, logical *
	stable, integer *nz, doublereal *dwork, integer *iwarn, integer *info,
	 ftnlen dico_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, k, k1, k2;
    static doublereal p1, pk1;
    static logical dicoc;
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer signum;


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

/*     To determine whether or not a given polynomial P(x) with real */
/*     coefficients is stable, either in the continuous-time or discrete- */
/*     time case. */

/*     A polynomial is said to be stable in the continuous-time case */
/*     if all its zeros lie in the left half-plane, and stable in the */
/*     discrete-time case if all its zeros lie inside the unit circle. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Indicates whether the stability test to be applied to */
/*             P(x) is in the continuous-time or discrete-time case as */
/*             follows: */
/*             = 'C':  Continuous-time case; */
/*             = 'D':  Discrete-time case. */

/*     Input/Output Parameters */

/*     DP      (input/output) INTEGER */
/*             On entry, the degree of the polynomial P(x).  DP >= 0. */
/*             On exit, if P(DP+1) = 0.0 on entry, then DP contains the */
/*             index of the highest power of x for which P(DP+1) <> 0.0. */

/*     P       (input) DOUBLE PRECISION array, dimension (DP+1) */
/*             This array must contain the coefficients of P(x) in */
/*             increasing powers of x. */

/*     STABLE  (output) LOGICAL */
/*             Contains the value .TRUE. if P(x) is stable and the value */
/*             .FALSE. otherwise (see also NUMERICAL ASPECTS). */

/*     NZ      (output) INTEGER */
/*             If INFO = 0, contains the number of unstable zeros - that */
/*             is, the number of zeros of P(x) in the right half-plane if */
/*             DICO = 'C' or the number of zeros of P(x) outside the unit */
/*             circle if DICO = 'D' (see also NUMERICAL ASPECTS). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (2*DP+2) */
/*             The leading (DP+1) elements of DWORK contain the Routh */
/*             coefficients, if DICO = 'C', or the constant terms of */
/*             the Schur-Cohn transforms, if DICO = 'D'. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = k:  if the degree of the polynomial P(x) has been */
/*                   reduced to (DB - k) because P(DB+1-j) = 0.0 on entry */
/*                   for j = 0, 1,..., k-1 and P(DB+1-k) <> 0.0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if on entry, P(x) is the zero polynomial; */
/*             = 2:  if the polynomial P(x) is most probably unstable, */
/*                   although it may be stable with one or more zeros */
/*                   very close to either the imaginary axis if */
/*                   DICO = 'C' or the unit circle if DICO = 'D'. */
/*                   The number of unstable zeros (NZ) is not determined. */

/*     METHOD */

/*     The stability of the real polynomial */
/*                                         2                DP */
/*        P(x) = p(0) + p(1) * x + p(2) * x  + ... + p(DP) x */

/*     is determined as follows. */

/*     In the continuous-time case (DICO = 'C') the Routh algorithm */
/*     (see [1]) is used. The routine computes the Routh coefficients and */
/*     if they are non-zero then the number of sign changes in the */
/*     sequence of the coefficients is equal to the number of zeros with */
/*     positive imaginary part. */

/*     In the discrete-time case (DICO = 'D') the Schur-Cohn */
/*     algorithm (see [2] and [3]) is applied to the reciprocal */
/*     polynomial */
/*                                                2               DP */
/*        Q(x) = p(DP) + p(DP-1) * x + p(DP-2) * x  + ... + p(0) x  . */

/*     The routine computes the constant terms of the Schur transforms */
/*     and if all of them are non-zero then the number of zeros of P(x) */
/*     with modulus greater than unity is obtained from the sequence of */
/*     constant terms. */

/*     REFERENCES */

/*     [1] Gantmacher, F.R. */
/*         Applications of the Theory of Matrices. */
/*         Interscience Publishers, New York, 1959. */

/*     [2] Kucera, V. */
/*         Discrete Linear Control. The Algorithmic Approach. */
/*         John Wiley & Sons, Chichester, 1979. */

/*     [3] Henrici, P. */
/*         Applied and Computational Complex Analysis (Vol. 1). */
/*         John Wiley & Sons, New York, 1974. */

/*     NUMERICAL ASPECTS */

/*     The algorithm used by the routine is numerically stable. */

/*     Note that if some of the Routh coefficients (DICO = 'C') or */
/*     some of the constant terms of the Schur-Cohn transforms (DICO = */
/*     'D') are small relative to EPS (the machine precision), then */
/*     the number of unstable zeros (and hence the value of STABLE) may */
/*     be incorrect. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01HD by F. Delebecque and */
/*     A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary polynomial operations, polynomial operations, */
/*     stability, stability criteria, zeros. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 180 "MC01TD.f"
    /* Parameter adjustments */
#line 180 "MC01TD.f"
    --dwork;
#line 180 "MC01TD.f"
    --p;
#line 180 "MC01TD.f"

#line 180 "MC01TD.f"
    /* Function Body */
#line 180 "MC01TD.f"
    *iwarn = 0;
#line 181 "MC01TD.f"
    *info = 0;
#line 182 "MC01TD.f"
    dicoc = lsame_(dico, "C", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 186 "MC01TD.f"
    if (! dicoc && ! lsame_(dico, "D", (ftnlen)1, (ftnlen)1)) {
#line 187 "MC01TD.f"
	*info = -1;
#line 188 "MC01TD.f"
    } else if (*dp < 0) {
#line 189 "MC01TD.f"
	*info = -2;
#line 190 "MC01TD.f"
    }

#line 192 "MC01TD.f"
    if (*info != 0) {

/*        Error return. */

#line 196 "MC01TD.f"
	i__1 = -(*info);
#line 196 "MC01TD.f"
	xerbla_("MC01TD", &i__1, (ftnlen)6);
#line 197 "MC01TD.f"
	return 0;
#line 198 "MC01TD.f"
    }

/*     WHILE (DP >= 0 and P(DP+1) = 0 ) DO */
#line 201 "MC01TD.f"
L20:
#line 201 "MC01TD.f"
    if (*dp >= 0) {
#line 202 "MC01TD.f"
	if (p[*dp + 1] == 0.) {
#line 203 "MC01TD.f"
	    --(*dp);
#line 204 "MC01TD.f"
	    ++(*iwarn);
#line 205 "MC01TD.f"
	    goto L20;
#line 206 "MC01TD.f"
	}
#line 207 "MC01TD.f"
    }
/*     END WHILE 20 */

#line 210 "MC01TD.f"
    if (*dp == -1) {
#line 211 "MC01TD.f"
	*info = 1;
#line 212 "MC01TD.f"
	return 0;
#line 213 "MC01TD.f"
    }

/*     P(x) is not the zero polynomial and its degree is exactly DP. */

#line 217 "MC01TD.f"
    if (dicoc) {

/*        Continuous-time case. */

/*        Compute the Routh coefficients and the number of sign changes. */

#line 223 "MC01TD.f"
	i__1 = *dp + 1;
#line 223 "MC01TD.f"
	dcopy_(&i__1, &p[1], &c__1, &dwork[1], &c__1);
#line 224 "MC01TD.f"
	*nz = 0;
#line 225 "MC01TD.f"
	k = *dp;
/*        WHILE ( K > 0 and DWORK(K) non-zero) DO */
#line 227 "MC01TD.f"
L40:
#line 227 "MC01TD.f"
	if (k > 0) {
#line 228 "MC01TD.f"
	    if (dwork[k] == 0.) {
#line 229 "MC01TD.f"
		*info = 2;
#line 230 "MC01TD.f"
	    } else {
#line 231 "MC01TD.f"
		alpha = dwork[k + 1] / dwork[k];
#line 232 "MC01TD.f"
		if (alpha < 0.) {
#line 232 "MC01TD.f"
		    ++(*nz);
#line 232 "MC01TD.f"
		}
#line 233 "MC01TD.f"
		--k;

#line 235 "MC01TD.f"
		for (i__ = k; i__ >= 2; i__ += -2) {
#line 236 "MC01TD.f"
		    dwork[i__] -= alpha * dwork[i__ - 1];
#line 237 "MC01TD.f"
/* L60: */
#line 237 "MC01TD.f"
		}

#line 239 "MC01TD.f"
		goto L40;
#line 240 "MC01TD.f"
	    }
#line 241 "MC01TD.f"
	}
/*        END WHILE 40 */
#line 243 "MC01TD.f"
    } else {

/*        Discrete-time case. */

/*        To apply [3], section 6.8, on the reciprocal of polynomial */
/*        P(x) the elements of the array P are copied in DWORK in */
/*        reverse order. */

#line 251 "MC01TD.f"
	i__1 = *dp + 1;
#line 251 "MC01TD.f"
	dcopy_(&i__1, &p[1], &c__1, &dwork[1], &c_n1);
/*                                                           K-1 */
/*        DWORK(K),...,DWORK(DP+1), are the coefficients of T   P(x) */
/*        scaled with a factor alpha(K) in order to avoid over- or */
/*        underflow, */
/*                                                    i-1 */
/*        DWORK(i), i = 1,...,K, contains alpha(i) * T   P(0). */

#line 259 "MC01TD.f"
	signum = 1;
#line 260 "MC01TD.f"
	*nz = 0;
#line 261 "MC01TD.f"
	k = 1;
/*        WHILE ( K <= DP and DWORK(K) non-zero ) DO */
#line 263 "MC01TD.f"
L80:
#line 263 "MC01TD.f"
	if (k <= *dp && *info == 0) {
/*                                        K */
/*           Compute the coefficients of T P(x). */

#line 267 "MC01TD.f"
	    k1 = *dp - k + 2;
#line 268 "MC01TD.f"
	    k2 = *dp + 2;
#line 269 "MC01TD.f"
	    alpha = dwork[k - 1 + idamax_(&k1, &dwork[k], &c__1)];
#line 270 "MC01TD.f"
	    if (alpha == 0.) {
#line 271 "MC01TD.f"
		*info = 2;
#line 272 "MC01TD.f"
	    } else {
#line 273 "MC01TD.f"
		dcopy_(&k1, &dwork[k], &c__1, &dwork[k2], &c__1);
#line 274 "MC01TD.f"
		drscl_(&k1, &alpha, &dwork[k2], &c__1);
#line 275 "MC01TD.f"
		p1 = dwork[k2];
#line 276 "MC01TD.f"
		pk1 = dwork[k2 + k1 - 1];

#line 278 "MC01TD.f"
		i__1 = k1 - 1;
#line 278 "MC01TD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 279 "MC01TD.f"
		    dwork[k + i__] = p1 * dwork[*dp + 1 + i__] - pk1 * dwork[
			    k2 + k1 - i__];
#line 280 "MC01TD.f"
/* L100: */
#line 280 "MC01TD.f"
		}

/*              Compute the number of unstable zeros. */

#line 284 "MC01TD.f"
		++k;
#line 285 "MC01TD.f"
		if (dwork[k] == 0.) {
#line 286 "MC01TD.f"
		    *info = 2;
#line 287 "MC01TD.f"
		} else {
#line 288 "MC01TD.f"
		    signum = (integer) (signum * d_sign(&c_b18, &dwork[k]));
#line 289 "MC01TD.f"
		    if ((doublereal) signum < 0.) {
#line 289 "MC01TD.f"
			++(*nz);
#line 289 "MC01TD.f"
		    }
#line 290 "MC01TD.f"
		}
#line 291 "MC01TD.f"
		goto L80;
#line 292 "MC01TD.f"
	    }
/*           END WHILE 80 */
#line 294 "MC01TD.f"
	}
#line 295 "MC01TD.f"
    }

#line 297 "MC01TD.f"
    if (*info == 0 && *nz == 0) {
#line 298 "MC01TD.f"
	*stable = TRUE_;
#line 299 "MC01TD.f"
    } else {
#line 300 "MC01TD.f"
	*stable = FALSE_;
#line 301 "MC01TD.f"
    }

#line 303 "MC01TD.f"
    return 0;
/* *** Last line of MC01TD *** */
} /* mc01td_ */

