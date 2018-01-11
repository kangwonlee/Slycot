#line 1 "FD01AD.f"
/* FD01AD.f -- translated by f2c (version 20100827).
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

#line 1 "FD01AD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int fd01ad_(char *jp, integer *l, doublereal *lambda, 
	doublereal *xin, doublereal *yin, doublereal *efor, doublereal *xf, 
	doublereal *epsbck, doublereal *cteta, doublereal *steta, doublereal *
	yq, doublereal *epos, doublereal *eout, doublereal *salph, integer *
	iwarn, integer *info, ftnlen jp_len)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal xfi, eps, yqi;
    static logical both;
    static doublereal temp, norm;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal fnode;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal ctemp;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), xerbla_(char *, 
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

/*     To solve the least-squares filtering problem recursively in time. */
/*     Each subroutine call implements one time update of the solution. */
/*     The algorithm uses a fast QR-decomposition based approach. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JP      CHARACTER*1 */
/*             Indicates whether the user wishes to apply both prediction */
/*             and filtering parts, as follows: */
/*             = 'B':  Both prediction and filtering parts are to be */
/*                     applied; */
/*             = 'P':  Only the prediction section is to be applied. */

/*     Input/Output Parameters */

/*     L       (input) INTEGER */
/*             The length of the impulse response of the equivalent */
/*             transversal filter model.  L >= 1. */

/*     LAMBDA  (input) DOUBLE PRECISION */
/*             Square root of the forgetting factor. */
/*             For tracking capabilities and exponentially stable error */
/*             propagation, LAMBDA < 1.0 (strict inequality) should */
/*             be used.  0.0 < LAMBDA <= 1.0. */

/*     XIN     (input) DOUBLE PRECISION */
/*             The input sample at instant n. */
/*             (The situation just before and just after the call of */
/*             the routine are denoted by instant (n-1) and instant n, */
/*             respectively.) */

/*     YIN     (input) DOUBLE PRECISION */
/*             If JP = 'B', then YIN must contain the reference sample */
/*             at instant n. */
/*             Otherwise, YIN is not referenced. */

/*     EFOR    (input/output) DOUBLE PRECISION */
/*             On entry, this parameter must contain the square root of */
/*             exponentially weighted forward prediction error energy */
/*             at instant (n-1).  EFOR >= 0.0. */
/*             On exit, this parameter contains the square root of the */
/*             exponentially weighted forward prediction error energy */
/*             at instant n. */

/*     XF      (input/output) DOUBLE PRECISION array, dimension (L) */
/*             On entry, this array must contain the transformed forward */
/*             prediction variables at instant (n-1). */
/*             On exit, this array contains the transformed forward */
/*             prediction variables at instant n. */

/*     EPSBCK  (input/output) DOUBLE PRECISION array, dimension (L+1) */
/*             On entry, the leading L elements of this array must */
/*             contain the normalized a posteriori backward prediction */
/*             error residuals of orders zero through L-1, respectively, */
/*             at instant (n-1), and EPSBCK(L+1) must contain the */
/*             square-root of the so-called "conversion factor" at */
/*             instant (n-1). */
/*             On exit, this array contains the normalized a posteriori */
/*             backward prediction error residuals, plus the square root */
/*             of the conversion factor at instant n. */

/*     CTETA   (input/output) DOUBLE PRECISION array, dimension (L) */
/*             On entry, this array must contain the cosines of the */
/*             rotation angles used in time updates, at instant (n-1). */
/*             On exit, this array contains the cosines of the rotation */
/*             angles at instant n. */

/*     STETA   (input/output) DOUBLE PRECISION array, dimension (L) */
/*             On entry, this array must contain the sines of the */
/*             rotation angles used in time updates, at instant (n-1). */
/*             On exit, this array contains the sines of the rotation */
/*             angles at instant n. */

/*     YQ      (input/output) DOUBLE PRECISION array, dimension (L) */
/*             On entry, if JP = 'B', then this array must contain the */
/*             orthogonally transformed reference vector at instant */
/*             (n-1). These elements are also the tap multipliers of an */
/*             equivalent normalized lattice least-squares filter. */
/*             Otherwise, YQ is not referenced and can be supplied as */
/*             a dummy array (i.e., declare this array to be YQ(1) in */
/*             the calling program). */
/*             On exit, if JP = 'B', then this array contains the */
/*             orthogonally transformed reference vector at instant n. */

/*     EPOS    (output) DOUBLE PRECISION */
/*             The a posteriori forward prediction error residual. */

/*     EOUT    (output) DOUBLE PRECISION */
/*             If JP = 'B', then EOUT contains the a posteriori output */
/*             error residual from the least-squares filter at instant n. */

/*     SALPH   (output) DOUBLE PRECISION array, dimension (L) */
/*             The element SALPH(i), i=1,...,L, contains the opposite of */
/*             the i-(th) reflection coefficient for the least-squares */
/*             normalized lattice predictor (whose value is -SALPH(i)). */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  an element to be annihilated by a rotation is less */
/*                   than the machine precision (see LAPACK Library */
/*                   routine DLAMCH). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The output error EOUT at instant n, denoted by EOUT(n), is the */
/*     reference sample minus a linear combination of L successive input */
/*     samples: */

/*                           L-1 */
/*        EOUT(n) = YIN(n) - SUM h_i * XIN(n-i), */
/*                           i=0 */

/*     where YIN(n) and XIN(n) are the scalar samples at instant n. */
/*     A least-squares filter uses those h_0,...,h_{L-1} which minimize */
/*     an exponentially weighted sum of successive output errors squared: */

/*         n */
/*        SUM [LAMBDA**(2(n-k)) * EOUT(k)**2]. */
/*        k=1 */

/*     Each subroutine call performs a time update of the least-squares */
/*     filter using a fast least-squares algorithm derived from a */
/*     QR decomposition, as described in references [1] and [2] (the */
/*     notation from [2] is followed in the naming of the arrays). */
/*     The algorithm does not compute the parameters h_0,...,h_{L-1} from */
/*     the above formula, but instead furnishes the parameters of an */
/*     equivalent normalized least-squares lattice filter, which are */
/*     available from the arrays SALPH (reflection coefficients) and YQ */
/*     (tap multipliers), as well as the exponentially weighted input */
/*     signal energy */

/*         n                                              L */
/*        SUM [LAMBDA**(2(n-k)) * XIN(k)**2] = EFOR**2 + SUM XF(i)**2. */
/*        k=1                                            i=1 */

/*     For more details on reflection coefficients and tap multipliers, */
/*     references [2] and [4] are recommended. */

/*     REFERENCES */

/*     [1]  Proudler, I. K., McWhirter, J. G., and Shepherd, T. J. */
/*          Fast QRD based algorithms for least-squares linear */
/*          prediction. */
/*          Proceedings IMA Conf. Mathematics in Signal Processing */
/*          Warwick, UK, December 1988. */

/*     [2]  Regalia, P. A., and Bellanger, M. G. */
/*          On the duality between QR methods and lattice methods in */
/*          least-squares adaptive filtering. */
/*          IEEE Trans. Signal Processing, SP-39, pp. 879-891, */
/*          April 1991. */

/*     [3]  Regalia, P. A. */
/*          Numerical stability properties of a QR-based fast */
/*          least-squares algorithm. */
/*          IEEE Trans. Signal Processing, SP-41, June 1993. */

/*     [4]  Lev-Ari, H., Kailath, T., and Cioffi, J. */
/*          Least-squares adaptive lattice and transversal filters: */
/*          A unified geometric theory. */
/*          IEEE Trans. Information Theory, IT-30, pp. 222-236, */
/*          March 1984. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O(L) operations for each subroutine call. */
/*     It is backward consistent for all input sequences XIN, and */
/*     backward stable for persistently exciting input sequences, */
/*     assuming LAMBDA < 1.0 (see [3]). */
/*     If the condition of the signal is very poor (IWARN = 1), then the */
/*     results are not guaranteed to be reliable. */

/*     FURTHER COMMENTS */

/*     1.  For tracking capabilities and exponentially stable error */
/*         propagation, LAMBDA < 1.0 should be used.  LAMBDA is typically */
/*         chosen slightly less than 1.0 so that "past" data are */
/*         exponentially forgotten. */
/*     2.  Prior to the first subroutine call, the variables must be */
/*         initialized. The following initial values are recommended: */

/*         XF(i) = 0.0,        i=1,...,L */
/*         EPSBCK(i) = 0.0     i=1,...,L */
/*         EPSBCK(L+1) = 1.0 */
/*         CTETA(i) = 1.0      i=1,...,L */
/*         STETA(i) = 0.0      i=1,...,L */
/*         YQ(i) = 0.0         i=1,...,L */

/*         EFOR = 0.0          (exact start) */
/*         EFOR = "small positive constant" (soft start). */

/*         Soft starts are numerically more reliable, but result in a */
/*         biased least-squares solution during the first few iterations. */
/*         This bias decays exponentially fast provided LAMBDA < 1.0. */
/*         If sigma is the standard deviation of the input sequence */
/*         XIN, then initializing EFOR = sigma*1.0E-02 usually works */
/*         well. */

/*     CONTRIBUTOR */

/*     P. A. Regalia (October 1994). */
/*     Release 4.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1999. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Kalman filtering, least-squares estimator, optimal filtering, */
/*     orthogonal transformation, recursive estimation, QR decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions */
/*     .. Executable statements .. */

/*     Test the input scalar arguments. */

#line 276 "FD01AD.f"
    /* Parameter adjustments */
#line 276 "FD01AD.f"
    --salph;
#line 276 "FD01AD.f"
    --yq;
#line 276 "FD01AD.f"
    --steta;
#line 276 "FD01AD.f"
    --cteta;
#line 276 "FD01AD.f"
    --epsbck;
#line 276 "FD01AD.f"
    --xf;
#line 276 "FD01AD.f"

#line 276 "FD01AD.f"
    /* Function Body */
#line 276 "FD01AD.f"
    both = lsame_(jp, "B", (ftnlen)1, (ftnlen)1);
#line 277 "FD01AD.f"
    *iwarn = 0;
#line 278 "FD01AD.f"
    *info = 0;

#line 280 "FD01AD.f"
    if (! both && ! lsame_(jp, "P", (ftnlen)1, (ftnlen)1)) {
#line 281 "FD01AD.f"
	*info = -1;
#line 282 "FD01AD.f"
    } else if (*l < 1) {
#line 283 "FD01AD.f"
	*info = -2;
#line 284 "FD01AD.f"
    } else if (*lambda <= 0. || *lambda > 1.) {
#line 285 "FD01AD.f"
	*info = -3;
#line 286 "FD01AD.f"
    }

#line 288 "FD01AD.f"
    if (*info != 0) {

/*        Error return. */

#line 292 "FD01AD.f"
	i__1 = -(*info);
#line 292 "FD01AD.f"
	xerbla_("FD01AD", &i__1, (ftnlen)6);
#line 293 "FD01AD.f"
	return 0;
#line 294 "FD01AD.f"
    }

/*     Computation of the machine precision EPS. */

#line 298 "FD01AD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);

/*     Forward prediction rotations. */

#line 302 "FD01AD.f"
    fnode = *xin;

#line 304 "FD01AD.f"
    i__1 = *l;
#line 304 "FD01AD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 305 "FD01AD.f"
	xfi = xf[i__] * *lambda;
#line 306 "FD01AD.f"
	xf[i__] = steta[i__] * fnode + cteta[i__] * xfi;
#line 307 "FD01AD.f"
	fnode = cteta[i__] * fnode - steta[i__] * xfi;
#line 308 "FD01AD.f"
/* L10: */
#line 308 "FD01AD.f"
    }

#line 310 "FD01AD.f"
    *epos = fnode * epsbck[*l + 1];

/*     Update the square root of the prediction energy. */

#line 314 "FD01AD.f"
    *efor *= *lambda;
#line 315 "FD01AD.f"
    temp = dlapy2_(&fnode, efor);
#line 316 "FD01AD.f"
    if (temp < eps) {
#line 317 "FD01AD.f"
	fnode = 0.;
#line 318 "FD01AD.f"
	*iwarn = 1;
#line 319 "FD01AD.f"
    } else {
#line 320 "FD01AD.f"
	fnode = fnode * epsbck[*l + 1] / temp;
#line 321 "FD01AD.f"
    }
#line 322 "FD01AD.f"
    *efor = temp;

/*     Calculate the reflection coefficients and the backward prediction */
/*     errors. */

#line 327 "FD01AD.f"
    for (i__ = *l; i__ >= 1; --i__) {
#line 328 "FD01AD.f"
	if ((d__1 = xf[i__], abs(d__1)) < eps) {
#line 328 "FD01AD.f"
	    *iwarn = 1;
#line 328 "FD01AD.f"
	}
#line 330 "FD01AD.f"
	dlartg_(&temp, &xf[i__], &ctemp, &salph[i__], &norm);
#line 331 "FD01AD.f"
	epsbck[i__ + 1] = ctemp * epsbck[i__] - salph[i__] * fnode;
#line 332 "FD01AD.f"
	fnode = ctemp * fnode + salph[i__] * epsbck[i__];
#line 333 "FD01AD.f"
	temp = norm;
#line 334 "FD01AD.f"
/* L20: */
#line 334 "FD01AD.f"
    }

#line 336 "FD01AD.f"
    epsbck[1] = fnode;

/*     Update to new rotation angles. */

#line 340 "FD01AD.f"
    norm = dnrm2_(l, &epsbck[1], &c__1);
#line 341 "FD01AD.f"
    temp = sqrt((norm + 1.) * (1. - norm));
#line 342 "FD01AD.f"
    epsbck[*l + 1] = temp;

#line 344 "FD01AD.f"
    for (i__ = *l; i__ >= 1; --i__) {
#line 345 "FD01AD.f"
	if ((d__1 = epsbck[i__], abs(d__1)) < eps) {
#line 345 "FD01AD.f"
	    *iwarn = 1;
#line 345 "FD01AD.f"
	}
#line 347 "FD01AD.f"
	dlartg_(&temp, &epsbck[i__], &cteta[i__], &steta[i__], &norm);
#line 348 "FD01AD.f"
	temp = norm;
#line 349 "FD01AD.f"
/* L30: */
#line 349 "FD01AD.f"
    }

/*     Joint process section. */

#line 353 "FD01AD.f"
    if (both) {
#line 354 "FD01AD.f"
	fnode = *yin;

#line 356 "FD01AD.f"
	i__1 = *l;
#line 356 "FD01AD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 357 "FD01AD.f"
	    yqi = yq[i__] * *lambda;
#line 358 "FD01AD.f"
	    yq[i__] = steta[i__] * fnode + cteta[i__] * yqi;
#line 359 "FD01AD.f"
	    fnode = cteta[i__] * fnode - steta[i__] * yqi;
#line 360 "FD01AD.f"
/* L40: */
#line 360 "FD01AD.f"
	}

#line 362 "FD01AD.f"
	*eout = fnode * epsbck[*l + 1];
#line 363 "FD01AD.f"
    }

#line 365 "FD01AD.f"
    return 0;
/* *** Last line of FD01AD *** */
} /* fd01ad_ */

