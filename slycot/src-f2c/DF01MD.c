#line 1 "DF01MD.f"
/* DF01MD.f -- translated by f2c (version 20100827).
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

#line 1 "DF01MD.f"
/* Subroutine */ int df01md_(char *sico, integer *n, doublereal *dt, 
	doublereal *a, doublereal *dwork, integer *info, ftnlen sico_len)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double atan(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, m;
    static doublereal a0;
    static integer i2;
    static doublereal w1, w2, w3;
    static integer md2, ind1, ind2;
    static logical lsig;
    extern /* Subroutine */ int dg01nd_(char *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical lsico;
    static doublereal pibym;
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

/*     To compute the sine transform or cosine transform of a real */
/*     signal. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SICO    CHARACTER*1 */
/*             Indicates whether the sine transform or cosine transform */
/*             is to be computed as follows: */
/*             = 'S':  The sine transform is computed; */
/*             = 'C':  The cosine transform is computed. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of samples.  N must be a power of 2 plus 1. */
/*             N >= 5. */

/*     DT      (input) DOUBLE PRECISION */
/*             The sampling time of the signal. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the signal to be */
/*             processed. */
/*             On exit, this array contains either the sine transform, if */
/*             SICO = 'S', or the cosine transform, if SICO = 'C', of the */
/*             given signal. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N+1) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Let A(1), A(2),..., A(N) be a real signal of N samples. */

/*     If SICO = 'S', the routine computes the sine transform of A as */
/*     follows. First, transform A(i), i = 1,2,...,N, into the complex */
/*     signal B(i), i = 1,2,...,(N+1)/2, where */

/*        B(1) = -2*A(2), */
/*        B(i) = {A(2i-2) - A(2i)} - j*A(2i-1) for i = 2,3,...,(N-1)/2, */
/*        B((N+1)/2) = 2*A(N-1) and j**2 = -1. */

/*     Next, perform a discrete inverse Fourier transform on B(i) by */
/*     calling SLICOT Library Routine DG01ND, to give the complex signal */
/*     Z(i), i = 1,2,...,(N-1)/2, from which the real signal C(i) may be */
/*     obtained as follows: */

/*        C(2i-1) = Re(Z(i)),  C(2i) = Im(Z(i)) for i = 1,2,...,(N-1)/2. */

/*     Finally, compute the sine transform coefficients S ,S ,...,S */
/*                                                       1  2      N */
/*     given by */

/*        S  = 0, */
/*         1 */
/*                {                     [C(k) + C(N+1-k)]     } */
/*        S  = DT*{[C(k) - C(N+1-k)] - -----------------------}, */
/*         k      {                    [2*sin(pi*(k-1)/(N-1))]} */

/*           for k = 2,3,...,N-1, and */

/*        S = 0. */
/*         N */

/*     If SICO = 'C', the routine computes the cosine transform of A as */
/*     follows. First, transform A(i), i = 1,2,...,N, into the complex */
/*     signal B(i), i = 1,2,...,(N+1)/2, where */

/*        B(1) = 2*A(1), */
/*        B(i) = 2*A(2i-1) + 2*j*{[A(2i-2) - A(2i)]} */
/*        for i = 2,3,...,(N-1)/2 and B((N+1)/2) = 2*A(N). */

/*     Next, perform a discrete inverse Fourier transform on B(i) by */
/*     calling SLICOT Library Routine DG01ND, to give the complex signal */
/*     Z(i), i = 1,2,...,(N-1)/2, from which the real signal D(i) may be */
/*     obtained as follows: */

/*        D(2i-1) = Re(Z(i)),  D(2i) = Im(Z(i)) for i = 1,2,...,(N-1)/2. */

/*     Finally, compute the cosine transform coefficients S ,S ,...,S */
/*                                                         1  2      N */
/*     given by */

/*        S  = 2*DT*[D(1) + A0], */
/*         1 */
/*                {                     [D(k) - D(N+1-k)]     } */
/*        S  = DT*{[D(k) + D(N+1-k)] - -----------------------}, */
/*         k      {                    [2*sin(pi*(k-1)/(N-1))]} */


/*           for k = 2,3,...,N-1, and */

/*        S  = 2*DT*[D(1) - A0], */
/*         N */
/*                 (N-1)/2 */
/*     where A0 = 2*SUM   A(2i). */
/*                  i=1 */

/*     REFERENCES */

/*     [1] Rabiner, L.R. and Rader, C.M. */
/*         Digital Signal Processing. */
/*         IEEE Press, 1972. */

/*     [2] Oppenheim, A.V. and Schafer, R.W. */
/*         Discrete-Time Signal Processing. */
/*         Prentice-Hall Signal Processing Series, 1989. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires 0( N*log(N) ) operations. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Supersedes Release 2.0 routine DF01AD by F. Dumortier, and */
/*     R.M.C. Dekeyser, State University of Gent, Belgium. */

/*     REVISIONS */

/*     V. Sima, Jan. 2003. */

/*     KEYWORDS */

/*     Digital signal processing, fast Fourier transform, complex */
/*     signals. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 185 "DF01MD.f"
    /* Parameter adjustments */
#line 185 "DF01MD.f"
    --dwork;
#line 185 "DF01MD.f"
    --a;
#line 185 "DF01MD.f"

#line 185 "DF01MD.f"
    /* Function Body */
#line 185 "DF01MD.f"
    *info = 0;
#line 186 "DF01MD.f"
    lsico = lsame_(sico, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 190 "DF01MD.f"
    if (! lsico && ! lsame_(sico, "C", (ftnlen)1, (ftnlen)1)) {
#line 191 "DF01MD.f"
	*info = -1;
#line 192 "DF01MD.f"
    } else {
#line 193 "DF01MD.f"
	m = 0;
#line 194 "DF01MD.f"
	if (*n > 4) {
#line 195 "DF01MD.f"
	    m = *n - 1;
/*           WHILE ( MOD( M, 2 ).EQ.0 ) DO */
#line 197 "DF01MD.f"
L10:
#line 198 "DF01MD.f"
	    if (m % 2 == 0) {
#line 199 "DF01MD.f"
		m /= 2;
#line 200 "DF01MD.f"
		goto L10;
#line 201 "DF01MD.f"
	    }
/*           END WHILE 10 */
#line 203 "DF01MD.f"
	}
#line 204 "DF01MD.f"
	if (m != 1) {
#line 204 "DF01MD.f"
	    *info = -2;
#line 204 "DF01MD.f"
	}
#line 205 "DF01MD.f"
    }

#line 207 "DF01MD.f"
    if (*info != 0) {

/*        Error return. */

#line 211 "DF01MD.f"
	i__1 = -(*info);
#line 211 "DF01MD.f"
	xerbla_("DF01MD", &i__1, (ftnlen)6);
#line 212 "DF01MD.f"
	return 0;
#line 213 "DF01MD.f"
    }

/*     Initialisation. */

#line 217 "DF01MD.f"
    m = *n - 1;
#line 218 "DF01MD.f"
    md2 = (*n + 1) / 2;
#line 219 "DF01MD.f"
    pibym = atan(1.) * 4. / (doublereal) m;
#line 220 "DF01MD.f"
    i2 = 1;
#line 221 "DF01MD.f"
    dwork[md2 + 1] = 0.;
#line 222 "DF01MD.f"
    dwork[md2 * 2] = 0.;

#line 224 "DF01MD.f"
    if (lsico) {

/*        Sine transform. */

#line 228 "DF01MD.f"
	lsig = TRUE_;
#line 229 "DF01MD.f"
	dwork[1] = a[2] * -2.;
#line 230 "DF01MD.f"
	dwork[md2] = a[m] * 2.;

#line 232 "DF01MD.f"
	i__1 = m;
#line 232 "DF01MD.f"
	for (i__ = 4; i__ <= i__1; i__ += 2) {
#line 233 "DF01MD.f"
	    ++i2;
#line 234 "DF01MD.f"
	    dwork[i2] = a[i__ - 2] - a[i__];
#line 235 "DF01MD.f"
	    dwork[md2 + i2] = -a[i__ - 1];
#line 236 "DF01MD.f"
/* L20: */
#line 236 "DF01MD.f"
	}

#line 238 "DF01MD.f"
    } else {

/*        Cosine transform. */

#line 242 "DF01MD.f"
	lsig = FALSE_;
#line 243 "DF01MD.f"
	dwork[1] = a[1] * 2.;
#line 244 "DF01MD.f"
	dwork[md2] = a[*n] * 2.;
#line 245 "DF01MD.f"
	a0 = a[2];

#line 247 "DF01MD.f"
	i__1 = m;
#line 247 "DF01MD.f"
	for (i__ = 4; i__ <= i__1; i__ += 2) {
#line 248 "DF01MD.f"
	    ++i2;
#line 249 "DF01MD.f"
	    dwork[i2] = a[i__ - 1] * 2.;
#line 250 "DF01MD.f"
	    dwork[md2 + i2] = (a[i__ - 2] - a[i__]) * 2.;
#line 251 "DF01MD.f"
	    a0 += a[i__];
#line 252 "DF01MD.f"
/* L30: */
#line 252 "DF01MD.f"
	}

#line 254 "DF01MD.f"
	a0 *= 2.;
#line 255 "DF01MD.f"
    }

/*     Inverse Fourier transform. */

#line 259 "DF01MD.f"
    i__1 = md2 - 1;
#line 259 "DF01MD.f"
    dg01nd_("Inverse", &i__1, &dwork[1], &dwork[md2 + 1], info, (ftnlen)7);

/*     Sine or cosine coefficients. */

#line 263 "DF01MD.f"
    if (lsico) {
#line 264 "DF01MD.f"
	a[1] = 0.;
#line 265 "DF01MD.f"
	a[*n] = 0.;
#line 266 "DF01MD.f"
    } else {
#line 267 "DF01MD.f"
	a[1] = *dt * 2. * (dwork[1] + a0);
#line 268 "DF01MD.f"
	a[*n] = *dt * 2. * (dwork[1] - a0);
#line 269 "DF01MD.f"
    }

#line 271 "DF01MD.f"
    ind1 = md2 + 1;
#line 272 "DF01MD.f"
    ind2 = *n;

#line 274 "DF01MD.f"
    i__1 = m - 1;
#line 274 "DF01MD.f"
    for (i__ = 1; i__ <= i__1; i__ += 2) {
#line 275 "DF01MD.f"
	w1 = dwork[ind1];
#line 276 "DF01MD.f"
	w2 = dwork[ind2];
#line 277 "DF01MD.f"
	if (lsig) {
#line 277 "DF01MD.f"
	    w2 = -w2;
#line 277 "DF01MD.f"
	}
#line 278 "DF01MD.f"
	w3 = sin(pibym * (doublereal) i__) * 2.;
#line 279 "DF01MD.f"
	a[i__ + 1] = *dt * (w1 + w2 - (w1 - w2) / w3);
#line 280 "DF01MD.f"
	++ind1;
#line 281 "DF01MD.f"
	--ind2;
#line 282 "DF01MD.f"
/* L40: */
#line 282 "DF01MD.f"
    }

#line 284 "DF01MD.f"
    ind1 = 2;
#line 285 "DF01MD.f"
    ind2 = md2 - 1;

#line 287 "DF01MD.f"
    i__1 = m - 2;
#line 287 "DF01MD.f"
    for (i__ = 2; i__ <= i__1; i__ += 2) {
#line 288 "DF01MD.f"
	w1 = dwork[ind1];
#line 289 "DF01MD.f"
	w2 = dwork[ind2];
#line 290 "DF01MD.f"
	if (lsig) {
#line 290 "DF01MD.f"
	    w2 = -w2;
#line 290 "DF01MD.f"
	}
#line 291 "DF01MD.f"
	w3 = sin(pibym * (doublereal) i__) * 2.;
#line 292 "DF01MD.f"
	a[i__ + 1] = *dt * (w1 + w2 - (w1 - w2) / w3);
#line 293 "DF01MD.f"
	++ind1;
#line 294 "DF01MD.f"
	--ind2;
#line 295 "DF01MD.f"
/* L50: */
#line 295 "DF01MD.f"
    }

#line 297 "DF01MD.f"
    return 0;
/* *** Last line of DF01MD *** */
} /* df01md_ */

