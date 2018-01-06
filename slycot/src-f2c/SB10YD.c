#line 1 "SB10YD.f"
/* SB10YD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10YD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__4096 = 4096;
static integer c__2048 = 2048;
static doublereal c_b21 = 2.44140625e-4;
static doublereal c_b48 = 0.;
static doublereal c_b49 = 1.;

/* Subroutine */ int sb10yd_(integer *discfl, integer *flag__, integer *
	lendat, doublereal *rfrdat, doublereal *ifrdat, doublereal *omega, 
	integer *n, doublereal *a, integer *lda, doublereal *b, doublereal *
	c__, doublereal *d__, doublereal *tol, integer *iwork, doublereal *
	dwork, integer *ldwork, doublecomplex *zwork, integer *lzwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double atan(doublereal), sqrt(doublereal), acos(doublereal), log(
	    doublereal), exp(doublereal), cos(doublereal), sin(doublereal), 
	    d_imag(doublecomplex *);

    /* Local variables */
    static integer i__, k, p, n1, n2;
    static doublereal p1, p2;
    static integer ii;
    static doublereal pi;
    static integer mn;
    static doublereal pw;
    static integer ip1, ip2, lw1, lw2, lw3, lw4;
    static doublereal rat;
    static integer iws, iwa0, iwab, rank;
    static doublereal tolb;
    static integer iwbp;
    static doublecomplex xhat[1024];
    static integer iwbx;
    static doublereal toll;
    static integer iwxi, iwxr, info2;
    extern /* Subroutine */ int ab04md_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     dg01md_(char *, integer *, doublereal *, doublereal *, integer *,
	     ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *)
	    ;
    static integer iwmag, iwdme;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), sb10zp_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *);
    static integer iwvar, istop;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlamch_(char *, 
	    ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer iwbmat;
    extern /* Subroutine */ int dgelsy_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *);
    static integer clwmax, dlwmax, iwymag, iwdomo, istart;


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

/*     To fit a supplied frequency response data with a stable, minimum */
/*     phase SISO (single-input single-output) system represented by its */
/*     matrices A, B, C, D. It handles both discrete- and continuous-time */
/*     cases. */

/*     ARGUMENTS */

/*     Input/Output parameters */

/*     DISCFL  (input) INTEGER */
/*             Indicates the type of the system, as follows: */
/*             = 0: continuous-time system; */
/*             = 1: discrete-time system. */

/*     FLAG    (input) INTEGER */
/*             If FLAG = 0, then the system zeros and poles are not */
/*             constrained. */
/*             If FLAG = 1, then the system zeros and poles will have */
/*             negative real parts in the continuous-time case, or moduli */
/*             less than 1 in the discrete-time case. Consequently, FLAG */
/*             must be equal to 1 in mu-synthesis routines. */

/*     LENDAT  (input) INTEGER */
/*             The length of the vectors RFRDAT, IFRDAT and OMEGA. */
/*             LENDAT >= 2. */

/*     RFRDAT  (input) DOUBLE PRECISION array, dimension (LENDAT) */
/*             The real part of the frequency data to be fitted. */

/*     IFRDAT  (input) DOUBLE PRECISION array, dimension (LENDAT) */
/*             The imaginary part of the frequency data to be fitted. */

/*     OMEGA   (input) DOUBLE PRECISION array, dimension (LENDAT) */
/*             The frequencies corresponding to RFRDAT and IFRDAT. */
/*             These values must be nonnegative and monotonically */
/*             increasing. Additionally, for discrete-time systems */
/*             they must be between 0 and PI. */

/*     N       (input/output) INTEGER */
/*             On entry, the desired order of the system to be fitted. */
/*             N <= LENDAT-1. */
/*             On exit, the order of the obtained system. The value of N */
/*             could only be modified if N > 0 and FLAG = 1. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             matrix A. If FLAG = 1, then A is in an upper Hessenberg */
/*             form, and corresponds to a minimal realization. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (output) DOUBLE PRECISION array, dimension (N) */
/*             The computed vector B. */

/*     C       (output) DOUBLE PRECISION array, dimension (N) */
/*             The computed vector C. If FLAG = 1, the first N-1 elements */
/*             are zero (for the exit value of N). */

/*     D       (output) DOUBLE PRECISION array, dimension (1) */
/*             The computed scalar D. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used for determining the effective */
/*             rank of matrices. If the user sets TOL > 0, then the given */
/*             value of TOL is used as a lower bound for the reciprocal */
/*             condition number;  a (sub)matrix whose estimated condition */
/*             number is less than 1/TOL is considered to be of full */
/*             rank.  If the user sets TOL <= 0, then an implicitly */
/*             computed, default tolerance, defined by TOLDEF = SIZE*EPS, */
/*             is used instead, where SIZE is the product of the matrix */
/*             dimensions, and EPS is the machine precision (see LAPACK */
/*             Library routine DLAMCH). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension max(2,2*N+1) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK and DWORK(2) contains the optimal value of */
/*             LZWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK = max( 2, LW1, LW2, LW3, LW4 ), where */
/*             LW1 = 2*LENDAT + 4*HNPTS;  HNPTS = 2048; */
/*             LW2 =   LENDAT + 6*HNPTS; */
/*             MN  = min( 2*LENDAT, 2*N+1 ) */
/*             LW3 = 2*LENDAT*(2*N+1) + max( 2*LENDAT, 2*N+1 ) + */
/*                   max( MN + 6*N + 4, 2*MN + 1 ), if N > 0; */
/*             LW3 = 4*LENDAT + 5                 , if N = 0; */
/*             LW4 = max( N*N + 5*N, 6*N + 1 + min( 1,N ) ), if FLAG = 1; */
/*             LW4 = 0,                                      if FLAG = 0. */
/*             For optimum performance LDWORK should be larger. */

/*     ZWORK   COMPLEX*16 array, dimension (LZWORK) */

/*     LZWORK  INTEGER */
/*             The length of the array ZWORK. */
/*             LZWORK = LENDAT*(2*N+3), if N > 0; */
/*             LZWORK = LENDAT,         if N = 0. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the discrete --> continuous transformation cannot */
/*                   be made; */
/*             = 2:  if the system poles cannot be found; */
/*             = 3:  if the inverse system cannot be found, i.e., D is */
/*                   (close to) zero; */
/*             = 4:  if the system zeros cannot be found; */
/*             = 5:  if the state-space representation of the new */
/*                   transfer function T(s) cannot be found; */
/*             = 6:  if the continuous --> discrete transformation cannot */
/*                   be made. */

/*     METHOD */

/*     First, if the given frequency data are corresponding to a */
/*     continuous-time system, they are changed to a discrete-time */
/*     system using a bilinear transformation with a scaled alpha. */
/*     Then, the magnitude is obtained from the supplied data. */
/*     Then, the frequency data are linearly interpolated around */
/*     the unit-disc. */
/*     Then, Oppenheim and Schafer complex cepstrum method is applied */
/*     to get frequency data corresponding to a stable, minimum- */
/*     phase system. This is done in the following steps: */
/*     - Obtain LOG (magnitude) */
/*     - Obtain IFFT of the result (DG01MD SLICOT subroutine); */
/*     - halve the data at 0; */
/*     - Obtain FFT of the halved data (DG01MD SLICOT subroutine); */
/*     - Obtain EXP of the result. */
/*     Then, the new frequency data are interpolated back to the */
/*     original frequency. */
/*     Then, based on these newly obtained data, the system matrices */
/*     A, B, C, D are constructed; the very identification is */
/*     performed by Least Squares Method using DGELSY LAPACK subroutine. */
/*     If needed, a discrete-to-continuous time transformation is */
/*     applied on the system matrices by AB04MD SLICOT subroutine. */
/*     Finally, if requested, the poles and zeros of the system are */
/*     checked. If some of them have positive real parts in the */
/*     continuous-time case (or are not inside the unit disk in the */
/*     complex plane in the discrete-time case), they are exchanged with */
/*     their negatives (or reciprocals, respectively), to preserve the */
/*     frequency response, while getting a minimum phase and stable */
/*     system. This is done by SB10ZP SLICOT subroutine. */

/*     REFERENCES */

/*     [1] Oppenheim, A.V. and Schafer, R.W. */
/*         Discrete-Time Signal Processing. */
/*         Prentice-Hall Signal Processing Series, 1989. */

/*     [2] Balas, G., Doyle, J., Glover, K., Packard, A., and Smith, R. */
/*         Mu-analysis and Synthesis toolbox - User's Guide, */
/*         The Mathworks Inc., Natick, MA, USA, 1998. */

/*     CONTRIBUTORS */

/*     Asparuh Markovski, Technical University of Sofia, July 2003. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 2003. */
/*     A. Markovski, Technical University of Sofia, October 2003. */

/*     KEYWORDS */

/*     Bilinear transformation, frequency response, least-squares */
/*     approximation, stability. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */

/*     Test input parameters and workspace. */

#line 247 "SB10YD.f"
    /* Parameter adjustments */
#line 247 "SB10YD.f"
    --rfrdat;
#line 247 "SB10YD.f"
    --ifrdat;
#line 247 "SB10YD.f"
    --omega;
#line 247 "SB10YD.f"
    a_dim1 = *lda;
#line 247 "SB10YD.f"
    a_offset = 1 + a_dim1;
#line 247 "SB10YD.f"
    a -= a_offset;
#line 247 "SB10YD.f"
    --b;
#line 247 "SB10YD.f"
    --c__;
#line 247 "SB10YD.f"
    --d__;
#line 247 "SB10YD.f"
    --iwork;
#line 247 "SB10YD.f"
    --dwork;
#line 247 "SB10YD.f"
    --zwork;
#line 247 "SB10YD.f"

#line 247 "SB10YD.f"
    /* Function Body */
#line 247 "SB10YD.f"
    pi = atan(1.) * 4.;
#line 248 "SB10YD.f"
    pw = omega[1];
#line 249 "SB10YD.f"
    n1 = *n + 1;
#line 250 "SB10YD.f"
    n2 = *n + n1;

#line 252 "SB10YD.f"
    *info = 0;
#line 253 "SB10YD.f"
    if (*discfl != 0 && *discfl != 1) {
#line 254 "SB10YD.f"
	*info = -1;
#line 255 "SB10YD.f"
    } else if (*flag__ != 0 && *flag__ != 1) {
#line 256 "SB10YD.f"
	*info = -2;
#line 257 "SB10YD.f"
    } else if (*lendat < 2) {
#line 258 "SB10YD.f"
	*info = -3;
#line 259 "SB10YD.f"
    } else if (pw < 0.) {
#line 260 "SB10YD.f"
	*info = -6;
#line 261 "SB10YD.f"
    } else if (*n > *lendat - 1) {
#line 262 "SB10YD.f"
	*info = -7;
#line 263 "SB10YD.f"
    } else if (*lda < max(1,*n)) {
#line 264 "SB10YD.f"
	*info = -9;
#line 265 "SB10YD.f"
    } else {

#line 267 "SB10YD.f"
	i__1 = *lendat;
#line 267 "SB10YD.f"
	for (k = 2; k <= i__1; ++k) {
#line 268 "SB10YD.f"
	    if (omega[k] < pw) {
#line 268 "SB10YD.f"
		*info = -6;
#line 268 "SB10YD.f"
	    }
#line 270 "SB10YD.f"
	    pw = omega[k];
#line 271 "SB10YD.f"
/* L10: */
#line 271 "SB10YD.f"
	}

#line 273 "SB10YD.f"
	if (*discfl == 1 && omega[*lendat] > pi) {
#line 273 "SB10YD.f"
	    *info = -6;
#line 273 "SB10YD.f"
	}
#line 275 "SB10YD.f"
    }

#line 277 "SB10YD.f"
    if (*info == 0) {

/*        Workspace. */

#line 281 "SB10YD.f"
	lw1 = (*lendat << 1) + 8192;
#line 282 "SB10YD.f"
	lw2 = *lendat + 12288;
/* Computing MIN */
#line 283 "SB10YD.f"
	i__1 = *lendat << 1;
#line 283 "SB10YD.f"
	mn = min(i__1,n2);

#line 285 "SB10YD.f"
	if (*n > 0) {
/* Computing MAX */
#line 286 "SB10YD.f"
	    i__1 = *lendat << 1;
/* Computing MAX */
#line 286 "SB10YD.f"
	    i__2 = mn + *n * 6 + 4, i__3 = (mn << 1) + 1;
#line 286 "SB10YD.f"
	    lw3 = (*lendat << 1) * n2 + max(i__1,n2) + max(i__2,i__3);
#line 288 "SB10YD.f"
	} else {
#line 289 "SB10YD.f"
	    lw3 = (*lendat << 2) + 5;
#line 290 "SB10YD.f"
	}

#line 292 "SB10YD.f"
	if (*flag__ == 0) {
#line 293 "SB10YD.f"
	    lw4 = 0;
#line 294 "SB10YD.f"
	} else {
/* Computing MAX */
#line 295 "SB10YD.f"
	    i__1 = *n * *n + *n * 5, i__2 = *n * 6 + 1 + min(1,*n);
#line 295 "SB10YD.f"
	    lw4 = max(i__1,i__2);
#line 296 "SB10YD.f"
	}

/* Computing MAX */
#line 298 "SB10YD.f"
	i__1 = max(2,lw1), i__1 = max(i__1,lw2), i__1 = max(i__1,lw3);
#line 298 "SB10YD.f"
	dlwmax = max(i__1,lw4);

#line 300 "SB10YD.f"
	if (*n > 0) {
#line 301 "SB10YD.f"
	    clwmax = *lendat * (n2 + 2);
#line 302 "SB10YD.f"
	} else {
#line 303 "SB10YD.f"
	    clwmax = *lendat;
#line 304 "SB10YD.f"
	}

#line 306 "SB10YD.f"
	if (*ldwork < dlwmax) {
#line 307 "SB10YD.f"
	    *info = -16;
#line 308 "SB10YD.f"
	} else if (*lzwork < clwmax) {
#line 309 "SB10YD.f"
	    *info = -18;
#line 310 "SB10YD.f"
	}
#line 311 "SB10YD.f"
    }

#line 313 "SB10YD.f"
    if (*info != 0) {

/*        Error return. */

#line 317 "SB10YD.f"
	i__1 = -(*info);
#line 317 "SB10YD.f"
	xerbla_("SB10YD", &i__1, (ftnlen)6);
#line 318 "SB10YD.f"
	return 0;
#line 319 "SB10YD.f"
    }

/*     Set tolerances. */

#line 323 "SB10YD.f"
    tolb = dlamch_("Epsilon", (ftnlen)7);
#line 324 "SB10YD.f"
    toll = *tol;
#line 325 "SB10YD.f"
    if (toll <= 0.) {
#line 325 "SB10YD.f"
	toll = (doublereal) (*lendat * *n) * 4. * tolb;
#line 325 "SB10YD.f"
    }

/*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     Workspace usage 1. */
/*     Workspace:  need  2*LENDAT + 4*HNPTS. */

#line 333 "SB10YD.f"
    iwdomo = 1;
#line 334 "SB10YD.f"
    iwdme = iwdomo + *lendat;
#line 335 "SB10YD.f"
    iwymag = iwdme + 4096;
#line 336 "SB10YD.f"
    iwmag = iwymag + 4096;

/*     Bilinear transformation. */

#line 340 "SB10YD.f"
    if (*discfl == 0) {
#line 341 "SB10YD.f"
	pw = sqrt(omega[1] * omega[*lendat] + sqrt(tolb));

#line 343 "SB10YD.f"
	i__1 = *lendat;
#line 343 "SB10YD.f"
	for (k = 1; k <= i__1; ++k) {
/* Computing 2nd power */
#line 344 "SB10YD.f"
	    d__1 = omega[k] / pw;
#line 344 "SB10YD.f"
	    dwork[iwdme + k - 1] = d__1 * d__1;
#line 345 "SB10YD.f"
	    dwork[iwdomo + k - 1] = acos((1. - dwork[iwdme + k - 1]) / (dwork[
		    iwdme + k - 1] + 1.));
#line 348 "SB10YD.f"
/* L20: */
#line 348 "SB10YD.f"
	}

#line 350 "SB10YD.f"
    } else {
#line 351 "SB10YD.f"
	dcopy_(lendat, &omega[1], &c__1, &dwork[iwdomo], &c__1);
#line 352 "SB10YD.f"
    }

/*     Linear interpolation. */

#line 356 "SB10YD.f"
    i__1 = *lendat;
#line 356 "SB10YD.f"
    for (k = 1; k <= i__1; ++k) {
#line 357 "SB10YD.f"
	dwork[iwmag + k - 1] = dlapy2_(&rfrdat[k], &ifrdat[k]);
#line 358 "SB10YD.f"
	dwork[iwmag + k - 1] = 1. / log(10.) * log(dwork[iwmag + k - 1]);
#line 359 "SB10YD.f"
/* L30: */
#line 359 "SB10YD.f"
    }

#line 361 "SB10YD.f"
    for (k = 1; k <= 2048; ++k) {
#line 362 "SB10YD.f"
	dwork[iwdme + k - 1] = (k - 1) * pi / 2048;
#line 363 "SB10YD.f"
	dwork[iwymag + k - 1] = 0.;

#line 365 "SB10YD.f"
	if (dwork[iwdme + k - 1] < dwork[iwdomo]) {
#line 366 "SB10YD.f"
	    dwork[iwymag + k - 1] = dwork[iwmag];
#line 367 "SB10YD.f"
	} else if (dwork[iwdme + k - 1] >= dwork[iwdomo + *lendat - 1]) {
#line 368 "SB10YD.f"
	    dwork[iwymag + k - 1] = dwork[iwmag + *lendat - 1];
#line 369 "SB10YD.f"
	}

#line 371 "SB10YD.f"
/* L40: */
#line 371 "SB10YD.f"
    }

#line 373 "SB10YD.f"
    i__1 = *lendat;
#line 373 "SB10YD.f"
    for (i__ = 2; i__ <= i__1; ++i__) {
#line 374 "SB10YD.f"
	p1 = dwork[iwdomo + i__ - 2] * 2048 / pi + 1.;

#line 376 "SB10YD.f"
	ip1 = (integer) p1;
#line 377 "SB10YD.f"
	if ((doublereal) ip1 != p1) {
#line 377 "SB10YD.f"
	    ++ip1;
#line 377 "SB10YD.f"
	}

#line 380 "SB10YD.f"
	p2 = dwork[iwdomo + i__ - 1] * 2048 / pi + 1.;

#line 382 "SB10YD.f"
	ip2 = (integer) p2;
#line 383 "SB10YD.f"
	if ((doublereal) ip2 != p2) {
#line 383 "SB10YD.f"
	    ++ip2;
#line 383 "SB10YD.f"
	}

#line 386 "SB10YD.f"
	i__2 = ip2 - 1;
#line 386 "SB10YD.f"
	for (p = ip1; p <= i__2; ++p) {
#line 387 "SB10YD.f"
	    rat = dwork[iwdme + p - 1] - dwork[iwdomo + i__ - 2];
#line 388 "SB10YD.f"
	    rat /= dwork[iwdomo + i__ - 1] - dwork[iwdomo + i__ - 2];
#line 389 "SB10YD.f"
	    dwork[iwymag + p - 1] = (1. - rat) * dwork[iwmag + i__ - 2] + rat 
		    * dwork[iwmag + i__ - 1];
#line 391 "SB10YD.f"
/* L50: */
#line 391 "SB10YD.f"
	}

#line 393 "SB10YD.f"
/* L60: */
#line 393 "SB10YD.f"
    }

#line 395 "SB10YD.f"
    for (k = 1; k <= 2048; ++k) {
#line 396 "SB10YD.f"
	dwork[iwymag + k - 1] = exp(log(10.) * dwork[iwymag + k - 1]);
#line 397 "SB10YD.f"
/* L70: */
#line 397 "SB10YD.f"
    }

/*     Duplicate data around disc. */

#line 401 "SB10YD.f"
    for (k = 1; k <= 2048; ++k) {
#line 402 "SB10YD.f"
	dwork[iwdme + 2048 + k - 1] = pi * 2. - dwork[iwdme + 2048 - k];
#line 403 "SB10YD.f"
	dwork[iwymag + 2048 + k - 1] = dwork[iwymag + 2048 - k];
#line 404 "SB10YD.f"
/* L80: */
#line 404 "SB10YD.f"
    }

/*     Complex cepstrum to get min phase: */
/*     LOG (Magnitude) */

#line 409 "SB10YD.f"
    for (k = 1; k <= 4096; ++k) {
#line 410 "SB10YD.f"
	dwork[iwymag + k - 1] = log(dwork[iwymag + k - 1]) * 2.;
#line 411 "SB10YD.f"
/* L90: */
#line 411 "SB10YD.f"
    }

/*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*     Workspace usage 2. */
/*     Workspace:  need  LENDAT + 6*HNPTS. */

#line 418 "SB10YD.f"
    iwxr = iwymag;
#line 419 "SB10YD.f"
    iwxi = iwmag;

#line 421 "SB10YD.f"
    for (k = 1; k <= 4096; ++k) {
#line 422 "SB10YD.f"
	dwork[iwxi + k - 1] = 0.;
#line 423 "SB10YD.f"
/* L100: */
#line 423 "SB10YD.f"
    }

/*     IFFT */

#line 427 "SB10YD.f"
    dg01md_("I", &c__4096, &dwork[iwxr], &dwork[iwxi], &info2, (ftnlen)1);

/*     Rescale, because DG01MD doesn't do it. */

#line 431 "SB10YD.f"
    dscal_(&c__2048, &c_b21, &dwork[iwxr], &c__1);
#line 432 "SB10YD.f"
    dscal_(&c__2048, &c_b21, &dwork[iwxi], &c__1);

/*     Halve the result at 0. */

#line 436 "SB10YD.f"
    dwork[iwxr] /= 2.;
#line 437 "SB10YD.f"
    dwork[iwxi] /= 2.;

/*     FFT */

#line 441 "SB10YD.f"
    dg01md_("D", &c__2048, &dwork[iwxr], &dwork[iwxi], &info2, (ftnlen)1);

/*     Get the EXP of the result. */

#line 445 "SB10YD.f"
    for (k = 1; k <= 1024; ++k) {
#line 446 "SB10YD.f"
	i__1 = k - 1;
#line 446 "SB10YD.f"
	d__1 = exp(dwork[iwxr + k - 1]);
#line 446 "SB10YD.f"
	d__2 = cos(dwork[iwxi + k - 1]);
#line 446 "SB10YD.f"
	d__3 = sin(dwork[iwxi + k - 1]);
#line 446 "SB10YD.f"
	z__2.r = d__2, z__2.i = d__3;
#line 446 "SB10YD.f"
	z__1.r = d__1 * z__2.r, z__1.i = d__1 * z__2.i;
#line 446 "SB10YD.f"
	xhat[i__1].r = z__1.r, xhat[i__1].i = z__1.i;
#line 448 "SB10YD.f"
	dwork[iwdme + k - 1] = dwork[iwdme + (k << 1) - 2];
#line 449 "SB10YD.f"
/* L110: */
#line 449 "SB10YD.f"
    }

/*     Interpolate back to original frequency data. */

#line 453 "SB10YD.f"
    istart = 1;
#line 454 "SB10YD.f"
    istop = *lendat;

#line 456 "SB10YD.f"
    i__1 = *lendat;
#line 456 "SB10YD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 457 "SB10YD.f"
	i__2 = i__;
#line 457 "SB10YD.f"
	zwork[i__2].r = 0., zwork[i__2].i = 0.;
#line 458 "SB10YD.f"
	if (dwork[iwdomo + i__ - 1] <= dwork[iwdme]) {
#line 459 "SB10YD.f"
	    i__2 = i__;
#line 459 "SB10YD.f"
	    zwork[i__2].r = xhat[0].r, zwork[i__2].i = xhat[0].i;
#line 460 "SB10YD.f"
	    istart = i__ + 1;
#line 461 "SB10YD.f"
	} else if (dwork[iwdomo + i__ - 1] >= dwork[iwdme + 1023]) {
#line 463 "SB10YD.f"
	    i__2 = i__;
#line 463 "SB10YD.f"
	    zwork[i__2].r = xhat[1023].r, zwork[i__2].i = xhat[1023].i;
#line 464 "SB10YD.f"
	    --istop;
#line 465 "SB10YD.f"
	}
#line 466 "SB10YD.f"
/* L120: */
#line 466 "SB10YD.f"
    }

#line 468 "SB10YD.f"
    i__1 = istop;
#line 468 "SB10YD.f"
    for (i__ = istart; i__ <= i__1; ++i__) {
#line 469 "SB10YD.f"
	ii = 1024;
#line 470 "SB10YD.f"
L130:
#line 471 "SB10YD.f"
	if (dwork[iwdme + ii - 1] >= dwork[iwdomo + i__ - 1]) {
#line 471 "SB10YD.f"
	    p = ii;
#line 471 "SB10YD.f"
	}
#line 473 "SB10YD.f"
	--ii;
#line 474 "SB10YD.f"
	if (ii > 0) {
#line 474 "SB10YD.f"
	    goto L130;
#line 474 "SB10YD.f"
	}
#line 476 "SB10YD.f"
	rat = (dwork[iwdomo + i__ - 1] - dwork[iwdme + p - 2]) / (dwork[iwdme 
		+ p - 1] - dwork[iwdme + p - 2]);
#line 478 "SB10YD.f"
	i__2 = i__;
#line 478 "SB10YD.f"
	i__3 = p - 1;
#line 478 "SB10YD.f"
	z__2.r = rat * xhat[i__3].r, z__2.i = rat * xhat[i__3].i;
#line 478 "SB10YD.f"
	d__1 = 1. - rat;
#line 478 "SB10YD.f"
	i__4 = p - 2;
#line 478 "SB10YD.f"
	z__3.r = d__1 * xhat[i__4].r, z__3.i = d__1 * xhat[i__4].i;
#line 478 "SB10YD.f"
	z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
#line 478 "SB10YD.f"
	zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 479 "SB10YD.f"
/* L140: */
#line 479 "SB10YD.f"
    }

/*     CASE N > 0. */
/*     This is the only allowed case in mu-synthesis subroutines. */

#line 484 "SB10YD.f"
    if (*n > 0) {

/*        Preparation for frequency identification. */

/*        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*        Complex workspace usage 1. */
/*        Complex workspace:  need  2*LENDAT + LENDAT*(N+1). */

#line 493 "SB10YD.f"
	iwa0 = *lendat + 1;
#line 494 "SB10YD.f"
	iwvar = iwa0 + *lendat * n1;

#line 496 "SB10YD.f"
	i__1 = *lendat;
#line 496 "SB10YD.f"
	for (k = 1; k <= i__1; ++k) {
#line 497 "SB10YD.f"
	    if (*discfl == 0) {
#line 498 "SB10YD.f"
		i__2 = iwvar + k - 1;
#line 498 "SB10YD.f"
		d__1 = cos(dwork[iwdomo + k - 1]);
#line 498 "SB10YD.f"
		d__2 = sin(dwork[iwdomo + k - 1]);
#line 498 "SB10YD.f"
		z__1.r = d__1, z__1.i = d__2;
#line 498 "SB10YD.f"
		zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 500 "SB10YD.f"
	    } else {
#line 501 "SB10YD.f"
		i__2 = iwvar + k - 1;
#line 501 "SB10YD.f"
		d__1 = cos(omega[k]);
#line 501 "SB10YD.f"
		d__2 = sin(omega[k]);
#line 501 "SB10YD.f"
		z__1.r = d__1, z__1.i = d__2;
#line 501 "SB10YD.f"
		zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 503 "SB10YD.f"
	    }
#line 504 "SB10YD.f"
/* L150: */
#line 504 "SB10YD.f"
	}

/*        Array for DGELSY. */

#line 508 "SB10YD.f"
	i__1 = n2;
#line 508 "SB10YD.f"
	for (k = 1; k <= i__1; ++k) {
#line 509 "SB10YD.f"
	    iwork[k] = 0;
#line 510 "SB10YD.f"
/* L160: */
#line 510 "SB10YD.f"
	}

/*        Constructing A0. */

#line 514 "SB10YD.f"
	i__1 = *lendat;
#line 514 "SB10YD.f"
	for (k = 1; k <= i__1; ++k) {
#line 515 "SB10YD.f"
	    i__2 = iwa0 + *n * *lendat + k - 1;
#line 515 "SB10YD.f"
	    zwork[i__2].r = 1., zwork[i__2].i = 0.;
#line 516 "SB10YD.f"
/* L170: */
#line 516 "SB10YD.f"
	}

#line 518 "SB10YD.f"
	i__1 = *n;
#line 518 "SB10YD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 519 "SB10YD.f"
	    i__2 = *lendat;
#line 519 "SB10YD.f"
	    for (k = 1; k <= i__2; ++k) {
#line 520 "SB10YD.f"
		i__3 = iwa0 + (*n - i__) * *lendat + k - 1;
#line 520 "SB10YD.f"
		i__4 = iwa0 + (n1 - i__) * *lendat + k - 1;
#line 520 "SB10YD.f"
		i__5 = iwvar + k - 1;
#line 520 "SB10YD.f"
		z__1.r = zwork[i__4].r * zwork[i__5].r - zwork[i__4].i * 
			zwork[i__5].i, z__1.i = zwork[i__4].r * zwork[i__5].i 
			+ zwork[i__4].i * zwork[i__5].r;
#line 520 "SB10YD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 522 "SB10YD.f"
/* L180: */
#line 522 "SB10YD.f"
	    }
#line 523 "SB10YD.f"
/* L190: */
#line 523 "SB10YD.f"
	}

/*        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*        Complex workspace usage 2. */
/*        Complex workspace:  need  2*LENDAT + LENDAT*(2*N+1). */

#line 530 "SB10YD.f"
	iwbp = iwvar;
#line 531 "SB10YD.f"
	iwab = iwbp + *lendat;

/*        Constructing BP. */

#line 535 "SB10YD.f"
	i__1 = *lendat;
#line 535 "SB10YD.f"
	for (k = 1; k <= i__1; ++k) {
#line 536 "SB10YD.f"
	    i__2 = iwbp + k - 1;
#line 536 "SB10YD.f"
	    i__3 = iwa0 + k - 1;
#line 536 "SB10YD.f"
	    i__4 = k;
#line 536 "SB10YD.f"
	    z__1.r = zwork[i__3].r * zwork[i__4].r - zwork[i__3].i * zwork[
		    i__4].i, z__1.i = zwork[i__3].r * zwork[i__4].i + zwork[
		    i__3].i * zwork[i__4].r;
#line 536 "SB10YD.f"
	    zwork[i__2].r = z__1.r, zwork[i__2].i = z__1.i;
#line 537 "SB10YD.f"
/* L200: */
#line 537 "SB10YD.f"
	}

/*        Constructing AB. */

#line 541 "SB10YD.f"
	i__1 = *n;
#line 541 "SB10YD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 542 "SB10YD.f"
	    i__2 = *lendat;
#line 542 "SB10YD.f"
	    for (k = 1; k <= i__2; ++k) {
#line 543 "SB10YD.f"
		i__3 = iwab + (i__ - 1) * *lendat + k - 1;
#line 543 "SB10YD.f"
		i__4 = k;
#line 543 "SB10YD.f"
		z__2.r = -zwork[i__4].r, z__2.i = -zwork[i__4].i;
#line 543 "SB10YD.f"
		i__5 = iwa0 + i__ * *lendat + k - 1;
#line 543 "SB10YD.f"
		z__1.r = z__2.r * zwork[i__5].r - z__2.i * zwork[i__5].i, 
			z__1.i = z__2.r * zwork[i__5].i + z__2.i * zwork[i__5]
			.r;
#line 543 "SB10YD.f"
		zwork[i__3].r = z__1.r, zwork[i__3].i = z__1.i;
#line 545 "SB10YD.f"
/* L210: */
#line 545 "SB10YD.f"
	    }
#line 546 "SB10YD.f"
/* L220: */
#line 546 "SB10YD.f"
	}

/*        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*        Workspace usage 3. */
/*        Workspace:  need  LW3 = 2*LENDAT*(2*N+1) + max(2*LENDAT,2*N+1). */

#line 553 "SB10YD.f"
	iwbx = (*lendat << 1) * n2 + 1;
/* Computing MAX */
#line 554 "SB10YD.f"
	i__1 = *lendat << 1;
#line 554 "SB10YD.f"
	iws = iwbx + max(i__1,n2);

/*        Constructing AX. */

#line 558 "SB10YD.f"
	i__1 = n1;
#line 558 "SB10YD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 559 "SB10YD.f"
	    i__2 = *lendat;
#line 559 "SB10YD.f"
	    for (k = 1; k <= i__2; ++k) {
#line 560 "SB10YD.f"
		i__3 = iwa0 + (i__ - 1) * *lendat + k - 1;
#line 560 "SB10YD.f"
		dwork[(i__ - 1 << 1) * *lendat + k] = zwork[i__3].r;
#line 562 "SB10YD.f"
		dwork[((i__ << 1) - 1) * *lendat + k] = d_imag(&zwork[iwa0 + (
			i__ - 1) * *lendat + k - 1]);
#line 564 "SB10YD.f"
/* L230: */
#line 564 "SB10YD.f"
	    }
#line 565 "SB10YD.f"
/* L240: */
#line 565 "SB10YD.f"
	}

#line 567 "SB10YD.f"
	i__1 = *n;
#line 567 "SB10YD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 568 "SB10YD.f"
	    i__2 = *lendat;
#line 568 "SB10YD.f"
	    for (k = 1; k <= i__2; ++k) {
#line 569 "SB10YD.f"
		i__3 = iwab + (i__ - 1) * *lendat + k - 1;
#line 569 "SB10YD.f"
		dwork[(n1 << 1) * *lendat + (i__ - 1 << 1) * *lendat + k] = 
			zwork[i__3].r;
#line 571 "SB10YD.f"
		dwork[(n1 << 1) * *lendat + ((i__ << 1) - 1) * *lendat + k] = 
			d_imag(&zwork[iwab + (i__ - 1) * *lendat + k - 1]);
#line 573 "SB10YD.f"
/* L250: */
#line 573 "SB10YD.f"
	    }
#line 574 "SB10YD.f"
/* L260: */
#line 574 "SB10YD.f"
	}

/*        Constructing BX. */

#line 578 "SB10YD.f"
	i__1 = *lendat;
#line 578 "SB10YD.f"
	for (k = 1; k <= i__1; ++k) {
#line 579 "SB10YD.f"
	    i__2 = iwbp + k - 1;
#line 579 "SB10YD.f"
	    dwork[iwbx + k - 1] = zwork[i__2].r;
#line 580 "SB10YD.f"
	    dwork[iwbx + *lendat + k - 1] = d_imag(&zwork[iwbp + k - 1]);
#line 581 "SB10YD.f"
/* L270: */
#line 581 "SB10YD.f"
	}

/*        Estimating X. */
/*        Workspace:  need    LW3 + max( MN+3*(2*N+1)+1, 2*MN+1 ), */
/*                            where MN = min( 2*LENDAT, 2*N+1 ); */
/*                            prefer  larger. */

#line 588 "SB10YD.f"
	i__1 = *lendat << 1;
#line 588 "SB10YD.f"
	i__2 = *lendat << 1;
/* Computing MAX */
#line 588 "SB10YD.f"
	i__4 = *lendat << 1;
#line 588 "SB10YD.f"
	i__3 = max(i__4,n2);
#line 588 "SB10YD.f"
	i__5 = *ldwork - iws + 1;
#line 588 "SB10YD.f"
	dgelsy_(&i__1, &n2, &c__1, &dwork[1], &i__2, &dwork[iwbx], &i__3, &
		iwork[1], &toll, &rank, &dwork[iws], &i__5, &info2);
/* Computing MAX */
#line 591 "SB10YD.f"
	i__1 = dlwmax, i__2 = (integer) (dwork[iws] + iws - 1);
#line 591 "SB10YD.f"
	dlwmax = max(i__1,i__2);

/*        Constructing A matrix. */

#line 595 "SB10YD.f"
	i__1 = *n;
#line 595 "SB10YD.f"
	for (k = 1; k <= i__1; ++k) {
#line 596 "SB10YD.f"
	    a[k + a_dim1] = -dwork[iwbx + n1 + k - 1];
#line 597 "SB10YD.f"
/* L280: */
#line 597 "SB10YD.f"
	}

#line 599 "SB10YD.f"
	if (*n > 1) {
#line 599 "SB10YD.f"
	    i__1 = *n - 1;
#line 599 "SB10YD.f"
	    dlaset_("Full", n, &i__1, &c_b48, &c_b49, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)4);
#line 599 "SB10YD.f"
	}

/*        Constructing B matrix. */

#line 604 "SB10YD.f"
	i__1 = *n;
#line 604 "SB10YD.f"
	for (k = 1; k <= i__1; ++k) {
#line 605 "SB10YD.f"
	    b[k] = dwork[iwbx + n1 + k - 1] * dwork[iwbx] - dwork[iwbx + k];
#line 606 "SB10YD.f"
/* L290: */
#line 606 "SB10YD.f"
	}

/*        Constructing C matrix. */

#line 610 "SB10YD.f"
	c__[1] = -1.;

#line 612 "SB10YD.f"
	i__1 = *n;
#line 612 "SB10YD.f"
	for (k = 2; k <= i__1; ++k) {
#line 613 "SB10YD.f"
	    c__[k] = 0.;
#line 614 "SB10YD.f"
/* L300: */
#line 614 "SB10YD.f"
	}

/*        Constructing D matrix. */

#line 618 "SB10YD.f"
	d__[1] = dwork[iwbx];

/*        Transform to continuous-time case, if needed. */
/*        Workspace:  need    max(1,N); */
/*                            prefer  larger. */

#line 624 "SB10YD.f"
	if (*discfl == 0) {
#line 625 "SB10YD.f"
	    ab04md_("D", n, &c__1, &c__1, &c_b49, &pw, &a[a_offset], lda, &b[
		    1], lda, &c__[1], &c__1, &d__[1], &c__1, &iwork[1], &
		    dwork[1], ldwork, &info2, (ftnlen)1);
#line 627 "SB10YD.f"
	    if (info2 != 0) {
#line 628 "SB10YD.f"
		*info = 1;
#line 629 "SB10YD.f"
		return 0;
#line 630 "SB10YD.f"
	    }
/* Computing MAX */
#line 631 "SB10YD.f"
	    i__1 = dlwmax, i__2 = (integer) dwork[1];
#line 631 "SB10YD.f"
	    dlwmax = max(i__1,i__2);
#line 632 "SB10YD.f"
	}

/*        Make all the real parts of the poles and the zeros negative. */

#line 636 "SB10YD.f"
	if (*flag__ == 1) {

/*           Workspace:  need    max(N*N + 5*N, 6*N + 1 + min(1,N)); */
/*                               prefer  larger. */
#line 640 "SB10YD.f"
	    sb10zp_(discfl, n, &a[a_offset], lda, &b[1], &c__[1], &d__[1], &
		    iwork[1], &dwork[1], ldwork, info);
#line 642 "SB10YD.f"
	    if (*info != 0) {
#line 642 "SB10YD.f"
		return 0;
#line 642 "SB10YD.f"
	    }
/* Computing MAX */
#line 644 "SB10YD.f"
	    i__1 = dlwmax, i__2 = (integer) dwork[1];
#line 644 "SB10YD.f"
	    dlwmax = max(i__1,i__2);
#line 645 "SB10YD.f"
	}

#line 647 "SB10YD.f"
    } else {

/*        CASE N = 0. */

/*        @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/*        Workspace usage 4. */
/*        Workspace:  need  4*LENDAT. */

#line 656 "SB10YD.f"
	iwbmat = (*lendat << 1) + 1;
#line 657 "SB10YD.f"
	iws = iwbmat + (*lendat << 1);

/*        Constructing AMAT and BMAT. */

#line 661 "SB10YD.f"
	i__1 = *lendat;
#line 661 "SB10YD.f"
	for (k = 1; k <= i__1; ++k) {
#line 662 "SB10YD.f"
	    dwork[k] = 1.;
#line 663 "SB10YD.f"
	    dwork[k + *lendat] = 0.;
#line 664 "SB10YD.f"
	    i__2 = k;
#line 664 "SB10YD.f"
	    dwork[iwbmat + k - 1] = zwork[i__2].r;
#line 665 "SB10YD.f"
	    dwork[iwbmat + *lendat + k - 1] = d_imag(&zwork[k]);
#line 666 "SB10YD.f"
/* L310: */
#line 666 "SB10YD.f"
	}

/*        Estimating D matrix. */
/*        Workspace:  need    4*LENDAT + 5; */
/*                            prefer  larger. */

#line 672 "SB10YD.f"
	iwork[1] = 0;
#line 673 "SB10YD.f"
	i__1 = *lendat << 1;
#line 673 "SB10YD.f"
	i__2 = *lendat << 1;
#line 673 "SB10YD.f"
	i__3 = *lendat << 1;
#line 673 "SB10YD.f"
	i__4 = *ldwork - iws + 1;
#line 673 "SB10YD.f"
	dgelsy_(&i__1, &c__1, &c__1, &dwork[1], &i__2, &dwork[iwbmat], &i__3, 
		&iwork[1], &toll, &rank, &dwork[iws], &i__4, &info2);
/* Computing MAX */
#line 676 "SB10YD.f"
	i__1 = dlwmax, i__2 = (integer) (dwork[iws] + iws - 1);
#line 676 "SB10YD.f"
	dlwmax = max(i__1,i__2);

#line 678 "SB10YD.f"
	d__[1] = dwork[iwbmat];

#line 680 "SB10YD.f"
    }

/*     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#line 684 "SB10YD.f"
    dwork[1] = (doublereal) dlwmax;
#line 685 "SB10YD.f"
    dwork[2] = (doublereal) clwmax;
#line 686 "SB10YD.f"
    return 0;

/* *** Last line of SB10YD *** */
} /* sb10yd_ */

