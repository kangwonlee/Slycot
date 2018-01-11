#line 1 "MB04YW.f"
/* MB04YW.f -- translated by f2c (version 20100827).
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

#line 1 "MB04YW.f"
/* Table of constant values */

static doublereal c_b16 = 1.;

/* Subroutine */ int mb04yw_(logical *qrit, logical *updatu, logical *updatv, 
	integer *m, integer *n, integer *l, integer *k, doublereal *shift, 
	doublereal *d__, doublereal *e, doublereal *u, integer *ldu, 
	doublereal *v, integer *ldv, doublereal *dwork)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal f, g, h__;
    static integer i__;
    static doublereal r__, cs, sn;
    static integer nm1, nm12, nm13, ncv;
    static doublereal cosl, sinl, cosr, sinr;
    static integer irot;
    static doublereal oldcs;
    extern /* Subroutine */ int dlasr_(char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static doublereal oldsn;
    extern /* Subroutine */ int dlartg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);


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

/*     To perform either one QR or QL iteration step onto the unreduced */
/*     bidiagonal submatrix Jk: */

/*              |D(l) E(l)    0  ...    0   | */
/*              | 0   D(l+1) E(l+1)     .   | */
/*         Jk = | .                     .   | */
/*              | .                     .   | */
/*              | .                   E(k-1)| */
/*              | 0   ...        ...   D(k) | */

/*     with k <= p and l >= 1, p = MIN(M,N), of the bidiagonal matrix J: */

/*              |D(1) E(1)  0    ...   0   | */
/*              | 0   D(2) E(2)        .   | */
/*          J = | .                    .   |. */
/*              | .                    .   | */
/*              | .                  E(p-1)| */
/*              | 0   ...        ...  D(p) | */

/*     Hereby, Jk is transformed to  S' Jk T with S and T products of */
/*     Givens rotations. These Givens rotations S (respectively, T) are */
/*     postmultiplied into U (respectively, V), if UPDATU (respectively, */
/*     UPDATV) is .TRUE.. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     QRIT    LOGICAL */
/*             Indicates whether a QR or QL iteration step is to be */
/*             taken (from larger end diagonal element towards smaller), */
/*             as follows: */
/*             = .TRUE. :  QR iteration step (chase bulge from top to */
/*                         bottom); */
/*             = .FALSE.:  QL iteration step (chase bulge from bottom to */
/*                         top). */

/*     UPDATU  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix U the left-hand Givens rotations S, as follows: */
/*             = .FALSE.:  Do not form U; */
/*             = .TRUE. :  The given matrix U is updated (postmultiplied) */
/*                         by the left-hand Givens rotations S. */

/*     UPDATV  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix V the right-hand Givens rotations S, as follows: */
/*             = .FALSE.:  Do not form V; */
/*             = .TRUE. :  The given matrix V is updated (postmultiplied) */
/*                         by the right-hand Givens rotations T. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix U.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of rows of the matrix V.  N >= 0. */

/*     L       (input) INTEGER */
/*             The index of the first diagonal entry of the considered */
/*             unreduced bidiagonal submatrix Jk of J. */

/*     K       (input) INTEGER */
/*             The index of the last diagonal entry of the considered */
/*             unreduced bidiagonal submatrix Jk of J. */

/*     SHIFT   (input) DOUBLE PRECISION */
/*             Value of the shift used in the QR or QL iteration step. */

/*     D       (input/output) DOUBLE PRECISION array, dimension (p) */
/*             where p = MIN(M,N) */
/*             On entry, D must contain the diagonal entries of the */
/*             bidiagonal matrix J. */
/*             On exit, D contains the diagonal entries of the */
/*             transformed bidiagonal matrix S' J T. */

/*     E       (input/output) DOUBLE PRECISION array, dimension (p-1) */
/*             On entry, E must contain the superdiagonal entries of J. */
/*             On exit, E contains the superdiagonal entries of the */
/*             transformed matrix S' J T. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,p) */
/*             On entry, if UPDATU = .TRUE., U must contain the M-by-p */
/*             left transformation matrix. */
/*             On exit, if UPDATU = .TRUE., the Givens rotations S on the */
/*             left have been postmultiplied into U, i.e., U * S is */
/*             returned. */
/*             U is not referenced if UPDATU = .FALSE.. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= max(1,M) if UPDATU = .TRUE.; */
/*             LDU >= 1        if UPDATU = .FALSE.. */

/*     V       (input/output) DOUBLE PRECISION array, dimension (LDV,p) */
/*             On entry, if UPDATV = .TRUE., V must contain the N-by-p */
/*             right transformation matrix. */
/*             On exit, if UPDATV = .TRUE., the Givens rotations T on the */
/*             right have been postmultiplied into V, i.e., V * T is */
/*             returned. */
/*             V is not referenced if UPDATV = .FALSE.. */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= max(1,N) if UPDATV = .TRUE.; */
/*             LDV >= 1        if UPDATV = .FALSE.. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(1,LDWORK)) */
/*             LDWORK >= 4*MIN(M,N)-4, if UPDATU = UPDATV = .TRUE.; */
/*             LDWORK >= 2*MIN(M,N)-2, if */
/*                             UPDATU = .TRUE. and UPDATV = .FALSE. or */
/*                             UPDATV = .TRUE. and UPDATU = .FALSE.; */
/*             LDWORK >= 1, if UPDATU = UPDATV = .FALSE.. */

/*     METHOD */

/*     QR iterations diagonalize the bidiagonal matrix by zeroing the */
/*     super-diagonal elements of Jk from bottom to top. */
/*     QL iterations diagonalize the bidiagonal matrix by zeroing the */
/*     super-diagonal elements of Jk from top to bottom. */
/*     The routine overwrites Jk with the bidiagonal matrix S' Jk T, */
/*     where S and T are products of Givens rotations. */
/*     T is essentially the orthogonal matrix that would be obtained by */
/*     applying one implicit symmetric shift QR (QL) step onto the matrix */
/*     Jk'Jk. This step factors the matrix (Jk'Jk - shift*I) into a */
/*     product of an orthogonal matrix T and a upper (lower) triangular */
/*     matrix. See [1,Sec.8.2-8.3] and [2] for more details. */

/*     REFERENCES */

/*     [1] Golub, G.H. and Van Loan, C.F. */
/*         Matrix Computations. */
/*         The Johns Hopkins University Press, Baltimore, Maryland, 1983. */

/*     [2] Bowdler, H., Martin, R.S. and Wilkinson, J.H. */
/*         The QR and QL algorithms for symmetric matrices. */
/*         Numer. Math., 11, pp. 293-306, 1968. */

/*     [3] Demmel, J. and Kahan, W. */
/*         Computing small singular values of bidiagonal matrices with */
/*         guaranteed high relative accuracy. */
/*         SIAM J. Sci. Statist. Comput., 11, pp. 873-912, 1990. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routines MB04QY and MB04QZ by S. Van */
/*     Huffel, Katholieke University Leuven, Belgium. */
/*     This subroutine is based on the QR/QL step implemented in LAPACK */
/*     routine DBDSQR. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Bidiagonal matrix, orthogonal transformation, singular values. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     For speed, no tests of the input scalar arguments are done. */

/*     Quick return if possible. */

#line 220 "MB04YW.f"
    /* Parameter adjustments */
#line 220 "MB04YW.f"
    --d__;
#line 220 "MB04YW.f"
    --e;
#line 220 "MB04YW.f"
    u_dim1 = *ldu;
#line 220 "MB04YW.f"
    u_offset = 1 + u_dim1;
#line 220 "MB04YW.f"
    u -= u_offset;
#line 220 "MB04YW.f"
    v_dim1 = *ldv;
#line 220 "MB04YW.f"
    v_offset = 1 + v_dim1;
#line 220 "MB04YW.f"
    v -= v_offset;
#line 220 "MB04YW.f"
    --dwork;
#line 220 "MB04YW.f"

#line 220 "MB04YW.f"
    /* Function Body */
#line 220 "MB04YW.f"
    ncv = min(*m,*n);
#line 221 "MB04YW.f"
    if (ncv <= 1 || *l == *k) {
#line 221 "MB04YW.f"
	return 0;
#line 221 "MB04YW.f"
    }

#line 224 "MB04YW.f"
    nm1 = ncv - 1;
#line 225 "MB04YW.f"
    nm12 = nm1 + nm1;
#line 226 "MB04YW.f"
    nm13 = nm12 + nm1;
#line 227 "MB04YW.f"
    if (! (*updatv)) {
#line 228 "MB04YW.f"
	nm12 = 0;
#line 229 "MB04YW.f"
	nm13 = nm1;
#line 230 "MB04YW.f"
    }

/*     If SHIFT = 0, do simplified QR iteration. */

#line 234 "MB04YW.f"
    if (*shift == 0.) {
#line 235 "MB04YW.f"
	if (*qrit) {

/*           Chase bulge from top to bottom. */
/*           Save cosines and sines for later U and/or V updates, */
/*           if needed. */

#line 241 "MB04YW.f"
	    cs = 1.;
#line 242 "MB04YW.f"
	    oldcs = 1.;
#line 243 "MB04YW.f"
	    d__1 = d__[*l] * cs;
#line 243 "MB04YW.f"
	    dlartg_(&d__1, &e[*l], &cs, &sn, &r__);
#line 244 "MB04YW.f"
	    d__1 = oldcs * r__;
#line 244 "MB04YW.f"
	    d__2 = d__[*l + 1] * sn;
#line 244 "MB04YW.f"
	    dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[*l]);
#line 245 "MB04YW.f"
	    if (*updatv) {
#line 246 "MB04YW.f"
		dwork[1] = cs;
#line 247 "MB04YW.f"
		dwork[nm1 + 1] = sn;
#line 248 "MB04YW.f"
	    }
#line 249 "MB04YW.f"
	    if (*updatu) {
#line 250 "MB04YW.f"
		dwork[nm12 + 1] = oldcs;
#line 251 "MB04YW.f"
		dwork[nm13 + 1] = oldsn;
#line 252 "MB04YW.f"
	    }
#line 253 "MB04YW.f"
	    irot = 1;

#line 255 "MB04YW.f"
	    i__1 = *k - 1;
#line 255 "MB04YW.f"
	    for (i__ = *l + 1; i__ <= i__1; ++i__) {
#line 256 "MB04YW.f"
		d__1 = d__[i__] * cs;
#line 256 "MB04YW.f"
		dlartg_(&d__1, &e[i__], &cs, &sn, &r__);
#line 257 "MB04YW.f"
		e[i__ - 1] = oldsn * r__;
#line 258 "MB04YW.f"
		d__1 = oldcs * r__;
#line 258 "MB04YW.f"
		d__2 = d__[i__ + 1] * sn;
#line 258 "MB04YW.f"
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 259 "MB04YW.f"
		++irot;
#line 260 "MB04YW.f"
		if (*updatv) {
#line 261 "MB04YW.f"
		    dwork[irot] = cs;
#line 262 "MB04YW.f"
		    dwork[irot + nm1] = sn;
#line 263 "MB04YW.f"
		}
#line 264 "MB04YW.f"
		if (*updatu) {
#line 265 "MB04YW.f"
		    dwork[irot + nm12] = oldcs;
#line 266 "MB04YW.f"
		    dwork[irot + nm13] = oldsn;
#line 267 "MB04YW.f"
		}
#line 268 "MB04YW.f"
/* L110: */
#line 268 "MB04YW.f"
	    }

#line 270 "MB04YW.f"
	    h__ = d__[*k] * cs;
#line 271 "MB04YW.f"
	    d__[*k] = h__ * oldcs;
#line 272 "MB04YW.f"
	    e[*k - 1] = h__ * oldsn;

/*           Update U and/or V. */

#line 276 "MB04YW.f"
	    if (*updatv) {
#line 276 "MB04YW.f"
		i__1 = *k - *l + 1;
#line 276 "MB04YW.f"
		dlasr_("R", "V", "F", n, &i__1, &dwork[1], &dwork[ncv], &v[*l 
			* v_dim1 + 1], ldv, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 276 "MB04YW.f"
	    }
#line 279 "MB04YW.f"
	    if (*updatu) {
#line 279 "MB04YW.f"
		i__1 = *k - *l + 1;
#line 279 "MB04YW.f"
		dlasr_("R", "V", "F", m, &i__1, &dwork[nm12 + 1], &dwork[nm13 
			+ 1], &u[*l * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 279 "MB04YW.f"
	    }

#line 283 "MB04YW.f"
	} else {

/*           Chase bulge from bottom to top. */
/*           Save cosines and sines for later U and/or V updates, */
/*           if needed. */

#line 289 "MB04YW.f"
	    cs = 1.;
#line 290 "MB04YW.f"
	    oldcs = 1.;
#line 291 "MB04YW.f"
	    d__1 = d__[*k] * cs;
#line 291 "MB04YW.f"
	    dlartg_(&d__1, &e[*k - 1], &cs, &sn, &r__);
#line 292 "MB04YW.f"
	    d__1 = oldcs * r__;
#line 292 "MB04YW.f"
	    d__2 = d__[*k - 1] * sn;
#line 292 "MB04YW.f"
	    dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[*k]);
#line 293 "MB04YW.f"
	    if (*updatv) {
#line 294 "MB04YW.f"
		dwork[*k - *l] = oldcs;
#line 295 "MB04YW.f"
		dwork[*k - *l + nm1] = -oldsn;
#line 296 "MB04YW.f"
	    }
#line 297 "MB04YW.f"
	    if (*updatu) {
#line 298 "MB04YW.f"
		dwork[*k - *l + nm12] = cs;
#line 299 "MB04YW.f"
		dwork[*k - *l + nm13] = -sn;
#line 300 "MB04YW.f"
	    }
#line 301 "MB04YW.f"
	    irot = *k - *l;

#line 303 "MB04YW.f"
	    i__1 = *l + 1;
#line 303 "MB04YW.f"
	    for (i__ = *k - 1; i__ >= i__1; --i__) {
#line 304 "MB04YW.f"
		d__1 = d__[i__] * cs;
#line 304 "MB04YW.f"
		dlartg_(&d__1, &e[i__ - 1], &cs, &sn, &r__);
#line 305 "MB04YW.f"
		e[i__] = oldsn * r__;
#line 306 "MB04YW.f"
		d__1 = oldcs * r__;
#line 306 "MB04YW.f"
		d__2 = d__[i__ - 1] * sn;
#line 306 "MB04YW.f"
		dlartg_(&d__1, &d__2, &oldcs, &oldsn, &d__[i__]);
#line 307 "MB04YW.f"
		--irot;
#line 308 "MB04YW.f"
		if (*updatv) {
#line 309 "MB04YW.f"
		    dwork[irot] = oldcs;
#line 310 "MB04YW.f"
		    dwork[irot + nm1] = -oldsn;
#line 311 "MB04YW.f"
		}
#line 312 "MB04YW.f"
		if (*updatu) {
#line 313 "MB04YW.f"
		    dwork[irot + nm12] = cs;
#line 314 "MB04YW.f"
		    dwork[irot + nm13] = -sn;
#line 315 "MB04YW.f"
		}
#line 316 "MB04YW.f"
/* L120: */
#line 316 "MB04YW.f"
	    }

#line 318 "MB04YW.f"
	    h__ = d__[*l] * cs;
#line 319 "MB04YW.f"
	    d__[*l] = h__ * oldcs;
#line 320 "MB04YW.f"
	    e[*l] = h__ * oldsn;

/*           Update U and/or V. */

#line 324 "MB04YW.f"
	    if (*updatv) {
#line 324 "MB04YW.f"
		i__1 = *k - *l + 1;
#line 324 "MB04YW.f"
		dlasr_("R", "V", "B", n, &i__1, &dwork[1], &dwork[ncv], &v[*l 
			* v_dim1 + 1], ldv, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 324 "MB04YW.f"
	    }
#line 327 "MB04YW.f"
	    if (*updatu) {
#line 327 "MB04YW.f"
		i__1 = *k - *l + 1;
#line 327 "MB04YW.f"
		dlasr_("R", "V", "B", m, &i__1, &dwork[nm12 + 1], &dwork[nm13 
			+ 1], &u[*l * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 327 "MB04YW.f"
	    }
#line 330 "MB04YW.f"
	}
#line 331 "MB04YW.f"
    } else {

/*        Use nonzero shift. */

#line 335 "MB04YW.f"
	if (*qrit) {

/*           Chase bulge from top to bottom. */
/*           Save cosines and sines for later U and/or V updates, */
/*           if needed. */

#line 341 "MB04YW.f"
	    f = ((d__1 = d__[*l], abs(d__1)) - *shift) * (d_sign(&c_b16, &d__[
		    *l]) + *shift / d__[*l]);
#line 343 "MB04YW.f"
	    g = e[*l];
#line 344 "MB04YW.f"
	    dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 345 "MB04YW.f"
	    f = cosr * d__[*l] + sinr * e[*l];
#line 346 "MB04YW.f"
	    e[*l] = cosr * e[*l] - sinr * d__[*l];
#line 347 "MB04YW.f"
	    g = sinr * d__[*l + 1];
#line 348 "MB04YW.f"
	    d__[*l + 1] = cosr * d__[*l + 1];
#line 349 "MB04YW.f"
	    dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 350 "MB04YW.f"
	    d__[*l] = r__;
#line 351 "MB04YW.f"
	    f = cosl * e[*l] + sinl * d__[*l + 1];
#line 352 "MB04YW.f"
	    d__[*l + 1] = cosl * d__[*l + 1] - sinl * e[*l];
#line 353 "MB04YW.f"
	    g = sinl * e[*l + 1];
#line 354 "MB04YW.f"
	    e[*l + 1] = cosl * e[*l + 1];
#line 355 "MB04YW.f"
	    if (*updatv) {
#line 356 "MB04YW.f"
		dwork[1] = cosr;
#line 357 "MB04YW.f"
		dwork[nm1 + 1] = sinr;
#line 358 "MB04YW.f"
	    }
#line 359 "MB04YW.f"
	    if (*updatu) {
#line 360 "MB04YW.f"
		dwork[nm12 + 1] = cosl;
#line 361 "MB04YW.f"
		dwork[nm13 + 1] = sinl;
#line 362 "MB04YW.f"
	    }
#line 363 "MB04YW.f"
	    irot = 1;

#line 365 "MB04YW.f"
	    i__1 = *k - 2;
#line 365 "MB04YW.f"
	    for (i__ = *l + 1; i__ <= i__1; ++i__) {
#line 366 "MB04YW.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 367 "MB04YW.f"
		e[i__ - 1] = r__;
#line 368 "MB04YW.f"
		f = cosr * d__[i__] + sinr * e[i__];
#line 369 "MB04YW.f"
		e[i__] = cosr * e[i__] - sinr * d__[i__];
#line 370 "MB04YW.f"
		g = sinr * d__[i__ + 1];
#line 371 "MB04YW.f"
		d__[i__ + 1] = cosr * d__[i__ + 1];
#line 372 "MB04YW.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 373 "MB04YW.f"
		d__[i__] = r__;
#line 374 "MB04YW.f"
		f = cosl * e[i__] + sinl * d__[i__ + 1];
#line 375 "MB04YW.f"
		d__[i__ + 1] = cosl * d__[i__ + 1] - sinl * e[i__];
#line 376 "MB04YW.f"
		g = sinl * e[i__ + 1];
#line 377 "MB04YW.f"
		e[i__ + 1] = cosl * e[i__ + 1];
#line 378 "MB04YW.f"
		++irot;
#line 379 "MB04YW.f"
		if (*updatv) {
#line 380 "MB04YW.f"
		    dwork[irot] = cosr;
#line 381 "MB04YW.f"
		    dwork[irot + nm1] = sinr;
#line 382 "MB04YW.f"
		}
#line 383 "MB04YW.f"
		if (*updatu) {
#line 384 "MB04YW.f"
		    dwork[irot + nm12] = cosl;
#line 385 "MB04YW.f"
		    dwork[irot + nm13] = sinl;
#line 386 "MB04YW.f"
		}
#line 387 "MB04YW.f"
/* L130: */
#line 387 "MB04YW.f"
	    }

#line 389 "MB04YW.f"
	    if (*l < *k - 1) {
#line 390 "MB04YW.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 391 "MB04YW.f"
		e[*k - 2] = r__;
#line 392 "MB04YW.f"
		f = cosr * d__[*k - 1] + sinr * e[*k - 1];
#line 393 "MB04YW.f"
		e[*k - 1] = cosr * e[*k - 1] - sinr * d__[*k - 1];
#line 394 "MB04YW.f"
		g = sinr * d__[*k];
#line 395 "MB04YW.f"
		d__[*k] = cosr * d__[*k];
#line 396 "MB04YW.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 397 "MB04YW.f"
		d__[*k - 1] = r__;
#line 398 "MB04YW.f"
		f = cosl * e[*k - 1] + sinl * d__[*k];
#line 399 "MB04YW.f"
		d__[*k] = cosl * d__[*k] - sinl * e[*k - 1];
#line 400 "MB04YW.f"
		++irot;
#line 401 "MB04YW.f"
		if (*updatv) {
#line 402 "MB04YW.f"
		    dwork[irot] = cosr;
#line 403 "MB04YW.f"
		    dwork[irot + nm1] = sinr;
#line 404 "MB04YW.f"
		}
#line 405 "MB04YW.f"
		if (*updatu) {
#line 406 "MB04YW.f"
		    dwork[irot + nm12] = cosl;
#line 407 "MB04YW.f"
		    dwork[irot + nm13] = sinl;
#line 408 "MB04YW.f"
		}
#line 409 "MB04YW.f"
	    }
#line 410 "MB04YW.f"
	    e[*k - 1] = f;

/*           Update U and/or V. */

#line 414 "MB04YW.f"
	    if (*updatv) {
#line 414 "MB04YW.f"
		i__1 = *k - *l + 1;
#line 414 "MB04YW.f"
		dlasr_("R", "V", "F", n, &i__1, &dwork[1], &dwork[ncv], &v[*l 
			* v_dim1 + 1], ldv, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 414 "MB04YW.f"
	    }
#line 417 "MB04YW.f"
	    if (*updatu) {
#line 417 "MB04YW.f"
		i__1 = *k - *l + 1;
#line 417 "MB04YW.f"
		dlasr_("R", "V", "F", m, &i__1, &dwork[nm12 + 1], &dwork[nm13 
			+ 1], &u[*l * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 417 "MB04YW.f"
	    }

#line 421 "MB04YW.f"
	} else {

/*           Chase bulge from bottom to top. */
/*           Save cosines and sines for later U and/or V updates, */
/*           if needed. */

#line 427 "MB04YW.f"
	    f = ((d__1 = d__[*k], abs(d__1)) - *shift) * (d_sign(&c_b16, &d__[
		    *k]) + *shift / d__[*k]);
#line 429 "MB04YW.f"
	    g = e[*k - 1];
#line 430 "MB04YW.f"
	    if (*l < *k - 1) {
#line 431 "MB04YW.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 432 "MB04YW.f"
		f = cosr * d__[*k] + sinr * e[*k - 1];
#line 433 "MB04YW.f"
		e[*k - 1] = cosr * e[*k - 1] - sinr * d__[*k];
#line 434 "MB04YW.f"
		g = sinr * d__[*k - 1];
#line 435 "MB04YW.f"
		d__[*k - 1] = cosr * d__[*k - 1];
#line 436 "MB04YW.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 437 "MB04YW.f"
		d__[*k] = r__;
#line 438 "MB04YW.f"
		f = cosl * e[*k - 1] + sinl * d__[*k - 1];
#line 439 "MB04YW.f"
		d__[*k - 1] = cosl * d__[*k - 1] - sinl * e[*k - 1];
#line 440 "MB04YW.f"
		g = sinl * e[*k - 2];
#line 441 "MB04YW.f"
		e[*k - 2] = cosl * e[*k - 2];
#line 442 "MB04YW.f"
		if (*updatv) {
#line 443 "MB04YW.f"
		    dwork[*k - *l] = cosl;
#line 444 "MB04YW.f"
		    dwork[*k - *l + nm1] = -sinl;
#line 445 "MB04YW.f"
		}
#line 446 "MB04YW.f"
		if (*updatu) {
#line 447 "MB04YW.f"
		    dwork[*k - *l + nm12] = cosr;
#line 448 "MB04YW.f"
		    dwork[*k - *l + nm13] = -sinr;
#line 449 "MB04YW.f"
		}
#line 450 "MB04YW.f"
		irot = *k - *l;
#line 451 "MB04YW.f"
	    } else {
#line 452 "MB04YW.f"
		irot = *k - *l + 1;
#line 453 "MB04YW.f"
	    }

#line 455 "MB04YW.f"
	    i__1 = *l + 2;
#line 455 "MB04YW.f"
	    for (i__ = *k - 1; i__ >= i__1; --i__) {
#line 456 "MB04YW.f"
		dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 457 "MB04YW.f"
		e[i__] = r__;
#line 458 "MB04YW.f"
		f = cosr * d__[i__] + sinr * e[i__ - 1];
#line 459 "MB04YW.f"
		e[i__ - 1] = cosr * e[i__ - 1] - sinr * d__[i__];
#line 460 "MB04YW.f"
		g = sinr * d__[i__ - 1];
#line 461 "MB04YW.f"
		d__[i__ - 1] = cosr * d__[i__ - 1];
#line 462 "MB04YW.f"
		dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 463 "MB04YW.f"
		d__[i__] = r__;
#line 464 "MB04YW.f"
		f = cosl * e[i__ - 1] + sinl * d__[i__ - 1];
#line 465 "MB04YW.f"
		d__[i__ - 1] = cosl * d__[i__ - 1] - sinl * e[i__ - 1];
#line 466 "MB04YW.f"
		g = sinl * e[i__ - 2];
#line 467 "MB04YW.f"
		e[i__ - 2] = cosl * e[i__ - 2];
#line 468 "MB04YW.f"
		--irot;
#line 469 "MB04YW.f"
		if (*updatv) {
#line 470 "MB04YW.f"
		    dwork[irot] = cosl;
#line 471 "MB04YW.f"
		    dwork[irot + nm1] = -sinl;
#line 472 "MB04YW.f"
		}
#line 473 "MB04YW.f"
		if (*updatu) {
#line 474 "MB04YW.f"
		    dwork[irot + nm12] = cosr;
#line 475 "MB04YW.f"
		    dwork[irot + nm13] = -sinr;
#line 476 "MB04YW.f"
		}
#line 477 "MB04YW.f"
/* L140: */
#line 477 "MB04YW.f"
	    }

#line 479 "MB04YW.f"
	    dlartg_(&f, &g, &cosr, &sinr, &r__);
#line 480 "MB04YW.f"
	    e[*l + 1] = r__;
#line 481 "MB04YW.f"
	    f = cosr * d__[*l + 1] + sinr * e[*l];
#line 482 "MB04YW.f"
	    e[*l] = cosr * e[*l] - sinr * d__[*l + 1];
#line 483 "MB04YW.f"
	    g = sinr * d__[*l];
#line 484 "MB04YW.f"
	    d__[*l] = cosr * d__[*l];
#line 485 "MB04YW.f"
	    dlartg_(&f, &g, &cosl, &sinl, &r__);
#line 486 "MB04YW.f"
	    d__[*l + 1] = r__;
#line 487 "MB04YW.f"
	    f = cosl * e[*l] + sinl * d__[*l];
#line 488 "MB04YW.f"
	    d__[*l] = cosl * d__[*l] - sinl * e[*l];
#line 489 "MB04YW.f"
	    --irot;
#line 490 "MB04YW.f"
	    if (*updatv) {
#line 491 "MB04YW.f"
		dwork[irot] = cosl;
#line 492 "MB04YW.f"
		dwork[irot + nm1] = -sinl;
#line 493 "MB04YW.f"
	    }
#line 494 "MB04YW.f"
	    if (*updatu) {
#line 495 "MB04YW.f"
		dwork[irot + nm12] = cosr;
#line 496 "MB04YW.f"
		dwork[irot + nm13] = -sinr;
#line 497 "MB04YW.f"
	    }
#line 498 "MB04YW.f"
	    e[*l] = f;

/*           Update U and/or V if desired. */

#line 502 "MB04YW.f"
	    if (*updatv) {
#line 502 "MB04YW.f"
		i__1 = *k - *l + 1;
#line 502 "MB04YW.f"
		dlasr_("R", "V", "B", n, &i__1, &dwork[1], &dwork[ncv], &v[*l 
			* v_dim1 + 1], ldv, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 502 "MB04YW.f"
	    }
#line 505 "MB04YW.f"
	    if (*updatu) {
#line 505 "MB04YW.f"
		i__1 = *k - *l + 1;
#line 505 "MB04YW.f"
		dlasr_("R", "V", "B", m, &i__1, &dwork[nm12 + 1], &dwork[nm13 
			+ 1], &u[*l * u_dim1 + 1], ldu, (ftnlen)1, (ftnlen)1, 
			(ftnlen)1);
#line 505 "MB04YW.f"
	    }
#line 508 "MB04YW.f"
	}
#line 509 "MB04YW.f"
    }

#line 511 "MB04YW.f"
    return 0;
/* *** Last line of MB04YW *** */
} /* mb04yw_ */

