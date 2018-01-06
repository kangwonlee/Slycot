#line 1 "SB04PX.f"
/* SB04PX.f -- translated by f2c (version 20100827).
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

#line 1 "SB04PX.f"
/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;

/* Subroutine */ int sb04px_(logical *ltranl, logical *ltranr, integer *isgn, 
	integer *n1, integer *n2, doublereal *tl, integer *ldtl, doublereal *
	tr, integer *ldtr, doublereal *b, integer *ldb, doublereal *scale, 
	doublereal *x, integer *ldx, doublereal *xnorm, integer *info)
{
    /* Initialized data */

    static integer locu12[4] = { 3,4,1,2 };
    static integer locl21[4] = { 2,1,4,3 };
    static integer locu22[4] = { 4,3,2,1 };
    static logical xswpiv[4] = { FALSE_,FALSE_,TRUE_,TRUE_ };
    static logical bswpiv[4] = { FALSE_,TRUE_,FALSE_,TRUE_ };

    /* System generated locals */
    integer b_dim1, b_offset, tl_dim1, tl_offset, tr_dim1, tr_offset, x_dim1, 
	    x_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;

    /* Local variables */
    static integer i__, j, k;
    static doublereal x2[2], l21, u11, u12;
    static integer ip, jp;
    static doublereal u22, t16[16]	/* was [4][4] */, gam, bet, eps, sgn, 
	    tmp[4], tau1, btmp[4], smin;
    static integer ipiv;
    static doublereal temp;
    static integer jpiv[4];
    static doublereal xmax;
    static integer ipsv, jpsv;
    static logical bswap;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical xswap;
    extern doublereal dlamch_(char *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal smlnum;


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

/*     To solve for the N1-by-N2 matrix X, 1 <= N1,N2 <= 2, in */

/*            op(TL)*X*op(TR) + ISGN*X = SCALE*B, */

/*     where TL is N1-by-N1, TR is N2-by-N2, B is N1-by-N2, and ISGN = 1 */
/*     or -1.  op(T) = T or T', where T' denotes the transpose of T. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LTRANL  LOGICAL */
/*             Specifies the form of op(TL) to be used, as follows: */
/*             = .FALSE.:  op(TL) = TL, */
/*             = .TRUE. :  op(TL) = TL'. */

/*     LTRANR  LOGICAL */
/*             Specifies the form of op(TR) to be used, as follows: */
/*             = .FALSE.:  op(TR) = TR, */
/*             = .TRUE. :  op(TR) = TR'. */

/*     ISGN    INTEGER */
/*             Specifies the sign of the equation as described before. */
/*             ISGN may only be 1 or -1. */

/*     Input/Output Parameters */

/*     N1      (input) INTEGER */
/*             The order of matrix TL.  N1 may only be 0, 1 or 2. */

/*     N2      (input) INTEGER */
/*             The order of matrix TR.  N2 may only be 0, 1 or 2. */

/*     TL      (input) DOUBLE PRECISION array, dimension (LDTL,N1) */
/*             The leading N1-by-N1 part of this array must contain the */
/*             matrix TL. */

/*     LDTL    INTEGER */
/*             The leading dimension of array TL.  LDTL >= MAX(1,N1). */

/*     TR      (input) DOUBLE PRECISION array, dimension (LDTR,N2) */
/*             The leading N2-by-N2 part of this array must contain the */
/*             matrix TR. */

/*     LDTR    INTEGER */
/*             The leading dimension of array TR.  LDTR >= MAX(1,N2). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,N2) */
/*             The leading N1-by-N2 part of this array must contain the */
/*             right-hand side of the equation. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N1). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor. SCALE is chosen less than or equal to 1 */
/*             to prevent the solution overflowing. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,N2) */
/*             The leading N1-by-N2 part of this array contains the */
/*             solution of the equation. */
/*             Note that X may be identified with B in the calling */
/*             statement. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= MAX(1,N1). */

/*     XNORM   (output) DOUBLE PRECISION */
/*             The infinity-norm of the solution. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if TL and -ISGN*TR have almost reciprocal */
/*                   eigenvalues, so TL or TR is perturbed to get a */
/*                   nonsingular equation. */

/*             NOTE: In the interests of speed, this routine does not */
/*                   check the inputs for errors. */

/*     METHOD */

/*     The equivalent linear algebraic system of equations is formed and */
/*     solved using Gaussian elimination with complete pivoting. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is stable and reliable, since Gaussian elimination */
/*     with complete pivoting is used. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 2000. */
/*     This is a modification and slightly more efficient version of */
/*     SLICOT Library routine SB03MU. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Discrete-time system, Sylvester equation, matrix algebra. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Data statements .. */
#line 177 "SB04PX.f"
    /* Parameter adjustments */
#line 177 "SB04PX.f"
    tl_dim1 = *ldtl;
#line 177 "SB04PX.f"
    tl_offset = 1 + tl_dim1;
#line 177 "SB04PX.f"
    tl -= tl_offset;
#line 177 "SB04PX.f"
    tr_dim1 = *ldtr;
#line 177 "SB04PX.f"
    tr_offset = 1 + tr_dim1;
#line 177 "SB04PX.f"
    tr -= tr_offset;
#line 177 "SB04PX.f"
    b_dim1 = *ldb;
#line 177 "SB04PX.f"
    b_offset = 1 + b_dim1;
#line 177 "SB04PX.f"
    b -= b_offset;
#line 177 "SB04PX.f"
    x_dim1 = *ldx;
#line 177 "SB04PX.f"
    x_offset = 1 + x_dim1;
#line 177 "SB04PX.f"
    x -= x_offset;
#line 177 "SB04PX.f"

#line 177 "SB04PX.f"
    /* Function Body */
/*     .. */
/*     .. Executable Statements .. */

/*     Do not check the input parameters for errors. */

#line 186 "SB04PX.f"
    *info = 0;
#line 187 "SB04PX.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 191 "SB04PX.f"
    if (*n1 == 0 || *n2 == 0) {
#line 192 "SB04PX.f"
	*xnorm = 0.;
#line 193 "SB04PX.f"
	return 0;
#line 194 "SB04PX.f"
    }

/*     Set constants to control overflow. */

#line 198 "SB04PX.f"
    eps = dlamch_("P", (ftnlen)1);
#line 199 "SB04PX.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 200 "SB04PX.f"
    sgn = (doublereal) (*isgn);

#line 202 "SB04PX.f"
    k = *n1 + *n1 + *n2 - 2;
#line 203 "SB04PX.f"
    switch (k) {
#line 203 "SB04PX.f"
	case 1:  goto L10;
#line 203 "SB04PX.f"
	case 2:  goto L20;
#line 203 "SB04PX.f"
	case 3:  goto L30;
#line 203 "SB04PX.f"
	case 4:  goto L50;
#line 203 "SB04PX.f"
    }

/*     1-by-1: TL11*X*TR11 + ISGN*X = B11. */

#line 207 "SB04PX.f"
L10:
#line 208 "SB04PX.f"
    tau1 = tl[tl_dim1 + 1] * tr[tr_dim1 + 1] + sgn;
#line 209 "SB04PX.f"
    bet = abs(tau1);
#line 210 "SB04PX.f"
    if (bet <= smlnum) {
#line 211 "SB04PX.f"
	tau1 = smlnum;
#line 212 "SB04PX.f"
	bet = smlnum;
#line 213 "SB04PX.f"
	*info = 1;
#line 214 "SB04PX.f"
    }

#line 216 "SB04PX.f"
    gam = (d__1 = b[b_dim1 + 1], abs(d__1));
#line 217 "SB04PX.f"
    if (smlnum * gam > bet) {
#line 217 "SB04PX.f"
	*scale = 1. / gam;
#line 217 "SB04PX.f"
    }

#line 220 "SB04PX.f"
    x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / tau1;
#line 221 "SB04PX.f"
    *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1));
#line 222 "SB04PX.f"
    return 0;

/*     1-by-2: */
/*     TL11*[X11 X12]*op[TR11 TR12] + ISGN*[X11 X12] = [B11 B12]. */
/*                      [TR21 TR22] */

#line 228 "SB04PX.f"
L20:

/* Computing MAX */
/* Computing MAX */
#line 230 "SB04PX.f"
    d__7 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__8 = (d__2 = tr[(tr_dim1 << 
	    1) + 1], abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = tr[
	    tr_dim1 + 2], abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = 
	    tr[(tr_dim1 << 1) + 2], abs(d__4));
#line 230 "SB04PX.f"
    d__6 = max(d__7,d__8) * (d__5 = tl[tl_dim1 + 1], abs(d__5)) * eps;
#line 230 "SB04PX.f"
    smin = max(d__6,smlnum);
#line 234 "SB04PX.f"
    tmp[0] = tl[tl_dim1 + 1] * tr[tr_dim1 + 1] + sgn;
#line 235 "SB04PX.f"
    tmp[3] = tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 2] + sgn;
#line 236 "SB04PX.f"
    if (*ltranr) {
#line 237 "SB04PX.f"
	tmp[1] = tl[tl_dim1 + 1] * tr[tr_dim1 + 2];
#line 238 "SB04PX.f"
	tmp[2] = tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 1];
#line 239 "SB04PX.f"
    } else {
#line 240 "SB04PX.f"
	tmp[1] = tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 1];
#line 241 "SB04PX.f"
	tmp[2] = tl[tl_dim1 + 1] * tr[tr_dim1 + 2];
#line 242 "SB04PX.f"
    }
#line 243 "SB04PX.f"
    btmp[0] = b[b_dim1 + 1];
#line 244 "SB04PX.f"
    btmp[1] = b[(b_dim1 << 1) + 1];
#line 245 "SB04PX.f"
    goto L40;

/*     2-by-1: */
/*     op[TL11 TL12]*[X11]*TR11 + ISGN*[X11] = [B11]. */
/*       [TL21 TL22] [X21]             [X21]   [B21] */

#line 251 "SB04PX.f"
L30:
/* Computing MAX */
/* Computing MAX */
#line 252 "SB04PX.f"
    d__7 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__8 = (d__2 = tl[(tl_dim1 << 
	    1) + 1], abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = tl[
	    tl_dim1 + 2], abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = 
	    tl[(tl_dim1 << 1) + 2], abs(d__4));
#line 252 "SB04PX.f"
    d__6 = max(d__7,d__8) * (d__5 = tr[tr_dim1 + 1], abs(d__5)) * eps;
#line 252 "SB04PX.f"
    smin = max(d__6,smlnum);
#line 256 "SB04PX.f"
    tmp[0] = tl[tl_dim1 + 1] * tr[tr_dim1 + 1] + sgn;
#line 257 "SB04PX.f"
    tmp[3] = tl[(tl_dim1 << 1) + 2] * tr[tr_dim1 + 1] + sgn;
#line 258 "SB04PX.f"
    if (*ltranl) {
#line 259 "SB04PX.f"
	tmp[1] = tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 1];
#line 260 "SB04PX.f"
	tmp[2] = tl[tl_dim1 + 2] * tr[tr_dim1 + 1];
#line 261 "SB04PX.f"
    } else {
#line 262 "SB04PX.f"
	tmp[1] = tl[tl_dim1 + 2] * tr[tr_dim1 + 1];
#line 263 "SB04PX.f"
	tmp[2] = tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 1];
#line 264 "SB04PX.f"
    }
#line 265 "SB04PX.f"
    btmp[0] = b[b_dim1 + 1];
#line 266 "SB04PX.f"
    btmp[1] = b[b_dim1 + 2];
#line 267 "SB04PX.f"
L40:

/*     Solve 2-by-2 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

#line 272 "SB04PX.f"
    ipiv = idamax_(&c__4, tmp, &c__1);
#line 273 "SB04PX.f"
    u11 = tmp[ipiv - 1];
#line 274 "SB04PX.f"
    if (abs(u11) <= smin) {
#line 275 "SB04PX.f"
	*info = 1;
#line 276 "SB04PX.f"
	u11 = smin;
#line 277 "SB04PX.f"
    }
#line 278 "SB04PX.f"
    u12 = tmp[locu12[ipiv - 1] - 1];
#line 279 "SB04PX.f"
    l21 = tmp[locl21[ipiv - 1] - 1] / u11;
#line 280 "SB04PX.f"
    u22 = tmp[locu22[ipiv - 1] - 1] - u12 * l21;
#line 281 "SB04PX.f"
    xswap = xswpiv[ipiv - 1];
#line 282 "SB04PX.f"
    bswap = bswpiv[ipiv - 1];
#line 283 "SB04PX.f"
    if (abs(u22) <= smin) {
#line 284 "SB04PX.f"
	*info = 1;
#line 285 "SB04PX.f"
	u22 = smin;
#line 286 "SB04PX.f"
    }
#line 287 "SB04PX.f"
    if (bswap) {
#line 288 "SB04PX.f"
	temp = btmp[1];
#line 289 "SB04PX.f"
	btmp[1] = btmp[0] - l21 * temp;
#line 290 "SB04PX.f"
	btmp[0] = temp;
#line 291 "SB04PX.f"
    } else {
#line 292 "SB04PX.f"
	btmp[1] -= l21 * btmp[0];
#line 293 "SB04PX.f"
    }
#line 294 "SB04PX.f"
    if (smlnum * 2. * abs(btmp[1]) > abs(u22) || smlnum * 2. * abs(btmp[0]) > 
	    abs(u11)) {
/* Computing MAX */
#line 296 "SB04PX.f"
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]);
#line 296 "SB04PX.f"
	*scale = .5 / max(d__1,d__2);
#line 297 "SB04PX.f"
	btmp[0] *= *scale;
#line 298 "SB04PX.f"
	btmp[1] *= *scale;
#line 299 "SB04PX.f"
    }
#line 300 "SB04PX.f"
    x2[1] = btmp[1] / u22;
#line 301 "SB04PX.f"
    x2[0] = btmp[0] / u11 - u12 / u11 * x2[1];
#line 302 "SB04PX.f"
    if (xswap) {
#line 303 "SB04PX.f"
	temp = x2[1];
#line 304 "SB04PX.f"
	x2[1] = x2[0];
#line 305 "SB04PX.f"
	x2[0] = temp;
#line 306 "SB04PX.f"
    }
#line 307 "SB04PX.f"
    x[x_dim1 + 1] = x2[0];
#line 308 "SB04PX.f"
    if (*n1 == 1) {
#line 309 "SB04PX.f"
	x[(x_dim1 << 1) + 1] = x2[1];
#line 310 "SB04PX.f"
	*xnorm = abs(x2[0]) + abs(x2[1]);
#line 311 "SB04PX.f"
    } else {
#line 312 "SB04PX.f"
	x[x_dim1 + 2] = x2[1];
/* Computing MAX */
#line 313 "SB04PX.f"
	d__1 = abs(x2[0]), d__2 = abs(x2[1]);
#line 313 "SB04PX.f"
	*xnorm = max(d__1,d__2);
#line 314 "SB04PX.f"
    }
#line 315 "SB04PX.f"
    return 0;

/*     2-by-2: */
/*     op[TL11 TL12]*[X11 X12]*op[TR11 TR12] + ISGN*[X11 X12] = [B11 B12] */
/*       [TL21 TL22] [X21 X22]   [TR21 TR22]        [X21 X22]   [B21 B22] */

/*     Solve equivalent 4-by-4 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

#line 324 "SB04PX.f"
L50:
/* Computing MAX */
#line 325 "SB04PX.f"
    d__5 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__6 = (d__2 = tr[(tr_dim1 << 
	    1) + 1], abs(d__2)), d__5 = max(d__5,d__6), d__6 = (d__3 = tr[
	    tr_dim1 + 2], abs(d__3)), d__5 = max(d__5,d__6), d__6 = (d__4 = 
	    tr[(tr_dim1 << 1) + 2], abs(d__4));
#line 325 "SB04PX.f"
    smin = max(d__5,d__6);
/* Computing MAX */
#line 327 "SB04PX.f"
    d__5 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__6 = (d__2 = tl[(tl_dim1 << 
	    1) + 1], abs(d__2)), d__5 = max(d__5,d__6), d__6 = (d__3 = tl[
	    tl_dim1 + 2], abs(d__3)), d__5 = max(d__5,d__6), d__6 = (d__4 = 
	    tl[(tl_dim1 << 1) + 2], abs(d__4));
#line 327 "SB04PX.f"
    smin = max(d__5,d__6) * smin;
/* Computing MAX */
#line 329 "SB04PX.f"
    d__1 = eps * smin;
#line 329 "SB04PX.f"
    smin = max(d__1,smlnum);
#line 330 "SB04PX.f"
    t16[0] = tl[tl_dim1 + 1] * tr[tr_dim1 + 1] + sgn;
#line 331 "SB04PX.f"
    t16[5] = tl[(tl_dim1 << 1) + 2] * tr[tr_dim1 + 1] + sgn;
#line 332 "SB04PX.f"
    t16[10] = tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 2] + sgn;
#line 333 "SB04PX.f"
    t16[15] = tl[(tl_dim1 << 1) + 2] * tr[(tr_dim1 << 1) + 2] + sgn;
#line 334 "SB04PX.f"
    if (*ltranl) {
#line 335 "SB04PX.f"
	t16[4] = tl[tl_dim1 + 2] * tr[tr_dim1 + 1];
#line 336 "SB04PX.f"
	t16[1] = tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 1];
#line 337 "SB04PX.f"
	t16[14] = tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 2];
#line 338 "SB04PX.f"
	t16[11] = tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 2];
#line 339 "SB04PX.f"
    } else {
#line 340 "SB04PX.f"
	t16[4] = tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 1];
#line 341 "SB04PX.f"
	t16[1] = tl[tl_dim1 + 2] * tr[tr_dim1 + 1];
#line 342 "SB04PX.f"
	t16[14] = tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 2];
#line 343 "SB04PX.f"
	t16[11] = tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 2];
#line 344 "SB04PX.f"
    }
#line 345 "SB04PX.f"
    if (*ltranr) {
#line 346 "SB04PX.f"
	t16[8] = tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 1];
#line 347 "SB04PX.f"
	t16[13] = tl[(tl_dim1 << 1) + 2] * tr[(tr_dim1 << 1) + 1];
#line 348 "SB04PX.f"
	t16[2] = tl[tl_dim1 + 1] * tr[tr_dim1 + 2];
#line 349 "SB04PX.f"
	t16[7] = tl[(tl_dim1 << 1) + 2] * tr[tr_dim1 + 2];
#line 350 "SB04PX.f"
    } else {
#line 351 "SB04PX.f"
	t16[8] = tl[tl_dim1 + 1] * tr[tr_dim1 + 2];
#line 352 "SB04PX.f"
	t16[13] = tl[(tl_dim1 << 1) + 2] * tr[tr_dim1 + 2];
#line 353 "SB04PX.f"
	t16[2] = tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 1];
#line 354 "SB04PX.f"
	t16[7] = tl[(tl_dim1 << 1) + 2] * tr[(tr_dim1 << 1) + 1];
#line 355 "SB04PX.f"
    }
#line 356 "SB04PX.f"
    if (*ltranl && *ltranr) {
#line 357 "SB04PX.f"
	t16[12] = tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 1];
#line 358 "SB04PX.f"
	t16[9] = tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 1];
#line 359 "SB04PX.f"
	t16[6] = tl[tl_dim1 + 2] * tr[tr_dim1 + 2];
#line 360 "SB04PX.f"
	t16[3] = tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 2];
#line 361 "SB04PX.f"
    } else if (*ltranl && ! (*ltranr)) {
#line 362 "SB04PX.f"
	t16[12] = tl[tl_dim1 + 2] * tr[tr_dim1 + 2];
#line 363 "SB04PX.f"
	t16[9] = tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 2];
#line 364 "SB04PX.f"
	t16[6] = tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 1];
#line 365 "SB04PX.f"
	t16[3] = tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 1];
#line 366 "SB04PX.f"
    } else if (! (*ltranl) && *ltranr) {
#line 367 "SB04PX.f"
	t16[12] = tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 1];
#line 368 "SB04PX.f"
	t16[9] = tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 1];
#line 369 "SB04PX.f"
	t16[6] = tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 2];
#line 370 "SB04PX.f"
	t16[3] = tl[tl_dim1 + 2] * tr[tr_dim1 + 2];
#line 371 "SB04PX.f"
    } else {
#line 372 "SB04PX.f"
	t16[12] = tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 2];
#line 373 "SB04PX.f"
	t16[9] = tl[tl_dim1 + 2] * tr[tr_dim1 + 2];
#line 374 "SB04PX.f"
	t16[6] = tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 1];
#line 375 "SB04PX.f"
	t16[3] = tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 1];
#line 376 "SB04PX.f"
    }
#line 377 "SB04PX.f"
    btmp[0] = b[b_dim1 + 1];
#line 378 "SB04PX.f"
    btmp[1] = b[b_dim1 + 2];
#line 379 "SB04PX.f"
    btmp[2] = b[(b_dim1 << 1) + 1];
#line 380 "SB04PX.f"
    btmp[3] = b[(b_dim1 << 1) + 2];

/*     Perform elimination. */

#line 384 "SB04PX.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 385 "SB04PX.f"
	xmax = 0.;

#line 387 "SB04PX.f"
	for (ip = i__; ip <= 4; ++ip) {

#line 389 "SB04PX.f"
	    for (jp = i__; jp <= 4; ++jp) {
#line 390 "SB04PX.f"
		if ((d__1 = t16[ip + (jp << 2) - 5], abs(d__1)) >= xmax) {
#line 391 "SB04PX.f"
		    xmax = (d__1 = t16[ip + (jp << 2) - 5], abs(d__1));
#line 392 "SB04PX.f"
		    ipsv = ip;
#line 393 "SB04PX.f"
		    jpsv = jp;
#line 394 "SB04PX.f"
		}
#line 395 "SB04PX.f"
/* L60: */
#line 395 "SB04PX.f"
	    }

#line 397 "SB04PX.f"
/* L70: */
#line 397 "SB04PX.f"
	}

#line 399 "SB04PX.f"
	if (ipsv != i__) {
#line 400 "SB04PX.f"
	    dswap_(&c__4, &t16[ipsv - 1], &c__4, &t16[i__ - 1], &c__4);
#line 401 "SB04PX.f"
	    temp = btmp[i__ - 1];
#line 402 "SB04PX.f"
	    btmp[i__ - 1] = btmp[ipsv - 1];
#line 403 "SB04PX.f"
	    btmp[ipsv - 1] = temp;
#line 404 "SB04PX.f"
	}
#line 405 "SB04PX.f"
	if (jpsv != i__) {
#line 405 "SB04PX.f"
	    dswap_(&c__4, &t16[(jpsv << 2) - 4], &c__1, &t16[(i__ << 2) - 4], 
		    &c__1);
#line 405 "SB04PX.f"
	}
#line 407 "SB04PX.f"
	jpiv[i__ - 1] = jpsv;
#line 408 "SB04PX.f"
	if ((d__1 = t16[i__ + (i__ << 2) - 5], abs(d__1)) < smin) {
#line 409 "SB04PX.f"
	    *info = 1;
#line 410 "SB04PX.f"
	    t16[i__ + (i__ << 2) - 5] = smin;
#line 411 "SB04PX.f"
	}

#line 413 "SB04PX.f"
	for (j = i__ + 1; j <= 4; ++j) {
#line 414 "SB04PX.f"
	    t16[j + (i__ << 2) - 5] /= t16[i__ + (i__ << 2) - 5];
#line 415 "SB04PX.f"
	    btmp[j - 1] -= t16[j + (i__ << 2) - 5] * btmp[i__ - 1];

#line 417 "SB04PX.f"
	    for (k = i__ + 1; k <= 4; ++k) {
#line 418 "SB04PX.f"
		t16[j + (k << 2) - 5] -= t16[j + (i__ << 2) - 5] * t16[i__ + (
			k << 2) - 5];
#line 419 "SB04PX.f"
/* L80: */
#line 419 "SB04PX.f"
	    }

#line 421 "SB04PX.f"
/* L90: */
#line 421 "SB04PX.f"
	}

#line 423 "SB04PX.f"
/* L100: */
#line 423 "SB04PX.f"
    }

#line 425 "SB04PX.f"
    if (abs(t16[15]) < smin) {
#line 425 "SB04PX.f"
	t16[15] = smin;
#line 425 "SB04PX.f"
    }
#line 427 "SB04PX.f"
    if (smlnum * 8. * abs(btmp[0]) > abs(t16[0]) || smlnum * 8. * abs(btmp[1])
	     > abs(t16[5]) || smlnum * 8. * abs(btmp[2]) > abs(t16[10]) || 
	    smlnum * 8. * abs(btmp[3]) > abs(t16[15])) {
/* Computing MAX */
#line 431 "SB04PX.f"
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]), d__1 = max(d__1,d__2), d__2 
		= abs(btmp[2]), d__1 = max(d__1,d__2), d__2 = abs(btmp[3]);
#line 431 "SB04PX.f"
	*scale = .125 / max(d__1,d__2);
#line 434 "SB04PX.f"
	btmp[0] *= *scale;
#line 435 "SB04PX.f"
	btmp[1] *= *scale;
#line 436 "SB04PX.f"
	btmp[2] *= *scale;
#line 437 "SB04PX.f"
	btmp[3] *= *scale;
#line 438 "SB04PX.f"
    }

#line 440 "SB04PX.f"
    for (i__ = 1; i__ <= 4; ++i__) {
#line 441 "SB04PX.f"
	k = 5 - i__;
#line 442 "SB04PX.f"
	temp = 1. / t16[k + (k << 2) - 5];
#line 443 "SB04PX.f"
	tmp[k - 1] = btmp[k - 1] * temp;

#line 445 "SB04PX.f"
	for (j = k + 1; j <= 4; ++j) {
#line 446 "SB04PX.f"
	    tmp[k - 1] -= temp * t16[k + (j << 2) - 5] * tmp[j - 1];
#line 447 "SB04PX.f"
/* L110: */
#line 447 "SB04PX.f"
	}

#line 449 "SB04PX.f"
/* L120: */
#line 449 "SB04PX.f"
    }

#line 451 "SB04PX.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 452 "SB04PX.f"
	if (jpiv[4 - i__ - 1] != 4 - i__) {
#line 453 "SB04PX.f"
	    temp = tmp[4 - i__ - 1];
#line 454 "SB04PX.f"
	    tmp[4 - i__ - 1] = tmp[jpiv[4 - i__ - 1] - 1];
#line 455 "SB04PX.f"
	    tmp[jpiv[4 - i__ - 1] - 1] = temp;
#line 456 "SB04PX.f"
	}
#line 457 "SB04PX.f"
/* L130: */
#line 457 "SB04PX.f"
    }

#line 459 "SB04PX.f"
    x[x_dim1 + 1] = tmp[0];
#line 460 "SB04PX.f"
    x[x_dim1 + 2] = tmp[1];
#line 461 "SB04PX.f"
    x[(x_dim1 << 1) + 1] = tmp[2];
#line 462 "SB04PX.f"
    x[(x_dim1 << 1) + 2] = tmp[3];
/* Computing MAX */
#line 463 "SB04PX.f"
    d__1 = abs(tmp[0]) + abs(tmp[2]), d__2 = abs(tmp[1]) + abs(tmp[3]);
#line 463 "SB04PX.f"
    *xnorm = max(d__1,d__2);

#line 466 "SB04PX.f"
    return 0;
/* *** Last line of SB04PX *** */
} /* sb04px_ */

