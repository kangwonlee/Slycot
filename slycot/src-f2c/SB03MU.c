#line 1 "SB03MU.f"
/* SB03MU.f -- translated by f2c (version 20100827).
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

#line 1 "SB03MU.f"
/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;

/* Subroutine */ int sb03mu_(logical *ltranl, logical *ltranr, integer *isgn, 
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

/*            ISGN*op(TL)*X*op(TR) - X = SCALE*B, */

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

/*     TL      (input) DOUBLE PRECISION array, dimension (LDTL,2) */
/*             The leading N1-by-N1 part of this array must contain the */
/*             matrix TL. */

/*     LDTL    INTEGER */
/*             The leading dimension of array TL.  LDTL >= MAX(1,N1). */

/*     TR      (input) DOUBLE PRECISION array, dimension (LDTR,2) */
/*             The leading N2-by-N2 part of this array must contain the */
/*             matrix TR. */

/*     LDTR    INTEGER */
/*             The leading dimension of array TR.  LDTR >= MAX(1,N2). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,2) */
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
/*             = 1:  if TL and TR have almost reciprocal eigenvalues, so */
/*                   TL or TR is perturbed to get a nonsingular equation. */

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

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Based on DLASD2 by P. Petkov, Tech. University of Sofia, September */
/*     1993. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999. */

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
#line 176 "SB03MU.f"
    /* Parameter adjustments */
#line 176 "SB03MU.f"
    tl_dim1 = *ldtl;
#line 176 "SB03MU.f"
    tl_offset = 1 + tl_dim1;
#line 176 "SB03MU.f"
    tl -= tl_offset;
#line 176 "SB03MU.f"
    tr_dim1 = *ldtr;
#line 176 "SB03MU.f"
    tr_offset = 1 + tr_dim1;
#line 176 "SB03MU.f"
    tr -= tr_offset;
#line 176 "SB03MU.f"
    b_dim1 = *ldb;
#line 176 "SB03MU.f"
    b_offset = 1 + b_dim1;
#line 176 "SB03MU.f"
    b -= b_offset;
#line 176 "SB03MU.f"
    x_dim1 = *ldx;
#line 176 "SB03MU.f"
    x_offset = 1 + x_dim1;
#line 176 "SB03MU.f"
    x -= x_offset;
#line 176 "SB03MU.f"

#line 176 "SB03MU.f"
    /* Function Body */
/*     .. */
/*     .. Executable Statements .. */

/*     Do not check the input parameters for errors. */

#line 185 "SB03MU.f"
    *info = 0;
#line 186 "SB03MU.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 190 "SB03MU.f"
    if (*n1 == 0 || *n2 == 0) {
#line 191 "SB03MU.f"
	*xnorm = 0.;
#line 192 "SB03MU.f"
	return 0;
#line 193 "SB03MU.f"
    }

/*     Set constants to control overflow. */

#line 197 "SB03MU.f"
    eps = dlamch_("P", (ftnlen)1);
#line 198 "SB03MU.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 199 "SB03MU.f"
    sgn = (doublereal) (*isgn);

#line 201 "SB03MU.f"
    k = *n1 + *n1 + *n2 - 2;
#line 202 "SB03MU.f"
    switch (k) {
#line 202 "SB03MU.f"
	case 1:  goto L10;
#line 202 "SB03MU.f"
	case 2:  goto L20;
#line 202 "SB03MU.f"
	case 3:  goto L30;
#line 202 "SB03MU.f"
	case 4:  goto L50;
#line 202 "SB03MU.f"
    }

/*     1-by-1: SGN*TL11*X*TR11 - X = B11. */

#line 206 "SB03MU.f"
L10:
#line 207 "SB03MU.f"
    tau1 = sgn * tl[tl_dim1 + 1] * tr[tr_dim1 + 1] - 1.;
#line 208 "SB03MU.f"
    bet = abs(tau1);
#line 209 "SB03MU.f"
    if (bet <= smlnum) {
#line 210 "SB03MU.f"
	tau1 = smlnum;
#line 211 "SB03MU.f"
	bet = smlnum;
#line 212 "SB03MU.f"
	*info = 1;
#line 213 "SB03MU.f"
    }

#line 215 "SB03MU.f"
    gam = (d__1 = b[b_dim1 + 1], abs(d__1));
#line 216 "SB03MU.f"
    if (smlnum * gam > bet) {
#line 216 "SB03MU.f"
	*scale = 1. / gam;
#line 216 "SB03MU.f"
    }

#line 219 "SB03MU.f"
    x[x_dim1 + 1] = b[b_dim1 + 1] * *scale / tau1;
#line 220 "SB03MU.f"
    *xnorm = (d__1 = x[x_dim1 + 1], abs(d__1));
#line 221 "SB03MU.f"
    return 0;

/*     1-by-2: */
/*     ISGN*TL11*[X11 X12]*op[TR11 TR12]  = [B11 B12]. */
/*                           [TR21 TR22] */

#line 227 "SB03MU.f"
L20:

/* Computing MAX */
/* Computing MAX */
#line 229 "SB03MU.f"
    d__7 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__8 = (d__2 = tr[(tr_dim1 << 
	    1) + 1], abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = tr[
	    tr_dim1 + 2], abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = 
	    tr[(tr_dim1 << 1) + 2], abs(d__4));
#line 229 "SB03MU.f"
    d__6 = max(d__7,d__8) * (d__5 = tl[tl_dim1 + 1], abs(d__5)) * eps;
#line 229 "SB03MU.f"
    smin = max(d__6,smlnum);
#line 233 "SB03MU.f"
    tmp[0] = sgn * tl[tl_dim1 + 1] * tr[tr_dim1 + 1] - 1.;
#line 234 "SB03MU.f"
    tmp[3] = sgn * tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 2] - 1.;
#line 235 "SB03MU.f"
    if (*ltranr) {
#line 236 "SB03MU.f"
	tmp[1] = sgn * tl[tl_dim1 + 1] * tr[tr_dim1 + 2];
#line 237 "SB03MU.f"
	tmp[2] = sgn * tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 1];
#line 238 "SB03MU.f"
    } else {
#line 239 "SB03MU.f"
	tmp[1] = sgn * tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 1];
#line 240 "SB03MU.f"
	tmp[2] = sgn * tl[tl_dim1 + 1] * tr[tr_dim1 + 2];
#line 241 "SB03MU.f"
    }
#line 242 "SB03MU.f"
    btmp[0] = b[b_dim1 + 1];
#line 243 "SB03MU.f"
    btmp[1] = b[(b_dim1 << 1) + 1];
#line 244 "SB03MU.f"
    goto L40;

/*     2-by-1: */
/*     ISGN*op[TL11 TL12]*[X11]*TR11  = [B11]. */
/*            [TL21 TL22] [X21]         [B21] */

#line 250 "SB03MU.f"
L30:
/* Computing MAX */
/* Computing MAX */
#line 251 "SB03MU.f"
    d__7 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__8 = (d__2 = tl[(tl_dim1 << 
	    1) + 1], abs(d__2)), d__7 = max(d__7,d__8), d__8 = (d__3 = tl[
	    tl_dim1 + 2], abs(d__3)), d__7 = max(d__7,d__8), d__8 = (d__4 = 
	    tl[(tl_dim1 << 1) + 2], abs(d__4));
#line 251 "SB03MU.f"
    d__6 = max(d__7,d__8) * (d__5 = tr[tr_dim1 + 1], abs(d__5)) * eps;
#line 251 "SB03MU.f"
    smin = max(d__6,smlnum);
#line 255 "SB03MU.f"
    tmp[0] = sgn * tl[tl_dim1 + 1] * tr[tr_dim1 + 1] - 1.;
#line 256 "SB03MU.f"
    tmp[3] = sgn * tl[(tl_dim1 << 1) + 2] * tr[tr_dim1 + 1] - 1.;
#line 257 "SB03MU.f"
    if (*ltranl) {
#line 258 "SB03MU.f"
	tmp[1] = sgn * tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 1];
#line 259 "SB03MU.f"
	tmp[2] = sgn * tl[tl_dim1 + 2] * tr[tr_dim1 + 1];
#line 260 "SB03MU.f"
    } else {
#line 261 "SB03MU.f"
	tmp[1] = sgn * tl[tl_dim1 + 2] * tr[tr_dim1 + 1];
#line 262 "SB03MU.f"
	tmp[2] = sgn * tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 1];
#line 263 "SB03MU.f"
    }
#line 264 "SB03MU.f"
    btmp[0] = b[b_dim1 + 1];
#line 265 "SB03MU.f"
    btmp[1] = b[b_dim1 + 2];
#line 266 "SB03MU.f"
L40:

/*     Solve 2-by-2 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

#line 271 "SB03MU.f"
    ipiv = idamax_(&c__4, tmp, &c__1);
#line 272 "SB03MU.f"
    u11 = tmp[ipiv - 1];
#line 273 "SB03MU.f"
    if (abs(u11) <= smin) {
#line 274 "SB03MU.f"
	*info = 1;
#line 275 "SB03MU.f"
	u11 = smin;
#line 276 "SB03MU.f"
    }
#line 277 "SB03MU.f"
    u12 = tmp[locu12[ipiv - 1] - 1];
#line 278 "SB03MU.f"
    l21 = tmp[locl21[ipiv - 1] - 1] / u11;
#line 279 "SB03MU.f"
    u22 = tmp[locu22[ipiv - 1] - 1] - u12 * l21;
#line 280 "SB03MU.f"
    xswap = xswpiv[ipiv - 1];
#line 281 "SB03MU.f"
    bswap = bswpiv[ipiv - 1];
#line 282 "SB03MU.f"
    if (abs(u22) <= smin) {
#line 283 "SB03MU.f"
	*info = 1;
#line 284 "SB03MU.f"
	u22 = smin;
#line 285 "SB03MU.f"
    }
#line 286 "SB03MU.f"
    if (bswap) {
#line 287 "SB03MU.f"
	temp = btmp[1];
#line 288 "SB03MU.f"
	btmp[1] = btmp[0] - l21 * temp;
#line 289 "SB03MU.f"
	btmp[0] = temp;
#line 290 "SB03MU.f"
    } else {
#line 291 "SB03MU.f"
	btmp[1] -= l21 * btmp[0];
#line 292 "SB03MU.f"
    }
#line 293 "SB03MU.f"
    if (smlnum * 2. * abs(btmp[1]) > abs(u22) || smlnum * 2. * abs(btmp[0]) > 
	    abs(u11)) {
/* Computing MAX */
#line 295 "SB03MU.f"
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]);
#line 295 "SB03MU.f"
	*scale = .5 / max(d__1,d__2);
#line 296 "SB03MU.f"
	btmp[0] *= *scale;
#line 297 "SB03MU.f"
	btmp[1] *= *scale;
#line 298 "SB03MU.f"
    }
#line 299 "SB03MU.f"
    x2[1] = btmp[1] / u22;
#line 300 "SB03MU.f"
    x2[0] = btmp[0] / u11 - u12 / u11 * x2[1];
#line 301 "SB03MU.f"
    if (xswap) {
#line 302 "SB03MU.f"
	temp = x2[1];
#line 303 "SB03MU.f"
	x2[1] = x2[0];
#line 304 "SB03MU.f"
	x2[0] = temp;
#line 305 "SB03MU.f"
    }
#line 306 "SB03MU.f"
    x[x_dim1 + 1] = x2[0];
#line 307 "SB03MU.f"
    if (*n1 == 1) {
#line 308 "SB03MU.f"
	x[(x_dim1 << 1) + 1] = x2[1];
#line 309 "SB03MU.f"
	*xnorm = abs(x2[0]) + abs(x2[1]);
#line 310 "SB03MU.f"
    } else {
#line 311 "SB03MU.f"
	x[x_dim1 + 2] = x2[1];
/* Computing MAX */
#line 312 "SB03MU.f"
	d__1 = abs(x2[0]), d__2 = abs(x2[1]);
#line 312 "SB03MU.f"
	*xnorm = max(d__1,d__2);
#line 313 "SB03MU.f"
    }
#line 314 "SB03MU.f"
    return 0;

/*     2-by-2: */
/*     ISGN*op[TL11 TL12]*[X11 X12]*op[TR11 TR12]-[X11 X12] = [B11 B12]. */
/*            [TL21 TL22] [X21 X22]   [TR21 TR22] [X21 X22]   [B21 B22] */

/*     Solve equivalent 4-by-4 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

#line 323 "SB03MU.f"
L50:
/* Computing MAX */
#line 324 "SB03MU.f"
    d__5 = (d__1 = tr[tr_dim1 + 1], abs(d__1)), d__6 = (d__2 = tr[(tr_dim1 << 
	    1) + 1], abs(d__2)), d__5 = max(d__5,d__6), d__6 = (d__3 = tr[
	    tr_dim1 + 2], abs(d__3)), d__5 = max(d__5,d__6), d__6 = (d__4 = 
	    tr[(tr_dim1 << 1) + 2], abs(d__4));
#line 324 "SB03MU.f"
    smin = max(d__5,d__6);
/* Computing MAX */
#line 326 "SB03MU.f"
    d__5 = (d__1 = tl[tl_dim1 + 1], abs(d__1)), d__6 = (d__2 = tl[(tl_dim1 << 
	    1) + 1], abs(d__2)), d__5 = max(d__5,d__6), d__6 = (d__3 = tl[
	    tl_dim1 + 2], abs(d__3)), d__5 = max(d__5,d__6), d__6 = (d__4 = 
	    tl[(tl_dim1 << 1) + 2], abs(d__4));
#line 326 "SB03MU.f"
    smin = max(d__5,d__6) * smin;
/* Computing MAX */
#line 328 "SB03MU.f"
    d__1 = eps * smin;
#line 328 "SB03MU.f"
    smin = max(d__1,smlnum);
#line 329 "SB03MU.f"
    t16[0] = sgn * tl[tl_dim1 + 1] * tr[tr_dim1 + 1] - 1.;
#line 330 "SB03MU.f"
    t16[5] = sgn * tl[(tl_dim1 << 1) + 2] * tr[tr_dim1 + 1] - 1.;
#line 331 "SB03MU.f"
    t16[10] = sgn * tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 2] - 1.;
#line 332 "SB03MU.f"
    t16[15] = sgn * tl[(tl_dim1 << 1) + 2] * tr[(tr_dim1 << 1) + 2] - 1.;
#line 333 "SB03MU.f"
    if (*ltranl) {
#line 334 "SB03MU.f"
	t16[4] = sgn * tl[tl_dim1 + 2] * tr[tr_dim1 + 1];
#line 335 "SB03MU.f"
	t16[1] = sgn * tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 1];
#line 336 "SB03MU.f"
	t16[14] = sgn * tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 2];
#line 337 "SB03MU.f"
	t16[11] = sgn * tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 2];
#line 338 "SB03MU.f"
    } else {
#line 339 "SB03MU.f"
	t16[4] = sgn * tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 1];
#line 340 "SB03MU.f"
	t16[1] = sgn * tl[tl_dim1 + 2] * tr[tr_dim1 + 1];
#line 341 "SB03MU.f"
	t16[14] = sgn * tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 2];
#line 342 "SB03MU.f"
	t16[11] = sgn * tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 2];
#line 343 "SB03MU.f"
    }
#line 344 "SB03MU.f"
    if (*ltranr) {
#line 345 "SB03MU.f"
	t16[8] = sgn * tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 1];
#line 346 "SB03MU.f"
	t16[13] = sgn * tl[(tl_dim1 << 1) + 2] * tr[(tr_dim1 << 1) + 1];
#line 347 "SB03MU.f"
	t16[2] = sgn * tl[tl_dim1 + 1] * tr[tr_dim1 + 2];
#line 348 "SB03MU.f"
	t16[7] = sgn * tl[(tl_dim1 << 1) + 2] * tr[tr_dim1 + 2];
#line 349 "SB03MU.f"
    } else {
#line 350 "SB03MU.f"
	t16[8] = sgn * tl[tl_dim1 + 1] * tr[tr_dim1 + 2];
#line 351 "SB03MU.f"
	t16[13] = sgn * tl[(tl_dim1 << 1) + 2] * tr[tr_dim1 + 2];
#line 352 "SB03MU.f"
	t16[2] = sgn * tl[tl_dim1 + 1] * tr[(tr_dim1 << 1) + 1];
#line 353 "SB03MU.f"
	t16[7] = sgn * tl[(tl_dim1 << 1) + 2] * tr[(tr_dim1 << 1) + 1];
#line 354 "SB03MU.f"
    }
#line 355 "SB03MU.f"
    if (*ltranl && *ltranr) {
#line 356 "SB03MU.f"
	t16[12] = sgn * tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 1];
#line 357 "SB03MU.f"
	t16[9] = sgn * tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 1];
#line 358 "SB03MU.f"
	t16[6] = sgn * tl[tl_dim1 + 2] * tr[tr_dim1 + 2];
#line 359 "SB03MU.f"
	t16[3] = sgn * tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 2];
#line 360 "SB03MU.f"
    } else if (*ltranl && ! (*ltranr)) {
#line 361 "SB03MU.f"
	t16[12] = sgn * tl[tl_dim1 + 2] * tr[tr_dim1 + 2];
#line 362 "SB03MU.f"
	t16[9] = sgn * tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 2];
#line 363 "SB03MU.f"
	t16[6] = sgn * tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 1];
#line 364 "SB03MU.f"
	t16[3] = sgn * tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 1];
#line 365 "SB03MU.f"
    } else if (! (*ltranl) && *ltranr) {
#line 366 "SB03MU.f"
	t16[12] = sgn * tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 1];
#line 367 "SB03MU.f"
	t16[9] = sgn * tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 1];
#line 368 "SB03MU.f"
	t16[6] = sgn * tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 2];
#line 369 "SB03MU.f"
	t16[3] = sgn * tl[tl_dim1 + 2] * tr[tr_dim1 + 2];
#line 370 "SB03MU.f"
    } else {
#line 371 "SB03MU.f"
	t16[12] = sgn * tl[(tl_dim1 << 1) + 1] * tr[tr_dim1 + 2];
#line 372 "SB03MU.f"
	t16[9] = sgn * tl[tl_dim1 + 2] * tr[tr_dim1 + 2];
#line 373 "SB03MU.f"
	t16[6] = sgn * tl[(tl_dim1 << 1) + 1] * tr[(tr_dim1 << 1) + 1];
#line 374 "SB03MU.f"
	t16[3] = sgn * tl[tl_dim1 + 2] * tr[(tr_dim1 << 1) + 1];
#line 375 "SB03MU.f"
    }
#line 376 "SB03MU.f"
    btmp[0] = b[b_dim1 + 1];
#line 377 "SB03MU.f"
    btmp[1] = b[b_dim1 + 2];
#line 378 "SB03MU.f"
    btmp[2] = b[(b_dim1 << 1) + 1];
#line 379 "SB03MU.f"
    btmp[3] = b[(b_dim1 << 1) + 2];

/*     Perform elimination */

#line 383 "SB03MU.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 384 "SB03MU.f"
	xmax = 0.;

#line 386 "SB03MU.f"
	for (ip = i__; ip <= 4; ++ip) {

#line 388 "SB03MU.f"
	    for (jp = i__; jp <= 4; ++jp) {
#line 389 "SB03MU.f"
		if ((d__1 = t16[ip + (jp << 2) - 5], abs(d__1)) >= xmax) {
#line 390 "SB03MU.f"
		    xmax = (d__1 = t16[ip + (jp << 2) - 5], abs(d__1));
#line 391 "SB03MU.f"
		    ipsv = ip;
#line 392 "SB03MU.f"
		    jpsv = jp;
#line 393 "SB03MU.f"
		}
#line 394 "SB03MU.f"
/* L60: */
#line 394 "SB03MU.f"
	    }

#line 396 "SB03MU.f"
/* L70: */
#line 396 "SB03MU.f"
	}

#line 398 "SB03MU.f"
	if (ipsv != i__) {
#line 399 "SB03MU.f"
	    dswap_(&c__4, &t16[ipsv - 1], &c__4, &t16[i__ - 1], &c__4);
#line 400 "SB03MU.f"
	    temp = btmp[i__ - 1];
#line 401 "SB03MU.f"
	    btmp[i__ - 1] = btmp[ipsv - 1];
#line 402 "SB03MU.f"
	    btmp[ipsv - 1] = temp;
#line 403 "SB03MU.f"
	}
#line 404 "SB03MU.f"
	if (jpsv != i__) {
#line 404 "SB03MU.f"
	    dswap_(&c__4, &t16[(jpsv << 2) - 4], &c__1, &t16[(i__ << 2) - 4], 
		    &c__1);
#line 404 "SB03MU.f"
	}
#line 406 "SB03MU.f"
	jpiv[i__ - 1] = jpsv;
#line 407 "SB03MU.f"
	if ((d__1 = t16[i__ + (i__ << 2) - 5], abs(d__1)) < smin) {
#line 408 "SB03MU.f"
	    *info = 1;
#line 409 "SB03MU.f"
	    t16[i__ + (i__ << 2) - 5] = smin;
#line 410 "SB03MU.f"
	}

#line 412 "SB03MU.f"
	for (j = i__ + 1; j <= 4; ++j) {
#line 413 "SB03MU.f"
	    t16[j + (i__ << 2) - 5] /= t16[i__ + (i__ << 2) - 5];
#line 414 "SB03MU.f"
	    btmp[j - 1] -= t16[j + (i__ << 2) - 5] * btmp[i__ - 1];

#line 416 "SB03MU.f"
	    for (k = i__ + 1; k <= 4; ++k) {
#line 417 "SB03MU.f"
		t16[j + (k << 2) - 5] -= t16[j + (i__ << 2) - 5] * t16[i__ + (
			k << 2) - 5];
#line 418 "SB03MU.f"
/* L80: */
#line 418 "SB03MU.f"
	    }

#line 420 "SB03MU.f"
/* L90: */
#line 420 "SB03MU.f"
	}

#line 422 "SB03MU.f"
/* L100: */
#line 422 "SB03MU.f"
    }

#line 424 "SB03MU.f"
    if (abs(t16[15]) < smin) {
#line 424 "SB03MU.f"
	t16[15] = smin;
#line 424 "SB03MU.f"
    }
#line 426 "SB03MU.f"
    if (smlnum * 8. * abs(btmp[0]) > abs(t16[0]) || smlnum * 8. * abs(btmp[1])
	     > abs(t16[5]) || smlnum * 8. * abs(btmp[2]) > abs(t16[10]) || 
	    smlnum * 8. * abs(btmp[3]) > abs(t16[15])) {
/* Computing MAX */
#line 430 "SB03MU.f"
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]), d__1 = max(d__1,d__2), d__2 
		= abs(btmp[2]), d__1 = max(d__1,d__2), d__2 = abs(btmp[3]);
#line 430 "SB03MU.f"
	*scale = .125 / max(d__1,d__2);
#line 433 "SB03MU.f"
	btmp[0] *= *scale;
#line 434 "SB03MU.f"
	btmp[1] *= *scale;
#line 435 "SB03MU.f"
	btmp[2] *= *scale;
#line 436 "SB03MU.f"
	btmp[3] *= *scale;
#line 437 "SB03MU.f"
    }

#line 439 "SB03MU.f"
    for (i__ = 1; i__ <= 4; ++i__) {
#line 440 "SB03MU.f"
	k = 5 - i__;
#line 441 "SB03MU.f"
	temp = 1. / t16[k + (k << 2) - 5];
#line 442 "SB03MU.f"
	tmp[k - 1] = btmp[k - 1] * temp;

#line 444 "SB03MU.f"
	for (j = k + 1; j <= 4; ++j) {
#line 445 "SB03MU.f"
	    tmp[k - 1] -= temp * t16[k + (j << 2) - 5] * tmp[j - 1];
#line 446 "SB03MU.f"
/* L110: */
#line 446 "SB03MU.f"
	}

#line 448 "SB03MU.f"
/* L120: */
#line 448 "SB03MU.f"
    }

#line 450 "SB03MU.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 451 "SB03MU.f"
	if (jpiv[4 - i__ - 1] != 4 - i__) {
#line 452 "SB03MU.f"
	    temp = tmp[4 - i__ - 1];
#line 453 "SB03MU.f"
	    tmp[4 - i__ - 1] = tmp[jpiv[4 - i__ - 1] - 1];
#line 454 "SB03MU.f"
	    tmp[jpiv[4 - i__ - 1] - 1] = temp;
#line 455 "SB03MU.f"
	}
#line 456 "SB03MU.f"
/* L130: */
#line 456 "SB03MU.f"
    }

#line 458 "SB03MU.f"
    x[x_dim1 + 1] = tmp[0];
#line 459 "SB03MU.f"
    x[x_dim1 + 2] = tmp[1];
#line 460 "SB03MU.f"
    x[(x_dim1 << 1) + 1] = tmp[2];
#line 461 "SB03MU.f"
    x[(x_dim1 << 1) + 2] = tmp[3];
/* Computing MAX */
#line 462 "SB03MU.f"
    d__1 = abs(tmp[0]) + abs(tmp[2]), d__2 = abs(tmp[1]) + abs(tmp[3]);
#line 462 "SB03MU.f"
    *xnorm = max(d__1,d__2);

#line 465 "SB03MU.f"
    return 0;
/* *** Last line of SB03MU *** */
} /* sb03mu_ */

