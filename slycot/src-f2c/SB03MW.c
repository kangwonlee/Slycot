#line 1 "SB03MW.f"
/* SB03MW.f -- translated by f2c (version 20100827).
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

#line 1 "SB03MW.f"
/* Table of constant values */

static integer c__3 = 3;
static integer c__1 = 1;

/* Subroutine */ int sb03mw_(logical *ltran, logical *lupper, doublereal *t, 
	integer *ldt, doublereal *b, integer *ldb, doublereal *scale, 
	doublereal *x, integer *ldx, doublereal *xnorm, integer *info)
{
    /* System generated locals */
    integer b_dim1, b_offset, t_dim1, t_offset, x_dim1, x_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Local variables */
    static integer i__, j, k;
    static doublereal t9[9]	/* was [3][3] */;
    static integer ip, jp;
    static doublereal eps, tmp[3], btmp[3], temp, smin;
    static integer jpiv[3];
    static doublereal xmax;
    static integer ipsv, jpsv;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
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

/*     To solve for the 2-by-2 symmetric matrix X in */

/*            op(T)'*X + X*op(T) = SCALE*B, */

/*     where T is 2-by-2, B is symmetric 2-by-2, and op(T) = T or T', */
/*     where T' denotes the transpose of T. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     LTRAN   LOGICAL */
/*             Specifies the form of op(T) to be used, as follows: */
/*             = .FALSE.:  op(T) = T, */
/*             = .TRUE. :  op(T) = T'. */

/*     LUPPER  LOGICAL */
/*             Specifies which triangle of the matrix B is used, and */
/*             which triangle of the matrix X is computed, as follows: */
/*             = .TRUE. :  The upper triangular part; */
/*             = .FALSE.:  The lower triangular part. */

/*     Input/Output Parameters */

/*     T       (input) DOUBLE PRECISION array, dimension (LDT,2) */
/*             The leading 2-by-2 part of this array must contain the */
/*             matrix T. */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= 2. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,2) */
/*             On entry with LUPPER = .TRUE., the leading 2-by-2 upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the symmetric matrix B and the strictly */
/*             lower triangular part of B is not referenced. */
/*             On entry with LUPPER = .FALSE., the leading 2-by-2 lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the symmetric matrix B and the strictly */
/*             upper triangular part of B is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= 2. */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor. SCALE is chosen less than or equal to 1 */
/*             to prevent the solution overflowing. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,2) */
/*             On exit with LUPPER = .TRUE., the leading 2-by-2 upper */
/*             triangular part of this array contains the upper */
/*             triangular part of the symmetric solution matrix X and the */
/*             strictly lower triangular part of X is not referenced. */
/*             On exit with LUPPER = .FALSE., the leading 2-by-2 lower */
/*             triangular part of this array contains the lower */
/*             triangular part of the symmetric solution matrix X and the */
/*             strictly upper triangular part of X is not referenced. */
/*             Note that X may be identified with B in the calling */
/*             statement. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= 2. */

/*     XNORM   (output) DOUBLE PRECISION */
/*             The infinity-norm of the solution. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if T and -T have too close eigenvalues, so T */
/*                   is perturbed to get a nonsingular equation. */

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
/*     Based on DLALY2 by P. Petkov, Tech. University of Sofia, September */
/*     1993. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Continuous-time system, Lyapunov equation, matrix algebra. */

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
/*     .. Executable Statements .. */

/*     Do not check the input parameters for errors */

#line 169 "SB03MW.f"
    /* Parameter adjustments */
#line 169 "SB03MW.f"
    t_dim1 = *ldt;
#line 169 "SB03MW.f"
    t_offset = 1 + t_dim1;
#line 169 "SB03MW.f"
    t -= t_offset;
#line 169 "SB03MW.f"
    b_dim1 = *ldb;
#line 169 "SB03MW.f"
    b_offset = 1 + b_dim1;
#line 169 "SB03MW.f"
    b -= b_offset;
#line 169 "SB03MW.f"
    x_dim1 = *ldx;
#line 169 "SB03MW.f"
    x_offset = 1 + x_dim1;
#line 169 "SB03MW.f"
    x -= x_offset;
#line 169 "SB03MW.f"

#line 169 "SB03MW.f"
    /* Function Body */
#line 169 "SB03MW.f"
    *info = 0;

/*     Set constants to control overflow */

#line 173 "SB03MW.f"
    eps = dlamch_("P", (ftnlen)1);
#line 174 "SB03MW.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;

/*     Solve equivalent 3-by-3 system using complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

/* Computing MAX */
/* Computing MAX */
#line 179 "SB03MW.f"
    d__6 = (d__1 = t[t_dim1 + 1], abs(d__1)), d__7 = (d__2 = t[(t_dim1 << 1) 
	    + 1], abs(d__2)), d__6 = max(d__6,d__7), d__7 = (d__3 = t[t_dim1 
	    + 2], abs(d__3)), d__6 = max(d__6,d__7), d__7 = (d__4 = t[(t_dim1 
	    << 1) + 2], abs(d__4));
#line 179 "SB03MW.f"
    d__5 = max(d__6,d__7) * eps;
#line 179 "SB03MW.f"
    smin = max(d__5,smlnum);
#line 182 "SB03MW.f"
    t9[6] = 0.;
#line 183 "SB03MW.f"
    t9[2] = 0.;
#line 184 "SB03MW.f"
    t9[0] = t[t_dim1 + 1];
#line 185 "SB03MW.f"
    t9[4] = t[t_dim1 + 1] + t[(t_dim1 << 1) + 2];
#line 186 "SB03MW.f"
    t9[8] = t[(t_dim1 << 1) + 2];
#line 187 "SB03MW.f"
    if (*ltran) {
#line 188 "SB03MW.f"
	t9[3] = t[(t_dim1 << 1) + 1];
#line 189 "SB03MW.f"
	t9[1] = t[t_dim1 + 2];
#line 190 "SB03MW.f"
	t9[7] = t[(t_dim1 << 1) + 1];
#line 191 "SB03MW.f"
	t9[5] = t[t_dim1 + 2];
#line 192 "SB03MW.f"
    } else {
#line 193 "SB03MW.f"
	t9[3] = t[t_dim1 + 2];
#line 194 "SB03MW.f"
	t9[1] = t[(t_dim1 << 1) + 1];
#line 195 "SB03MW.f"
	t9[7] = t[t_dim1 + 2];
#line 196 "SB03MW.f"
	t9[5] = t[(t_dim1 << 1) + 1];
#line 197 "SB03MW.f"
    }
#line 198 "SB03MW.f"
    btmp[0] = b[b_dim1 + 1] / 2.;
#line 199 "SB03MW.f"
    if (*lupper) {
#line 200 "SB03MW.f"
	btmp[1] = b[(b_dim1 << 1) + 1];
#line 201 "SB03MW.f"
    } else {
#line 202 "SB03MW.f"
	btmp[1] = b[b_dim1 + 2];
#line 203 "SB03MW.f"
    }
#line 204 "SB03MW.f"
    btmp[2] = b[(b_dim1 << 1) + 2] / 2.;

/*     Perform elimination */

#line 208 "SB03MW.f"
    for (i__ = 1; i__ <= 2; ++i__) {
#line 209 "SB03MW.f"
	xmax = 0.;

#line 211 "SB03MW.f"
	for (ip = i__; ip <= 3; ++ip) {

#line 213 "SB03MW.f"
	    for (jp = i__; jp <= 3; ++jp) {
#line 214 "SB03MW.f"
		if ((d__1 = t9[ip + jp * 3 - 4], abs(d__1)) >= xmax) {
#line 215 "SB03MW.f"
		    xmax = (d__1 = t9[ip + jp * 3 - 4], abs(d__1));
#line 216 "SB03MW.f"
		    ipsv = ip;
#line 217 "SB03MW.f"
		    jpsv = jp;
#line 218 "SB03MW.f"
		}
#line 219 "SB03MW.f"
/* L10: */
#line 219 "SB03MW.f"
	    }

#line 221 "SB03MW.f"
/* L20: */
#line 221 "SB03MW.f"
	}

#line 223 "SB03MW.f"
	if (ipsv != i__) {
#line 224 "SB03MW.f"
	    dswap_(&c__3, &t9[ipsv - 1], &c__3, &t9[i__ - 1], &c__3);
#line 225 "SB03MW.f"
	    temp = btmp[i__ - 1];
#line 226 "SB03MW.f"
	    btmp[i__ - 1] = btmp[ipsv - 1];
#line 227 "SB03MW.f"
	    btmp[ipsv - 1] = temp;
#line 228 "SB03MW.f"
	}
#line 229 "SB03MW.f"
	if (jpsv != i__) {
#line 229 "SB03MW.f"
	    dswap_(&c__3, &t9[jpsv * 3 - 3], &c__1, &t9[i__ * 3 - 3], &c__1);
#line 229 "SB03MW.f"
	}
#line 231 "SB03MW.f"
	jpiv[i__ - 1] = jpsv;
#line 232 "SB03MW.f"
	if ((d__1 = t9[i__ + i__ * 3 - 4], abs(d__1)) < smin) {
#line 233 "SB03MW.f"
	    *info = 1;
#line 234 "SB03MW.f"
	    t9[i__ + i__ * 3 - 4] = smin;
#line 235 "SB03MW.f"
	}

#line 237 "SB03MW.f"
	for (j = i__ + 1; j <= 3; ++j) {
#line 238 "SB03MW.f"
	    t9[j + i__ * 3 - 4] /= t9[i__ + i__ * 3 - 4];
#line 239 "SB03MW.f"
	    btmp[j - 1] -= t9[j + i__ * 3 - 4] * btmp[i__ - 1];

#line 241 "SB03MW.f"
	    for (k = i__ + 1; k <= 3; ++k) {
#line 242 "SB03MW.f"
		t9[j + k * 3 - 4] -= t9[j + i__ * 3 - 4] * t9[i__ + k * 3 - 4]
			;
#line 243 "SB03MW.f"
/* L30: */
#line 243 "SB03MW.f"
	    }

#line 245 "SB03MW.f"
/* L40: */
#line 245 "SB03MW.f"
	}

#line 247 "SB03MW.f"
/* L50: */
#line 247 "SB03MW.f"
    }

#line 249 "SB03MW.f"
    if (abs(t9[8]) < smin) {
#line 249 "SB03MW.f"
	t9[8] = smin;
#line 249 "SB03MW.f"
    }
#line 251 "SB03MW.f"
    *scale = 1.;
#line 252 "SB03MW.f"
    if (smlnum * 4. * abs(btmp[0]) > abs(t9[0]) || smlnum * 4. * abs(btmp[1]) 
	    > abs(t9[4]) || smlnum * 4. * abs(btmp[2]) > abs(t9[8])) {
/* Computing MAX */
#line 255 "SB03MW.f"
	d__1 = abs(btmp[0]), d__2 = abs(btmp[1]), d__1 = max(d__1,d__2), d__2 
		= abs(btmp[2]);
#line 255 "SB03MW.f"
	*scale = .25 / max(d__1,d__2);
#line 257 "SB03MW.f"
	btmp[0] *= *scale;
#line 258 "SB03MW.f"
	btmp[1] *= *scale;
#line 259 "SB03MW.f"
	btmp[2] *= *scale;
#line 260 "SB03MW.f"
    }

#line 262 "SB03MW.f"
    for (i__ = 1; i__ <= 3; ++i__) {
#line 263 "SB03MW.f"
	k = 4 - i__;
#line 264 "SB03MW.f"
	temp = 1. / t9[k + k * 3 - 4];
#line 265 "SB03MW.f"
	tmp[k - 1] = btmp[k - 1] * temp;

#line 267 "SB03MW.f"
	for (j = k + 1; j <= 3; ++j) {
#line 268 "SB03MW.f"
	    tmp[k - 1] -= temp * t9[k + j * 3 - 4] * tmp[j - 1];
#line 269 "SB03MW.f"
/* L60: */
#line 269 "SB03MW.f"
	}

#line 271 "SB03MW.f"
/* L70: */
#line 271 "SB03MW.f"
    }

#line 273 "SB03MW.f"
    for (i__ = 1; i__ <= 2; ++i__) {
#line 274 "SB03MW.f"
	if (jpiv[3 - i__ - 1] != 3 - i__) {
#line 275 "SB03MW.f"
	    temp = tmp[3 - i__ - 1];
#line 276 "SB03MW.f"
	    tmp[3 - i__ - 1] = tmp[jpiv[3 - i__ - 1] - 1];
#line 277 "SB03MW.f"
	    tmp[jpiv[3 - i__ - 1] - 1] = temp;
#line 278 "SB03MW.f"
	}
#line 279 "SB03MW.f"
/* L80: */
#line 279 "SB03MW.f"
    }

#line 281 "SB03MW.f"
    x[x_dim1 + 1] = tmp[0];
#line 282 "SB03MW.f"
    if (*lupper) {
#line 283 "SB03MW.f"
	x[(x_dim1 << 1) + 1] = tmp[1];
#line 284 "SB03MW.f"
    } else {
#line 285 "SB03MW.f"
	x[x_dim1 + 2] = tmp[1];
#line 286 "SB03MW.f"
    }
#line 287 "SB03MW.f"
    x[(x_dim1 << 1) + 2] = tmp[2];
/* Computing MAX */
#line 288 "SB03MW.f"
    d__1 = abs(tmp[0]) + abs(tmp[1]), d__2 = abs(tmp[1]) + abs(tmp[2]);
#line 288 "SB03MW.f"
    *xnorm = max(d__1,d__2);

#line 291 "SB03MW.f"
    return 0;
/* *** Last line of SB03MW *** */
} /* sb03mw_ */

