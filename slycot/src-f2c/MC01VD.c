#line 1 "MC01VD.f"
/* MC01VD.f -- translated by f2c (version 20100827).
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

#line 1 "MC01VD.f"
/* Table of constant values */

static doublereal c_b4 = 1.;

/* Subroutine */ int mc01vd_(doublereal *a, doublereal *b, doublereal *c__, 
	doublereal *z1re, doublereal *z1im, doublereal *z2re, doublereal *
	z2im, integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal w, m1, m2;
    static integer ea, eb, ec, ed;
    static doublereal ma, mb, mc, md;
    static integer eb2;
    static doublereal big, absa, absb, absc;
    static integer beta;
    extern /* Subroutine */ int mc01sw_(doublereal *, integer *, doublereal *,
	     integer *);
    static doublereal sfmin;
    extern /* Subroutine */ int mc01sy_(doublereal *, integer *, integer *, 
	    doublereal *, logical *);
    extern doublereal dlamch_(char *, ftnlen);
    static integer eaplec;
    static logical ovflow;


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

/*     To compute the roots of a quadratic equation with real */
/*     coefficients. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     A       (input) DOUBLE PRECISION */
/*             The value of the coefficient of the quadratic term. */

/*     B       (input) DOUBLE PRECISION */
/*             The value of the coefficient of the linear term. */

/*     C       (input) DOUBLE PRECISION */
/*             The value of the coefficient of the constant term. */

/*     Z1RE    (output) DOUBLE PRECISION */
/*     Z1IM    (output) DOUBLE PRECISION */
/*             The real and imaginary parts, respectively, of the largest */
/*             root in magnitude. */

/*     Z2RE    (output) DOUBLE PRECISION */
/*     Z2IM    (output) DOUBLE PRECISION */
/*             The real and imaginary parts, respectively, of the */
/*             smallest root in magnitude. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if on entry, either A = B = 0.0 or A = 0.0 and the */
/*                   root -C/B overflows; in this case Z1RE, Z1IM, Z2RE */
/*                   and Z2IM are unassigned; */
/*             = 2:  if on entry, A = 0.0; in this case Z1RE contains */
/*                   BIG and Z1IM contains zero, where BIG is a */
/*                   representable number near the overflow threshold */
/*                   of the machine (see LAPACK Library Routine DLAMCH); */
/*             = 3:  if on entry, either C = 0.0 and the root -B/A */
/*                   overflows or A, B and C are non-zero and the largest */
/*                   real root in magnitude cannot be computed without */
/*                   overflow; in this case Z1RE contains BIG and Z1IM */
/*                   contains zero; */
/*             = 4:  if the roots cannot be computed without overflow; in */
/*                   this case Z1RE, Z1IM, Z2RE and Z2IM are unassigned. */

/*     METHOD */

/*     The routine computes the roots (r1 and r2) of the real quadratic */
/*     equation */
/*             2 */
/*        a * x  + b * x + c = 0 */

/*     as */
/*             - b - SIGN(b) * SQRT(b * b - 4 * a * c)             c */
/*        r1 = ---------------------------------------  and r2 = ------ */
/*                              2 * a                            a * r1 */

/*     unless a = 0, in which case */

/*             -c */
/*        r1 = --. */
/*              b */

/*     Precautions are taken to avoid overflow and underflow wherever */
/*     possible. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MC01JD by A.J. Geurts. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Quadratic equation, zeros. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Detect special cases. */

#line 130 "MC01VD.f"
    *info = 0;
#line 131 "MC01VD.f"
    beta = (integer) dlamch_("Base", (ftnlen)4);
#line 132 "MC01VD.f"
    sfmin = dlamch_("Safe minimum", (ftnlen)12);
#line 133 "MC01VD.f"
    big = 1. / sfmin;
#line 134 "MC01VD.f"
    if (*a == 0.) {
#line 135 "MC01VD.f"
	if (*b == 0.) {
#line 136 "MC01VD.f"
	    *info = 1;
#line 137 "MC01VD.f"
	} else {
#line 138 "MC01VD.f"
	    ovflow = FALSE_;
#line 139 "MC01VD.f"
	    *z2re = 0.;
#line 140 "MC01VD.f"
	    if (*c__ != 0.) {
#line 141 "MC01VD.f"
		absb = abs(*b);
#line 142 "MC01VD.f"
		if (absb >= 1.) {
#line 143 "MC01VD.f"
		    if (abs(*c__) >= absb * sfmin) {
#line 143 "MC01VD.f"
			*z2re = -(*c__) / *b;
#line 143 "MC01VD.f"
		    }
#line 144 "MC01VD.f"
		} else {
#line 145 "MC01VD.f"
		    if (abs(*c__) <= absb * big) {
#line 146 "MC01VD.f"
			*z2re = -(*c__) / *b;
#line 147 "MC01VD.f"
		    } else {
#line 148 "MC01VD.f"
			ovflow = TRUE_;
#line 149 "MC01VD.f"
			*z2re = big;
#line 150 "MC01VD.f"
			if (d_sign(&c_b4, b) * d_sign(&c_b4, c__) > 0.) {
#line 150 "MC01VD.f"
			    *z2re = -big;
#line 150 "MC01VD.f"
			}
#line 152 "MC01VD.f"
		    }
#line 153 "MC01VD.f"
		}
#line 154 "MC01VD.f"
	    }
#line 155 "MC01VD.f"
	    if (ovflow) {
#line 156 "MC01VD.f"
		*info = 1;
#line 157 "MC01VD.f"
	    } else {
#line 158 "MC01VD.f"
		*z1re = big;
#line 159 "MC01VD.f"
		*z1im = 0.;
#line 160 "MC01VD.f"
		*z2im = 0.;
#line 161 "MC01VD.f"
		*info = 2;
#line 162 "MC01VD.f"
	    }
#line 163 "MC01VD.f"
	}
#line 164 "MC01VD.f"
	return 0;
#line 165 "MC01VD.f"
    }

#line 167 "MC01VD.f"
    if (*c__ == 0.) {
#line 168 "MC01VD.f"
	ovflow = FALSE_;
#line 169 "MC01VD.f"
	*z1re = 0.;
#line 170 "MC01VD.f"
	if (*b != 0.) {
#line 171 "MC01VD.f"
	    absa = abs(*a);
#line 172 "MC01VD.f"
	    if (absa >= 1.) {
#line 173 "MC01VD.f"
		if (abs(*b) >= absa * sfmin) {
#line 173 "MC01VD.f"
		    *z1re = -(*b) / *a;
#line 173 "MC01VD.f"
		}
#line 174 "MC01VD.f"
	    } else {
#line 175 "MC01VD.f"
		if (abs(*b) <= absa * big) {
#line 176 "MC01VD.f"
		    *z1re = -(*b) / *a;
#line 177 "MC01VD.f"
		} else {
#line 178 "MC01VD.f"
		    ovflow = TRUE_;
#line 179 "MC01VD.f"
		    *z1re = big;
#line 180 "MC01VD.f"
		}
#line 181 "MC01VD.f"
	    }
#line 182 "MC01VD.f"
	}
#line 183 "MC01VD.f"
	if (ovflow) {
#line 183 "MC01VD.f"
	    *info = 3;
#line 183 "MC01VD.f"
	}
#line 184 "MC01VD.f"
	*z1im = 0.;
#line 185 "MC01VD.f"
	*z2re = 0.;
#line 186 "MC01VD.f"
	*z2im = 0.;
#line 187 "MC01VD.f"
	return 0;
#line 188 "MC01VD.f"
    }

/*     A and C are non-zero. */

#line 192 "MC01VD.f"
    if (*b == 0.) {
#line 193 "MC01VD.f"
	ovflow = FALSE_;
#line 194 "MC01VD.f"
	absc = sqrt((abs(*c__)));
#line 195 "MC01VD.f"
	absa = sqrt((abs(*a)));
#line 196 "MC01VD.f"
	w = 0.;
#line 197 "MC01VD.f"
	if (absa >= 1.) {
#line 198 "MC01VD.f"
	    if (absc >= absa * sfmin) {
#line 198 "MC01VD.f"
		w = absc / absa;
#line 198 "MC01VD.f"
	    }
#line 199 "MC01VD.f"
	} else {
#line 200 "MC01VD.f"
	    if (absc <= absa * big) {
#line 201 "MC01VD.f"
		w = absc / absa;
#line 202 "MC01VD.f"
	    } else {
#line 203 "MC01VD.f"
		ovflow = TRUE_;
#line 204 "MC01VD.f"
		w = big;
#line 205 "MC01VD.f"
	    }
#line 206 "MC01VD.f"
	}
#line 207 "MC01VD.f"
	if (ovflow) {
#line 208 "MC01VD.f"
	    *info = 4;
#line 209 "MC01VD.f"
	} else {
#line 210 "MC01VD.f"
	    if (d_sign(&c_b4, a) * d_sign(&c_b4, c__) > 0.) {
#line 211 "MC01VD.f"
		*z1re = 0.;
#line 212 "MC01VD.f"
		*z2re = 0.;
#line 213 "MC01VD.f"
		*z1im = w;
#line 214 "MC01VD.f"
		*z2im = -w;
#line 215 "MC01VD.f"
	    } else {
#line 216 "MC01VD.f"
		*z1re = w;
#line 217 "MC01VD.f"
		*z2re = -w;
#line 218 "MC01VD.f"
		*z1im = 0.;
#line 219 "MC01VD.f"
		*z2im = 0.;
#line 220 "MC01VD.f"
	    }
#line 221 "MC01VD.f"
	}
#line 222 "MC01VD.f"
	return 0;
#line 223 "MC01VD.f"
    }

/*     A, B and C are non-zero. */

#line 227 "MC01VD.f"
    mc01sw_(a, &beta, &ma, &ea);
#line 228 "MC01VD.f"
    mc01sw_(b, &beta, &mb, &eb);
#line 229 "MC01VD.f"
    mc01sw_(c__, &beta, &mc, &ec);

/*     Compute a 'near' floating-point representation of the discriminant */
/*     D = MD * BETA**ED. */

#line 234 "MC01VD.f"
    eaplec = ea + ec;
#line 235 "MC01VD.f"
    eb2 = eb << 1;
#line 236 "MC01VD.f"
    if (eaplec > eb2) {
#line 237 "MC01VD.f"
	d__1 = mb * mb;
#line 237 "MC01VD.f"
	i__1 = eb2 - eaplec;
#line 237 "MC01VD.f"
	mc01sy_(&d__1, &i__1, &beta, &w, &ovflow);
#line 238 "MC01VD.f"
	w -= ma * 4. * mc;
#line 239 "MC01VD.f"
	mc01sw_(&w, &beta, &md, &ed);
#line 240 "MC01VD.f"
	ed += eaplec;
#line 241 "MC01VD.f"
    } else {
#line 242 "MC01VD.f"
	d__1 = ma * 4. * mc;
#line 242 "MC01VD.f"
	i__1 = eaplec - eb2;
#line 242 "MC01VD.f"
	mc01sy_(&d__1, &i__1, &beta, &w, &ovflow);
#line 243 "MC01VD.f"
	w = mb * mb - w;
#line 244 "MC01VD.f"
	mc01sw_(&w, &beta, &md, &ed);
#line 245 "MC01VD.f"
	ed += eb2;
#line 246 "MC01VD.f"
    }

#line 248 "MC01VD.f"
    if (ed % 2 != 0) {
#line 249 "MC01VD.f"
	++ed;
#line 250 "MC01VD.f"
	md /= beta;
#line 251 "MC01VD.f"
    }

/*     Complex roots. */

#line 255 "MC01VD.f"
    if (md < 0.) {
#line 256 "MC01VD.f"
	d__1 = -mb / (ma * 2);
#line 256 "MC01VD.f"
	i__1 = eb - ea;
#line 256 "MC01VD.f"
	mc01sy_(&d__1, &i__1, &beta, z1re, &ovflow);
#line 257 "MC01VD.f"
	if (ovflow) {
#line 258 "MC01VD.f"
	    *info = 4;
#line 259 "MC01VD.f"
	} else {
#line 260 "MC01VD.f"
	    d__1 = sqrt(-md) / (ma * 2);
#line 260 "MC01VD.f"
	    i__1 = ed / 2 - ea;
#line 260 "MC01VD.f"
	    mc01sy_(&d__1, &i__1, &beta, z1im, &ovflow);
#line 262 "MC01VD.f"
	    if (ovflow) {
#line 263 "MC01VD.f"
		*info = 4;
#line 264 "MC01VD.f"
	    } else {
#line 265 "MC01VD.f"
		*z2re = *z1re;
#line 266 "MC01VD.f"
		*z2im = -(*z1im);
#line 267 "MC01VD.f"
	    }
#line 268 "MC01VD.f"
	}
#line 269 "MC01VD.f"
	return 0;
#line 270 "MC01VD.f"
    }

/*     Real roots. */

#line 274 "MC01VD.f"
    md = sqrt(md);
#line 275 "MC01VD.f"
    ed /= 2;
#line 276 "MC01VD.f"
    if (ed > eb) {
#line 277 "MC01VD.f"
	d__1 = abs(mb);
#line 277 "MC01VD.f"
	i__1 = eb - ed;
#line 277 "MC01VD.f"
	mc01sy_(&d__1, &i__1, &beta, &w, &ovflow);
#line 278 "MC01VD.f"
	w += md;
#line 279 "MC01VD.f"
	m1 = -d_sign(&c_b4, &mb) * w / (ma * 2);
#line 280 "MC01VD.f"
	i__1 = ed - ea;
#line 280 "MC01VD.f"
	mc01sy_(&m1, &i__1, &beta, z1re, &ovflow);
#line 281 "MC01VD.f"
	if (ovflow) {
#line 282 "MC01VD.f"
	    *z1re = big;
#line 283 "MC01VD.f"
	    *info = 3;
#line 284 "MC01VD.f"
	}
#line 285 "MC01VD.f"
	m2 = -d_sign(&c_b4, &mb) * 2 * mc / w;
#line 286 "MC01VD.f"
	i__1 = ec - ed;
#line 286 "MC01VD.f"
	mc01sy_(&m2, &i__1, &beta, z2re, &ovflow);
#line 287 "MC01VD.f"
    } else {
#line 288 "MC01VD.f"
	i__1 = ed - eb;
#line 288 "MC01VD.f"
	mc01sy_(&md, &i__1, &beta, &w, &ovflow);
#line 289 "MC01VD.f"
	w += abs(mb);
#line 290 "MC01VD.f"
	m1 = -d_sign(&c_b4, &mb) * w / (ma * 2);
#line 291 "MC01VD.f"
	i__1 = eb - ea;
#line 291 "MC01VD.f"
	mc01sy_(&m1, &i__1, &beta, z1re, &ovflow);
#line 292 "MC01VD.f"
	if (ovflow) {
#line 293 "MC01VD.f"
	    *z1re = big;
#line 294 "MC01VD.f"
	    *info = 3;
#line 295 "MC01VD.f"
	}
#line 296 "MC01VD.f"
	m2 = -d_sign(&c_b4, &mb) * 2 * mc / w;
#line 297 "MC01VD.f"
	i__1 = ec - eb;
#line 297 "MC01VD.f"
	mc01sy_(&m2, &i__1, &beta, z2re, &ovflow);
#line 298 "MC01VD.f"
    }
#line 299 "MC01VD.f"
    *z1im = 0.;
#line 300 "MC01VD.f"
    *z2im = 0.;

#line 302 "MC01VD.f"
    return 0;
/* *** Last line of MC01VD *** */
} /* mc01vd_ */

