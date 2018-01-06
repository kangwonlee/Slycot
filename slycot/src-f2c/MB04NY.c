#line 1 "MB04NY.f"
/* MB04NY.f -- translated by f2c (version 20100827).
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

#line 1 "MB04NY.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b15 = 1.;

/* Subroutine */ int mb04ny_(integer *m, integer *n, doublereal *v, integer *
	incv, doublereal *tau, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *dwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, 
	    v6, v7, v8, v9;
    static integer iv;
    static doublereal sum;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), daxpy_(integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *);


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

/*     To apply a real elementary reflector H to a real m-by-(n+1) */
/*     matrix C = [ A  B ], from the right, where A has one column. H is */
/*     represented in the form */
/*                                        ( 1 ) */
/*           H = I - tau * u *u',    u  = (   ), */
/*                                        ( v ) */
/*     where tau is a real scalar and v is a real n-vector. */

/*     If tau = 0, then H is taken to be the unit matrix. */

/*     In-line code is used if H has order < 11. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices A and B.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix B.  N >= 0. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*             (1+(N-1)*ABS( INCV )) */
/*             The vector v in the representation of H. */

/*     INCV    (input) INTEGER */
/*             The increment between the elements of v.  INCV <> 0. */

/*     TAU     (input) DOUBLE PRECISION */
/*             The scalar factor of the elementary reflector H. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,1) */
/*             On entry, the leading M-by-1 part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading M-by-1 part of this array contains */
/*             the updated matrix A (the first column of C * H). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the updated matrix B (the last n columns of C * H). */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (M) */
/*             DWORK is not referenced if H has order less than 11. */

/*     METHOD */

/*     The routine applies the elementary reflector H, taking the special */
/*     structure of C into account. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1998. */
/*     Based on LAPACK routines DLARFX and DLATZM. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, elementary reflector, orthogonal */
/*     transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

#line 120 "MB04NY.f"
    /* Parameter adjustments */
#line 120 "MB04NY.f"
    --v;
#line 120 "MB04NY.f"
    a_dim1 = *lda;
#line 120 "MB04NY.f"
    a_offset = 1 + a_dim1;
#line 120 "MB04NY.f"
    a -= a_offset;
#line 120 "MB04NY.f"
    b_dim1 = *ldb;
#line 120 "MB04NY.f"
    b_offset = 1 + b_dim1;
#line 120 "MB04NY.f"
    b -= b_offset;
#line 120 "MB04NY.f"
    --dwork;
#line 120 "MB04NY.f"

#line 120 "MB04NY.f"
    /* Function Body */
#line 120 "MB04NY.f"
    if (*tau == 0.) {
#line 120 "MB04NY.f"
	return 0;
#line 120 "MB04NY.f"
    }

/*     Form  C * H, where H has order n+1. */

#line 125 "MB04NY.f"
    switch (*n + 1) {
#line 125 "MB04NY.f"
	case 1:  goto L10;
#line 125 "MB04NY.f"
	case 2:  goto L30;
#line 125 "MB04NY.f"
	case 3:  goto L50;
#line 125 "MB04NY.f"
	case 4:  goto L70;
#line 125 "MB04NY.f"
	case 5:  goto L90;
#line 125 "MB04NY.f"
	case 6:  goto L110;
#line 125 "MB04NY.f"
	case 7:  goto L130;
#line 125 "MB04NY.f"
	case 8:  goto L150;
#line 125 "MB04NY.f"
	case 9:  goto L170;
#line 125 "MB04NY.f"
	case 10:  goto L190;
#line 125 "MB04NY.f"
    }

/*     Code for general N. Compute */

/*     w := C*u,  C := C - tau * w * u'. */

#line 132 "MB04NY.f"
    dcopy_(m, &a[a_offset], &c__1, &dwork[1], &c__1);
#line 133 "MB04NY.f"
    dgemv_("No transpose", m, n, &c_b15, &b[b_offset], ldb, &v[1], incv, &
	    c_b15, &dwork[1], &c__1, (ftnlen)12);
#line 135 "MB04NY.f"
    d__1 = -(*tau);
#line 135 "MB04NY.f"
    daxpy_(m, &d__1, &dwork[1], &c__1, &a[a_offset], &c__1);
#line 136 "MB04NY.f"
    d__1 = -(*tau);
#line 136 "MB04NY.f"
    dger_(m, n, &d__1, &dwork[1], &c__1, &v[1], incv, &b[b_offset], ldb);
#line 137 "MB04NY.f"
    goto L210;
#line 138 "MB04NY.f"
L10:

/*     Special code for 1 x 1 Householder */

#line 142 "MB04NY.f"
    t1 = 1. - *tau;
#line 143 "MB04NY.f"
    i__1 = *m;
#line 143 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 144 "MB04NY.f"
	a[j + a_dim1] = t1 * a[j + a_dim1];
#line 145 "MB04NY.f"
/* L20: */
#line 145 "MB04NY.f"
    }
#line 146 "MB04NY.f"
    goto L210;
#line 147 "MB04NY.f"
L30:

/*     Special code for 2 x 2 Householder */

#line 151 "MB04NY.f"
    iv = 1;
#line 152 "MB04NY.f"
    if (*incv < 0) {
#line 152 "MB04NY.f"
	iv = (-(*n) + 1) * *incv + 1;
#line 152 "MB04NY.f"
    }
#line 154 "MB04NY.f"
    v1 = v[iv];
#line 155 "MB04NY.f"
    t1 = *tau * v1;
#line 156 "MB04NY.f"
    i__1 = *m;
#line 156 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 157 "MB04NY.f"
	sum = a[j + a_dim1] + v1 * b[j + b_dim1];
#line 158 "MB04NY.f"
	a[j + a_dim1] -= sum * *tau;
#line 159 "MB04NY.f"
	b[j + b_dim1] -= sum * t1;
#line 160 "MB04NY.f"
/* L40: */
#line 160 "MB04NY.f"
    }
#line 161 "MB04NY.f"
    goto L210;
#line 162 "MB04NY.f"
L50:

/*     Special code for 3 x 3 Householder */

#line 166 "MB04NY.f"
    iv = 1;
#line 167 "MB04NY.f"
    if (*incv < 0) {
#line 167 "MB04NY.f"
	iv = (-(*n) + 1) * *incv + 1;
#line 167 "MB04NY.f"
    }
#line 169 "MB04NY.f"
    v1 = v[iv];
#line 170 "MB04NY.f"
    t1 = *tau * v1;
#line 171 "MB04NY.f"
    iv += *incv;
#line 172 "MB04NY.f"
    v2 = v[iv];
#line 173 "MB04NY.f"
    t2 = *tau * v2;
#line 174 "MB04NY.f"
    i__1 = *m;
#line 174 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 175 "MB04NY.f"
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)];
#line 176 "MB04NY.f"
	a[j + a_dim1] -= sum * *tau;
#line 177 "MB04NY.f"
	b[j + b_dim1] -= sum * t1;
#line 178 "MB04NY.f"
	b[j + (b_dim1 << 1)] -= sum * t2;
#line 179 "MB04NY.f"
/* L60: */
#line 179 "MB04NY.f"
    }
#line 180 "MB04NY.f"
    goto L210;
#line 181 "MB04NY.f"
L70:

/*     Special code for 4 x 4 Householder */

#line 185 "MB04NY.f"
    iv = 1;
#line 186 "MB04NY.f"
    if (*incv < 0) {
#line 186 "MB04NY.f"
	iv = (-(*n) + 1) * *incv + 1;
#line 186 "MB04NY.f"
    }
#line 188 "MB04NY.f"
    v1 = v[iv];
#line 189 "MB04NY.f"
    t1 = *tau * v1;
#line 190 "MB04NY.f"
    iv += *incv;
#line 191 "MB04NY.f"
    v2 = v[iv];
#line 192 "MB04NY.f"
    t2 = *tau * v2;
#line 193 "MB04NY.f"
    iv += *incv;
#line 194 "MB04NY.f"
    v3 = v[iv];
#line 195 "MB04NY.f"
    t3 = *tau * v3;
#line 196 "MB04NY.f"
    i__1 = *m;
#line 196 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 197 "MB04NY.f"
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3];
#line 198 "MB04NY.f"
	a[j + a_dim1] -= sum * *tau;
#line 199 "MB04NY.f"
	b[j + b_dim1] -= sum * t1;
#line 200 "MB04NY.f"
	b[j + (b_dim1 << 1)] -= sum * t2;
#line 201 "MB04NY.f"
	b[j + b_dim1 * 3] -= sum * t3;
#line 202 "MB04NY.f"
/* L80: */
#line 202 "MB04NY.f"
    }
#line 203 "MB04NY.f"
    goto L210;
#line 204 "MB04NY.f"
L90:

/*     Special code for 5 x 5 Householder */

#line 208 "MB04NY.f"
    iv = 1;
#line 209 "MB04NY.f"
    if (*incv < 0) {
#line 209 "MB04NY.f"
	iv = (-(*n) + 1) * *incv + 1;
#line 209 "MB04NY.f"
    }
#line 211 "MB04NY.f"
    v1 = v[iv];
#line 212 "MB04NY.f"
    t1 = *tau * v1;
#line 213 "MB04NY.f"
    iv += *incv;
#line 214 "MB04NY.f"
    v2 = v[iv];
#line 215 "MB04NY.f"
    t2 = *tau * v2;
#line 216 "MB04NY.f"
    iv += *incv;
#line 217 "MB04NY.f"
    v3 = v[iv];
#line 218 "MB04NY.f"
    t3 = *tau * v3;
#line 219 "MB04NY.f"
    iv += *incv;
#line 220 "MB04NY.f"
    v4 = v[iv];
#line 221 "MB04NY.f"
    t4 = *tau * v4;
#line 222 "MB04NY.f"
    i__1 = *m;
#line 222 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 223 "MB04NY.f"
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)];
#line 225 "MB04NY.f"
	a[j + a_dim1] -= sum * *tau;
#line 226 "MB04NY.f"
	b[j + b_dim1] -= sum * t1;
#line 227 "MB04NY.f"
	b[j + (b_dim1 << 1)] -= sum * t2;
#line 228 "MB04NY.f"
	b[j + b_dim1 * 3] -= sum * t3;
#line 229 "MB04NY.f"
	b[j + (b_dim1 << 2)] -= sum * t4;
#line 230 "MB04NY.f"
/* L100: */
#line 230 "MB04NY.f"
    }
#line 231 "MB04NY.f"
    goto L210;
#line 232 "MB04NY.f"
L110:

/*     Special code for 6 x 6 Householder */

#line 236 "MB04NY.f"
    iv = 1;
#line 237 "MB04NY.f"
    if (*incv < 0) {
#line 237 "MB04NY.f"
	iv = (-(*n) + 1) * *incv + 1;
#line 237 "MB04NY.f"
    }
#line 239 "MB04NY.f"
    v1 = v[iv];
#line 240 "MB04NY.f"
    t1 = *tau * v1;
#line 241 "MB04NY.f"
    iv += *incv;
#line 242 "MB04NY.f"
    v2 = v[iv];
#line 243 "MB04NY.f"
    t2 = *tau * v2;
#line 244 "MB04NY.f"
    iv += *incv;
#line 245 "MB04NY.f"
    v3 = v[iv];
#line 246 "MB04NY.f"
    t3 = *tau * v3;
#line 247 "MB04NY.f"
    iv += *incv;
#line 248 "MB04NY.f"
    v4 = v[iv];
#line 249 "MB04NY.f"
    t4 = *tau * v4;
#line 250 "MB04NY.f"
    iv += *incv;
#line 251 "MB04NY.f"
    v5 = v[iv];
#line 252 "MB04NY.f"
    t5 = *tau * v5;
#line 253 "MB04NY.f"
    i__1 = *m;
#line 253 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 254 "MB04NY.f"
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5];
#line 256 "MB04NY.f"
	a[j + a_dim1] -= sum * *tau;
#line 257 "MB04NY.f"
	b[j + b_dim1] -= sum * t1;
#line 258 "MB04NY.f"
	b[j + (b_dim1 << 1)] -= sum * t2;
#line 259 "MB04NY.f"
	b[j + b_dim1 * 3] -= sum * t3;
#line 260 "MB04NY.f"
	b[j + (b_dim1 << 2)] -= sum * t4;
#line 261 "MB04NY.f"
	b[j + b_dim1 * 5] -= sum * t5;
#line 262 "MB04NY.f"
/* L120: */
#line 262 "MB04NY.f"
    }
#line 263 "MB04NY.f"
    goto L210;
#line 264 "MB04NY.f"
L130:

/*     Special code for 7 x 7 Householder */

#line 268 "MB04NY.f"
    iv = 1;
#line 269 "MB04NY.f"
    if (*incv < 0) {
#line 269 "MB04NY.f"
	iv = (-(*n) + 1) * *incv + 1;
#line 269 "MB04NY.f"
    }
#line 271 "MB04NY.f"
    v1 = v[iv];
#line 272 "MB04NY.f"
    t1 = *tau * v1;
#line 273 "MB04NY.f"
    iv += *incv;
#line 274 "MB04NY.f"
    v2 = v[iv];
#line 275 "MB04NY.f"
    t2 = *tau * v2;
#line 276 "MB04NY.f"
    iv += *incv;
#line 277 "MB04NY.f"
    v3 = v[iv];
#line 278 "MB04NY.f"
    t3 = *tau * v3;
#line 279 "MB04NY.f"
    iv += *incv;
#line 280 "MB04NY.f"
    v4 = v[iv];
#line 281 "MB04NY.f"
    t4 = *tau * v4;
#line 282 "MB04NY.f"
    iv += *incv;
#line 283 "MB04NY.f"
    v5 = v[iv];
#line 284 "MB04NY.f"
    t5 = *tau * v5;
#line 285 "MB04NY.f"
    iv += *incv;
#line 286 "MB04NY.f"
    v6 = v[iv];
#line 287 "MB04NY.f"
    t6 = *tau * v6;
#line 288 "MB04NY.f"
    i__1 = *m;
#line 288 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 289 "MB04NY.f"
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5] + v6 * b[j + b_dim1 * 6];
#line 291 "MB04NY.f"
	a[j + a_dim1] -= sum * *tau;
#line 292 "MB04NY.f"
	b[j + b_dim1] -= sum * t1;
#line 293 "MB04NY.f"
	b[j + (b_dim1 << 1)] -= sum * t2;
#line 294 "MB04NY.f"
	b[j + b_dim1 * 3] -= sum * t3;
#line 295 "MB04NY.f"
	b[j + (b_dim1 << 2)] -= sum * t4;
#line 296 "MB04NY.f"
	b[j + b_dim1 * 5] -= sum * t5;
#line 297 "MB04NY.f"
	b[j + b_dim1 * 6] -= sum * t6;
#line 298 "MB04NY.f"
/* L140: */
#line 298 "MB04NY.f"
    }
#line 299 "MB04NY.f"
    goto L210;
#line 300 "MB04NY.f"
L150:

/*     Special code for 8 x 8 Householder */

#line 304 "MB04NY.f"
    iv = 1;
#line 305 "MB04NY.f"
    if (*incv < 0) {
#line 305 "MB04NY.f"
	iv = (-(*n) + 1) * *incv + 1;
#line 305 "MB04NY.f"
    }
#line 307 "MB04NY.f"
    v1 = v[iv];
#line 308 "MB04NY.f"
    t1 = *tau * v1;
#line 309 "MB04NY.f"
    iv += *incv;
#line 310 "MB04NY.f"
    v2 = v[iv];
#line 311 "MB04NY.f"
    t2 = *tau * v2;
#line 312 "MB04NY.f"
    iv += *incv;
#line 313 "MB04NY.f"
    v3 = v[iv];
#line 314 "MB04NY.f"
    t3 = *tau * v3;
#line 315 "MB04NY.f"
    iv += *incv;
#line 316 "MB04NY.f"
    v4 = v[iv];
#line 317 "MB04NY.f"
    t4 = *tau * v4;
#line 318 "MB04NY.f"
    iv += *incv;
#line 319 "MB04NY.f"
    v5 = v[iv];
#line 320 "MB04NY.f"
    t5 = *tau * v5;
#line 321 "MB04NY.f"
    iv += *incv;
#line 322 "MB04NY.f"
    v6 = v[iv];
#line 323 "MB04NY.f"
    t6 = *tau * v6;
#line 324 "MB04NY.f"
    iv += *incv;
#line 325 "MB04NY.f"
    v7 = v[iv];
#line 326 "MB04NY.f"
    t7 = *tau * v7;
#line 327 "MB04NY.f"
    i__1 = *m;
#line 327 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 328 "MB04NY.f"
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5] + v6 * b[j + b_dim1 * 6] + v7 * b[j + b_dim1 *
		 7];
#line 331 "MB04NY.f"
	a[j + a_dim1] -= sum * *tau;
#line 332 "MB04NY.f"
	b[j + b_dim1] -= sum * t1;
#line 333 "MB04NY.f"
	b[j + (b_dim1 << 1)] -= sum * t2;
#line 334 "MB04NY.f"
	b[j + b_dim1 * 3] -= sum * t3;
#line 335 "MB04NY.f"
	b[j + (b_dim1 << 2)] -= sum * t4;
#line 336 "MB04NY.f"
	b[j + b_dim1 * 5] -= sum * t5;
#line 337 "MB04NY.f"
	b[j + b_dim1 * 6] -= sum * t6;
#line 338 "MB04NY.f"
	b[j + b_dim1 * 7] -= sum * t7;
#line 339 "MB04NY.f"
/* L160: */
#line 339 "MB04NY.f"
    }
#line 340 "MB04NY.f"
    goto L210;
#line 341 "MB04NY.f"
L170:

/*     Special code for 9 x 9 Householder */

#line 345 "MB04NY.f"
    iv = 1;
#line 346 "MB04NY.f"
    if (*incv < 0) {
#line 346 "MB04NY.f"
	iv = (-(*n) + 1) * *incv + 1;
#line 346 "MB04NY.f"
    }
#line 348 "MB04NY.f"
    v1 = v[iv];
#line 349 "MB04NY.f"
    t1 = *tau * v1;
#line 350 "MB04NY.f"
    iv += *incv;
#line 351 "MB04NY.f"
    v2 = v[iv];
#line 352 "MB04NY.f"
    t2 = *tau * v2;
#line 353 "MB04NY.f"
    iv += *incv;
#line 354 "MB04NY.f"
    v3 = v[iv];
#line 355 "MB04NY.f"
    t3 = *tau * v3;
#line 356 "MB04NY.f"
    iv += *incv;
#line 357 "MB04NY.f"
    v4 = v[iv];
#line 358 "MB04NY.f"
    t4 = *tau * v4;
#line 359 "MB04NY.f"
    iv += *incv;
#line 360 "MB04NY.f"
    v5 = v[iv];
#line 361 "MB04NY.f"
    t5 = *tau * v5;
#line 362 "MB04NY.f"
    iv += *incv;
#line 363 "MB04NY.f"
    v6 = v[iv];
#line 364 "MB04NY.f"
    t6 = *tau * v6;
#line 365 "MB04NY.f"
    iv += *incv;
#line 366 "MB04NY.f"
    v7 = v[iv];
#line 367 "MB04NY.f"
    t7 = *tau * v7;
#line 368 "MB04NY.f"
    iv += *incv;
#line 369 "MB04NY.f"
    v8 = v[iv];
#line 370 "MB04NY.f"
    t8 = *tau * v8;
#line 371 "MB04NY.f"
    i__1 = *m;
#line 371 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 372 "MB04NY.f"
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5] + v6 * b[j + b_dim1 * 6] + v7 * b[j + b_dim1 *
		 7] + v8 * b[j + (b_dim1 << 3)];
#line 375 "MB04NY.f"
	a[j + a_dim1] -= sum * *tau;
#line 376 "MB04NY.f"
	b[j + b_dim1] -= sum * t1;
#line 377 "MB04NY.f"
	b[j + (b_dim1 << 1)] -= sum * t2;
#line 378 "MB04NY.f"
	b[j + b_dim1 * 3] -= sum * t3;
#line 379 "MB04NY.f"
	b[j + (b_dim1 << 2)] -= sum * t4;
#line 380 "MB04NY.f"
	b[j + b_dim1 * 5] -= sum * t5;
#line 381 "MB04NY.f"
	b[j + b_dim1 * 6] -= sum * t6;
#line 382 "MB04NY.f"
	b[j + b_dim1 * 7] -= sum * t7;
#line 383 "MB04NY.f"
	b[j + (b_dim1 << 3)] -= sum * t8;
#line 384 "MB04NY.f"
/* L180: */
#line 384 "MB04NY.f"
    }
#line 385 "MB04NY.f"
    goto L210;
#line 386 "MB04NY.f"
L190:

/*     Special code for 10 x 10 Householder */

#line 390 "MB04NY.f"
    iv = 1;
#line 391 "MB04NY.f"
    if (*incv < 0) {
#line 391 "MB04NY.f"
	iv = (-(*n) + 1) * *incv + 1;
#line 391 "MB04NY.f"
    }
#line 393 "MB04NY.f"
    v1 = v[iv];
#line 394 "MB04NY.f"
    t1 = *tau * v1;
#line 395 "MB04NY.f"
    iv += *incv;
#line 396 "MB04NY.f"
    v2 = v[iv];
#line 397 "MB04NY.f"
    t2 = *tau * v2;
#line 398 "MB04NY.f"
    iv += *incv;
#line 399 "MB04NY.f"
    v3 = v[iv];
#line 400 "MB04NY.f"
    t3 = *tau * v3;
#line 401 "MB04NY.f"
    iv += *incv;
#line 402 "MB04NY.f"
    v4 = v[iv];
#line 403 "MB04NY.f"
    t4 = *tau * v4;
#line 404 "MB04NY.f"
    iv += *incv;
#line 405 "MB04NY.f"
    v5 = v[iv];
#line 406 "MB04NY.f"
    t5 = *tau * v5;
#line 407 "MB04NY.f"
    iv += *incv;
#line 408 "MB04NY.f"
    v6 = v[iv];
#line 409 "MB04NY.f"
    t6 = *tau * v6;
#line 410 "MB04NY.f"
    iv += *incv;
#line 411 "MB04NY.f"
    v7 = v[iv];
#line 412 "MB04NY.f"
    t7 = *tau * v7;
#line 413 "MB04NY.f"
    iv += *incv;
#line 414 "MB04NY.f"
    v8 = v[iv];
#line 415 "MB04NY.f"
    t8 = *tau * v8;
#line 416 "MB04NY.f"
    iv += *incv;
#line 417 "MB04NY.f"
    v9 = v[iv];
#line 418 "MB04NY.f"
    t9 = *tau * v9;
#line 419 "MB04NY.f"
    i__1 = *m;
#line 419 "MB04NY.f"
    for (j = 1; j <= i__1; ++j) {
#line 420 "MB04NY.f"
	sum = a[j + a_dim1] + v1 * b[j + b_dim1] + v2 * b[j + (b_dim1 << 1)] 
		+ v3 * b[j + b_dim1 * 3] + v4 * b[j + (b_dim1 << 2)] + v5 * b[
		j + b_dim1 * 5] + v6 * b[j + b_dim1 * 6] + v7 * b[j + b_dim1 *
		 7] + v8 * b[j + (b_dim1 << 3)] + v9 * b[j + b_dim1 * 9];
#line 423 "MB04NY.f"
	a[j + a_dim1] -= sum * *tau;
#line 424 "MB04NY.f"
	b[j + b_dim1] -= sum * t1;
#line 425 "MB04NY.f"
	b[j + (b_dim1 << 1)] -= sum * t2;
#line 426 "MB04NY.f"
	b[j + b_dim1 * 3] -= sum * t3;
#line 427 "MB04NY.f"
	b[j + (b_dim1 << 2)] -= sum * t4;
#line 428 "MB04NY.f"
	b[j + b_dim1 * 5] -= sum * t5;
#line 429 "MB04NY.f"
	b[j + b_dim1 * 6] -= sum * t6;
#line 430 "MB04NY.f"
	b[j + b_dim1 * 7] -= sum * t7;
#line 431 "MB04NY.f"
	b[j + (b_dim1 << 3)] -= sum * t8;
#line 432 "MB04NY.f"
	b[j + b_dim1 * 9] -= sum * t9;
#line 433 "MB04NY.f"
/* L200: */
#line 433 "MB04NY.f"
    }
#line 434 "MB04NY.f"
L210:
#line 435 "MB04NY.f"
    return 0;
/* *** Last line of MB04NY *** */
} /* mb04ny_ */

