#line 1 "MB04PY.f"
/* MB04PY.f -- translated by f2c (version 20100827).
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

#line 1 "MB04PY.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b15 = 1.;

/* Subroutine */ int mb04py_(char *side, integer *m, integer *n, doublereal *
	v, doublereal *tau, doublereal *c__, integer *ldc, doublereal *dwork, 
	ftnlen side_len)
{
    /* System generated locals */
    integer c_dim1, c_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, 
	    v6, v7, v8, v9, sum;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), daxpy_(integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *)
	    ;


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

/*     To apply a real elementary reflector H to a real m-by-n matrix */
/*     C, from either the left or the right. H is represented in the form */
/*                                        ( 1 ) */
/*           H = I - tau * u *u',    u  = (   ), */
/*                                        ( v ) */
/*     where tau is a real scalar and v is a real vector. */

/*     If tau = 0, then H is taken to be the unit matrix. */

/*     In-line code is used if H has order < 11. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     SIDE    CHARACTER*1 */
/*             Indicates whether the elementary reflector should be */
/*             applied from the left or from the right, as follows: */
/*             = 'L':  Compute H * C; */
/*             = 'R':  Compute C * H. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix C.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix C.  N >= 0. */

/*     V       (input) DOUBLE PRECISION array, dimension */
/*             (M-1), if SIDE = 'L', or */
/*             (N-1), if SIDE = 'R'. */
/*             The vector v in the representation of H. */

/*     TAU     (input) DOUBLE PRECISION */
/*             The scalar factor of the elementary reflector H. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix C. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the matrix H * C, if SIDE = 'L', or C * H, if SIDE = 'R'. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N), if SIDE = 'L', or */
/*                                               (M), if SIDE = 'R'. */
/*             DWORK is not referenced if H has order less than 11. */

/*     METHOD */

/*     The routine applies the elementary reflector H, taking its special */
/*     structure into account. The multiplications by the first component */
/*     of u (which is 1) are avoided, to increase the efficiency. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1999. */
/*     This is a modification of LAPACK Library routine DLARFX. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, elementary reflector, orthogonal */
/*     transformation. */

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
/*     .. Executable Statements .. */

#line 127 "MB04PY.f"
    /* Parameter adjustments */
#line 127 "MB04PY.f"
    --v;
#line 127 "MB04PY.f"
    c_dim1 = *ldc;
#line 127 "MB04PY.f"
    c_offset = 1 + c_dim1;
#line 127 "MB04PY.f"
    c__ -= c_offset;
#line 127 "MB04PY.f"
    --dwork;
#line 127 "MB04PY.f"

#line 127 "MB04PY.f"
    /* Function Body */
#line 127 "MB04PY.f"
    if (*tau == 0.) {
#line 127 "MB04PY.f"
	return 0;
#line 127 "MB04PY.f"
    }
#line 129 "MB04PY.f"
    if (lsame_(side, "L", (ftnlen)1, (ftnlen)1)) {

/*        Form  H * C, where H has order m. */

#line 133 "MB04PY.f"
	switch (*m) {
#line 133 "MB04PY.f"
	    case 1:  goto L10;
#line 133 "MB04PY.f"
	    case 2:  goto L30;
#line 133 "MB04PY.f"
	    case 3:  goto L50;
#line 133 "MB04PY.f"
	    case 4:  goto L70;
#line 133 "MB04PY.f"
	    case 5:  goto L90;
#line 133 "MB04PY.f"
	    case 6:  goto L110;
#line 133 "MB04PY.f"
	    case 7:  goto L130;
#line 133 "MB04PY.f"
	    case 8:  goto L150;
#line 133 "MB04PY.f"
	    case 9:  goto L170;
#line 133 "MB04PY.f"
	    case 10:  goto L190;
#line 133 "MB04PY.f"
	}

/*        Code for general M. */

/*        w := C'*u. */

#line 140 "MB04PY.f"
	dcopy_(n, &c__[c_offset], ldc, &dwork[1], &c__1);
#line 141 "MB04PY.f"
	i__1 = *m - 1;
#line 141 "MB04PY.f"
	dgemv_("Transpose", &i__1, n, &c_b15, &c__[c_dim1 + 2], ldc, &v[1], &
		c__1, &c_b15, &dwork[1], &c__1, (ftnlen)9);

/*        C := C - tau * u * w'. */

#line 146 "MB04PY.f"
	d__1 = -(*tau);
#line 146 "MB04PY.f"
	daxpy_(n, &d__1, &dwork[1], &c__1, &c__[c_offset], ldc);
#line 147 "MB04PY.f"
	i__1 = *m - 1;
#line 147 "MB04PY.f"
	d__1 = -(*tau);
#line 147 "MB04PY.f"
	dger_(&i__1, n, &d__1, &v[1], &c__1, &dwork[1], &c__1, &c__[c_dim1 + 
		2], ldc);
#line 148 "MB04PY.f"
	goto L410;
#line 149 "MB04PY.f"
L10:

/*        Special code for 1 x 1 Householder. */

#line 153 "MB04PY.f"
	t1 = 1. - *tau;
#line 154 "MB04PY.f"
	i__1 = *n;
#line 154 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 155 "MB04PY.f"
	    c__[j * c_dim1 + 1] = t1 * c__[j * c_dim1 + 1];
#line 156 "MB04PY.f"
/* L20: */
#line 156 "MB04PY.f"
	}
#line 157 "MB04PY.f"
	goto L410;
#line 158 "MB04PY.f"
L30:

/*        Special code for 2 x 2 Householder. */

#line 162 "MB04PY.f"
	v1 = v[1];
#line 163 "MB04PY.f"
	t1 = *tau * v1;
#line 164 "MB04PY.f"
	i__1 = *n;
#line 164 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 165 "MB04PY.f"
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2];
#line 166 "MB04PY.f"
	    c__[j * c_dim1 + 1] -= sum * *tau;
#line 167 "MB04PY.f"
	    c__[j * c_dim1 + 2] -= sum * t1;
#line 168 "MB04PY.f"
/* L40: */
#line 168 "MB04PY.f"
	}
#line 169 "MB04PY.f"
	goto L410;
#line 170 "MB04PY.f"
L50:

/*        Special code for 3 x 3 Householder. */

#line 174 "MB04PY.f"
	v1 = v[1];
#line 175 "MB04PY.f"
	t1 = *tau * v1;
#line 176 "MB04PY.f"
	v2 = v[2];
#line 177 "MB04PY.f"
	t2 = *tau * v2;
#line 178 "MB04PY.f"
	i__1 = *n;
#line 178 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 179 "MB04PY.f"
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3];
#line 180 "MB04PY.f"
	    c__[j * c_dim1 + 1] -= sum * *tau;
#line 181 "MB04PY.f"
	    c__[j * c_dim1 + 2] -= sum * t1;
#line 182 "MB04PY.f"
	    c__[j * c_dim1 + 3] -= sum * t2;
#line 183 "MB04PY.f"
/* L60: */
#line 183 "MB04PY.f"
	}
#line 184 "MB04PY.f"
	goto L410;
#line 185 "MB04PY.f"
L70:

/*        Special code for 4 x 4 Householder. */

#line 189 "MB04PY.f"
	v1 = v[1];
#line 190 "MB04PY.f"
	t1 = *tau * v1;
#line 191 "MB04PY.f"
	v2 = v[2];
#line 192 "MB04PY.f"
	t2 = *tau * v2;
#line 193 "MB04PY.f"
	v3 = v[3];
#line 194 "MB04PY.f"
	t3 = *tau * v3;
#line 195 "MB04PY.f"
	i__1 = *n;
#line 195 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 196 "MB04PY.f"
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4];
#line 198 "MB04PY.f"
	    c__[j * c_dim1 + 1] -= sum * *tau;
#line 199 "MB04PY.f"
	    c__[j * c_dim1 + 2] -= sum * t1;
#line 200 "MB04PY.f"
	    c__[j * c_dim1 + 3] -= sum * t2;
#line 201 "MB04PY.f"
	    c__[j * c_dim1 + 4] -= sum * t3;
#line 202 "MB04PY.f"
/* L80: */
#line 202 "MB04PY.f"
	}
#line 203 "MB04PY.f"
	goto L410;
#line 204 "MB04PY.f"
L90:

/*        Special code for 5 x 5 Householder. */

#line 208 "MB04PY.f"
	v1 = v[1];
#line 209 "MB04PY.f"
	t1 = *tau * v1;
#line 210 "MB04PY.f"
	v2 = v[2];
#line 211 "MB04PY.f"
	t2 = *tau * v2;
#line 212 "MB04PY.f"
	v3 = v[3];
#line 213 "MB04PY.f"
	t3 = *tau * v3;
#line 214 "MB04PY.f"
	v4 = v[4];
#line 215 "MB04PY.f"
	t4 = *tau * v4;
#line 216 "MB04PY.f"
	i__1 = *n;
#line 216 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 217 "MB04PY.f"
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5];
#line 219 "MB04PY.f"
	    c__[j * c_dim1 + 1] -= sum * *tau;
#line 220 "MB04PY.f"
	    c__[j * c_dim1 + 2] -= sum * t1;
#line 221 "MB04PY.f"
	    c__[j * c_dim1 + 3] -= sum * t2;
#line 222 "MB04PY.f"
	    c__[j * c_dim1 + 4] -= sum * t3;
#line 223 "MB04PY.f"
	    c__[j * c_dim1 + 5] -= sum * t4;
#line 224 "MB04PY.f"
/* L100: */
#line 224 "MB04PY.f"
	}
#line 225 "MB04PY.f"
	goto L410;
#line 226 "MB04PY.f"
L110:

/*        Special code for 6 x 6 Householder. */

#line 230 "MB04PY.f"
	v1 = v[1];
#line 231 "MB04PY.f"
	t1 = *tau * v1;
#line 232 "MB04PY.f"
	v2 = v[2];
#line 233 "MB04PY.f"
	t2 = *tau * v2;
#line 234 "MB04PY.f"
	v3 = v[3];
#line 235 "MB04PY.f"
	t3 = *tau * v3;
#line 236 "MB04PY.f"
	v4 = v[4];
#line 237 "MB04PY.f"
	t4 = *tau * v4;
#line 238 "MB04PY.f"
	v5 = v[5];
#line 239 "MB04PY.f"
	t5 = *tau * v5;
#line 240 "MB04PY.f"
	i__1 = *n;
#line 240 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 241 "MB04PY.f"
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6];
#line 243 "MB04PY.f"
	    c__[j * c_dim1 + 1] -= sum * *tau;
#line 244 "MB04PY.f"
	    c__[j * c_dim1 + 2] -= sum * t1;
#line 245 "MB04PY.f"
	    c__[j * c_dim1 + 3] -= sum * t2;
#line 246 "MB04PY.f"
	    c__[j * c_dim1 + 4] -= sum * t3;
#line 247 "MB04PY.f"
	    c__[j * c_dim1 + 5] -= sum * t4;
#line 248 "MB04PY.f"
	    c__[j * c_dim1 + 6] -= sum * t5;
#line 249 "MB04PY.f"
/* L120: */
#line 249 "MB04PY.f"
	}
#line 250 "MB04PY.f"
	goto L410;
#line 251 "MB04PY.f"
L130:

/*        Special code for 7 x 7 Householder. */

#line 255 "MB04PY.f"
	v1 = v[1];
#line 256 "MB04PY.f"
	t1 = *tau * v1;
#line 257 "MB04PY.f"
	v2 = v[2];
#line 258 "MB04PY.f"
	t2 = *tau * v2;
#line 259 "MB04PY.f"
	v3 = v[3];
#line 260 "MB04PY.f"
	t3 = *tau * v3;
#line 261 "MB04PY.f"
	v4 = v[4];
#line 262 "MB04PY.f"
	t4 = *tau * v4;
#line 263 "MB04PY.f"
	v5 = v[5];
#line 264 "MB04PY.f"
	t5 = *tau * v5;
#line 265 "MB04PY.f"
	v6 = v[6];
#line 266 "MB04PY.f"
	t6 = *tau * v6;
#line 267 "MB04PY.f"
	i__1 = *n;
#line 267 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 268 "MB04PY.f"
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6] + v6 * c__[j * 
		    c_dim1 + 7];
#line 271 "MB04PY.f"
	    c__[j * c_dim1 + 1] -= sum * *tau;
#line 272 "MB04PY.f"
	    c__[j * c_dim1 + 2] -= sum * t1;
#line 273 "MB04PY.f"
	    c__[j * c_dim1 + 3] -= sum * t2;
#line 274 "MB04PY.f"
	    c__[j * c_dim1 + 4] -= sum * t3;
#line 275 "MB04PY.f"
	    c__[j * c_dim1 + 5] -= sum * t4;
#line 276 "MB04PY.f"
	    c__[j * c_dim1 + 6] -= sum * t5;
#line 277 "MB04PY.f"
	    c__[j * c_dim1 + 7] -= sum * t6;
#line 278 "MB04PY.f"
/* L140: */
#line 278 "MB04PY.f"
	}
#line 279 "MB04PY.f"
	goto L410;
#line 280 "MB04PY.f"
L150:

/*        Special code for 8 x 8 Householder. */

#line 284 "MB04PY.f"
	v1 = v[1];
#line 285 "MB04PY.f"
	t1 = *tau * v1;
#line 286 "MB04PY.f"
	v2 = v[2];
#line 287 "MB04PY.f"
	t2 = *tau * v2;
#line 288 "MB04PY.f"
	v3 = v[3];
#line 289 "MB04PY.f"
	t3 = *tau * v3;
#line 290 "MB04PY.f"
	v4 = v[4];
#line 291 "MB04PY.f"
	t4 = *tau * v4;
#line 292 "MB04PY.f"
	v5 = v[5];
#line 293 "MB04PY.f"
	t5 = *tau * v5;
#line 294 "MB04PY.f"
	v6 = v[6];
#line 295 "MB04PY.f"
	t6 = *tau * v6;
#line 296 "MB04PY.f"
	v7 = v[7];
#line 297 "MB04PY.f"
	t7 = *tau * v7;
#line 298 "MB04PY.f"
	i__1 = *n;
#line 298 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 299 "MB04PY.f"
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6] + v6 * c__[j * 
		    c_dim1 + 7] + v7 * c__[j * c_dim1 + 8];
#line 302 "MB04PY.f"
	    c__[j * c_dim1 + 1] -= sum * *tau;
#line 303 "MB04PY.f"
	    c__[j * c_dim1 + 2] -= sum * t1;
#line 304 "MB04PY.f"
	    c__[j * c_dim1 + 3] -= sum * t2;
#line 305 "MB04PY.f"
	    c__[j * c_dim1 + 4] -= sum * t3;
#line 306 "MB04PY.f"
	    c__[j * c_dim1 + 5] -= sum * t4;
#line 307 "MB04PY.f"
	    c__[j * c_dim1 + 6] -= sum * t5;
#line 308 "MB04PY.f"
	    c__[j * c_dim1 + 7] -= sum * t6;
#line 309 "MB04PY.f"
	    c__[j * c_dim1 + 8] -= sum * t7;
#line 310 "MB04PY.f"
/* L160: */
#line 310 "MB04PY.f"
	}
#line 311 "MB04PY.f"
	goto L410;
#line 312 "MB04PY.f"
L170:

/*        Special code for 9 x 9 Householder. */

#line 316 "MB04PY.f"
	v1 = v[1];
#line 317 "MB04PY.f"
	t1 = *tau * v1;
#line 318 "MB04PY.f"
	v2 = v[2];
#line 319 "MB04PY.f"
	t2 = *tau * v2;
#line 320 "MB04PY.f"
	v3 = v[3];
#line 321 "MB04PY.f"
	t3 = *tau * v3;
#line 322 "MB04PY.f"
	v4 = v[4];
#line 323 "MB04PY.f"
	t4 = *tau * v4;
#line 324 "MB04PY.f"
	v5 = v[5];
#line 325 "MB04PY.f"
	t5 = *tau * v5;
#line 326 "MB04PY.f"
	v6 = v[6];
#line 327 "MB04PY.f"
	t6 = *tau * v6;
#line 328 "MB04PY.f"
	v7 = v[7];
#line 329 "MB04PY.f"
	t7 = *tau * v7;
#line 330 "MB04PY.f"
	v8 = v[8];
#line 331 "MB04PY.f"
	t8 = *tau * v8;
#line 332 "MB04PY.f"
	i__1 = *n;
#line 332 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 333 "MB04PY.f"
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6] + v6 * c__[j * 
		    c_dim1 + 7] + v7 * c__[j * c_dim1 + 8] + v8 * c__[j * 
		    c_dim1 + 9];
#line 336 "MB04PY.f"
	    c__[j * c_dim1 + 1] -= sum * *tau;
#line 337 "MB04PY.f"
	    c__[j * c_dim1 + 2] -= sum * t1;
#line 338 "MB04PY.f"
	    c__[j * c_dim1 + 3] -= sum * t2;
#line 339 "MB04PY.f"
	    c__[j * c_dim1 + 4] -= sum * t3;
#line 340 "MB04PY.f"
	    c__[j * c_dim1 + 5] -= sum * t4;
#line 341 "MB04PY.f"
	    c__[j * c_dim1 + 6] -= sum * t5;
#line 342 "MB04PY.f"
	    c__[j * c_dim1 + 7] -= sum * t6;
#line 343 "MB04PY.f"
	    c__[j * c_dim1 + 8] -= sum * t7;
#line 344 "MB04PY.f"
	    c__[j * c_dim1 + 9] -= sum * t8;
#line 345 "MB04PY.f"
/* L180: */
#line 345 "MB04PY.f"
	}
#line 346 "MB04PY.f"
	goto L410;
#line 347 "MB04PY.f"
L190:

/*        Special code for 10 x 10 Householder. */

#line 351 "MB04PY.f"
	v1 = v[1];
#line 352 "MB04PY.f"
	t1 = *tau * v1;
#line 353 "MB04PY.f"
	v2 = v[2];
#line 354 "MB04PY.f"
	t2 = *tau * v2;
#line 355 "MB04PY.f"
	v3 = v[3];
#line 356 "MB04PY.f"
	t3 = *tau * v3;
#line 357 "MB04PY.f"
	v4 = v[4];
#line 358 "MB04PY.f"
	t4 = *tau * v4;
#line 359 "MB04PY.f"
	v5 = v[5];
#line 360 "MB04PY.f"
	t5 = *tau * v5;
#line 361 "MB04PY.f"
	v6 = v[6];
#line 362 "MB04PY.f"
	t6 = *tau * v6;
#line 363 "MB04PY.f"
	v7 = v[7];
#line 364 "MB04PY.f"
	t7 = *tau * v7;
#line 365 "MB04PY.f"
	v8 = v[8];
#line 366 "MB04PY.f"
	t8 = *tau * v8;
#line 367 "MB04PY.f"
	v9 = v[9];
#line 368 "MB04PY.f"
	t9 = *tau * v9;
#line 369 "MB04PY.f"
	i__1 = *n;
#line 369 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 370 "MB04PY.f"
	    sum = c__[j * c_dim1 + 1] + v1 * c__[j * c_dim1 + 2] + v2 * c__[j 
		    * c_dim1 + 3] + v3 * c__[j * c_dim1 + 4] + v4 * c__[j * 
		    c_dim1 + 5] + v5 * c__[j * c_dim1 + 6] + v6 * c__[j * 
		    c_dim1 + 7] + v7 * c__[j * c_dim1 + 8] + v8 * c__[j * 
		    c_dim1 + 9] + v9 * c__[j * c_dim1 + 10];
#line 374 "MB04PY.f"
	    c__[j * c_dim1 + 1] -= sum * *tau;
#line 375 "MB04PY.f"
	    c__[j * c_dim1 + 2] -= sum * t1;
#line 376 "MB04PY.f"
	    c__[j * c_dim1 + 3] -= sum * t2;
#line 377 "MB04PY.f"
	    c__[j * c_dim1 + 4] -= sum * t3;
#line 378 "MB04PY.f"
	    c__[j * c_dim1 + 5] -= sum * t4;
#line 379 "MB04PY.f"
	    c__[j * c_dim1 + 6] -= sum * t5;
#line 380 "MB04PY.f"
	    c__[j * c_dim1 + 7] -= sum * t6;
#line 381 "MB04PY.f"
	    c__[j * c_dim1 + 8] -= sum * t7;
#line 382 "MB04PY.f"
	    c__[j * c_dim1 + 9] -= sum * t8;
#line 383 "MB04PY.f"
	    c__[j * c_dim1 + 10] -= sum * t9;
#line 384 "MB04PY.f"
/* L200: */
#line 384 "MB04PY.f"
	}
#line 385 "MB04PY.f"
	goto L410;
#line 386 "MB04PY.f"
    } else {

/*        Form  C * H, where H has order n. */

#line 390 "MB04PY.f"
	switch (*n) {
#line 390 "MB04PY.f"
	    case 1:  goto L210;
#line 390 "MB04PY.f"
	    case 2:  goto L230;
#line 390 "MB04PY.f"
	    case 3:  goto L250;
#line 390 "MB04PY.f"
	    case 4:  goto L270;
#line 390 "MB04PY.f"
	    case 5:  goto L290;
#line 390 "MB04PY.f"
	    case 6:  goto L310;
#line 390 "MB04PY.f"
	    case 7:  goto L330;
#line 390 "MB04PY.f"
	    case 8:  goto L350;
#line 390 "MB04PY.f"
	    case 9:  goto L370;
#line 390 "MB04PY.f"
	    case 10:  goto L390;
#line 390 "MB04PY.f"
	}

/*        Code for general N. */

/*        w := C * u. */

#line 397 "MB04PY.f"
	dcopy_(m, &c__[c_offset], &c__1, &dwork[1], &c__1);
#line 398 "MB04PY.f"
	i__1 = *n - 1;
#line 398 "MB04PY.f"
	dgemv_("No transpose", m, &i__1, &c_b15, &c__[(c_dim1 << 1) + 1], ldc,
		 &v[1], &c__1, &c_b15, &dwork[1], &c__1, (ftnlen)12);

/*        C := C - tau * w * u'. */

#line 403 "MB04PY.f"
	d__1 = -(*tau);
#line 403 "MB04PY.f"
	daxpy_(m, &d__1, &dwork[1], &c__1, &c__[c_offset], &c__1);
#line 404 "MB04PY.f"
	i__1 = *n - 1;
#line 404 "MB04PY.f"
	d__1 = -(*tau);
#line 404 "MB04PY.f"
	dger_(m, &i__1, &d__1, &dwork[1], &c__1, &v[1], &c__1, &c__[(c_dim1 <<
		 1) + 1], ldc);
#line 405 "MB04PY.f"
	goto L410;
#line 406 "MB04PY.f"
L210:

/*        Special code for 1 x 1 Householder. */

#line 410 "MB04PY.f"
	t1 = 1. - *tau;
#line 411 "MB04PY.f"
	i__1 = *m;
#line 411 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 412 "MB04PY.f"
	    c__[j + c_dim1] = t1 * c__[j + c_dim1];
#line 413 "MB04PY.f"
/* L220: */
#line 413 "MB04PY.f"
	}
#line 414 "MB04PY.f"
	goto L410;
#line 415 "MB04PY.f"
L230:

/*        Special code for 2 x 2 Householder. */

#line 419 "MB04PY.f"
	v1 = v[1];
#line 420 "MB04PY.f"
	t1 = *tau * v1;
#line 421 "MB04PY.f"
	i__1 = *m;
#line 421 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 422 "MB04PY.f"
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)];
#line 423 "MB04PY.f"
	    c__[j + c_dim1] -= sum * *tau;
#line 424 "MB04PY.f"
	    c__[j + (c_dim1 << 1)] -= sum * t1;
#line 425 "MB04PY.f"
/* L240: */
#line 425 "MB04PY.f"
	}
#line 426 "MB04PY.f"
	goto L410;
#line 427 "MB04PY.f"
L250:

/*        Special code for 3 x 3 Householder. */

#line 431 "MB04PY.f"
	v1 = v[1];
#line 432 "MB04PY.f"
	t1 = *tau * v1;
#line 433 "MB04PY.f"
	v2 = v[2];
#line 434 "MB04PY.f"
	t2 = *tau * v2;
#line 435 "MB04PY.f"
	i__1 = *m;
#line 435 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 436 "MB04PY.f"
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3];
#line 437 "MB04PY.f"
	    c__[j + c_dim1] -= sum * *tau;
#line 438 "MB04PY.f"
	    c__[j + (c_dim1 << 1)] -= sum * t1;
#line 439 "MB04PY.f"
	    c__[j + c_dim1 * 3] -= sum * t2;
#line 440 "MB04PY.f"
/* L260: */
#line 440 "MB04PY.f"
	}
#line 441 "MB04PY.f"
	goto L410;
#line 442 "MB04PY.f"
L270:

/*        Special code for 4 x 4 Householder. */

#line 446 "MB04PY.f"
	v1 = v[1];
#line 447 "MB04PY.f"
	t1 = *tau * v1;
#line 448 "MB04PY.f"
	v2 = v[2];
#line 449 "MB04PY.f"
	t2 = *tau * v2;
#line 450 "MB04PY.f"
	v3 = v[3];
#line 451 "MB04PY.f"
	t3 = *tau * v3;
#line 452 "MB04PY.f"
	i__1 = *m;
#line 452 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 453 "MB04PY.f"
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)];
#line 455 "MB04PY.f"
	    c__[j + c_dim1] -= sum * *tau;
#line 456 "MB04PY.f"
	    c__[j + (c_dim1 << 1)] -= sum * t1;
#line 457 "MB04PY.f"
	    c__[j + c_dim1 * 3] -= sum * t2;
#line 458 "MB04PY.f"
	    c__[j + (c_dim1 << 2)] -= sum * t3;
#line 459 "MB04PY.f"
/* L280: */
#line 459 "MB04PY.f"
	}
#line 460 "MB04PY.f"
	goto L410;
#line 461 "MB04PY.f"
L290:

/*        Special code for 5 x 5 Householder. */

#line 465 "MB04PY.f"
	v1 = v[1];
#line 466 "MB04PY.f"
	t1 = *tau * v1;
#line 467 "MB04PY.f"
	v2 = v[2];
#line 468 "MB04PY.f"
	t2 = *tau * v2;
#line 469 "MB04PY.f"
	v3 = v[3];
#line 470 "MB04PY.f"
	t3 = *tau * v3;
#line 471 "MB04PY.f"
	v4 = v[4];
#line 472 "MB04PY.f"
	t4 = *tau * v4;
#line 473 "MB04PY.f"
	i__1 = *m;
#line 473 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 474 "MB04PY.f"
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5];
#line 476 "MB04PY.f"
	    c__[j + c_dim1] -= sum * *tau;
#line 477 "MB04PY.f"
	    c__[j + (c_dim1 << 1)] -= sum * t1;
#line 478 "MB04PY.f"
	    c__[j + c_dim1 * 3] -= sum * t2;
#line 479 "MB04PY.f"
	    c__[j + (c_dim1 << 2)] -= sum * t3;
#line 480 "MB04PY.f"
	    c__[j + c_dim1 * 5] -= sum * t4;
#line 481 "MB04PY.f"
/* L300: */
#line 481 "MB04PY.f"
	}
#line 482 "MB04PY.f"
	goto L410;
#line 483 "MB04PY.f"
L310:

/*        Special code for 6 x 6 Householder. */

#line 487 "MB04PY.f"
	v1 = v[1];
#line 488 "MB04PY.f"
	t1 = *tau * v1;
#line 489 "MB04PY.f"
	v2 = v[2];
#line 490 "MB04PY.f"
	t2 = *tau * v2;
#line 491 "MB04PY.f"
	v3 = v[3];
#line 492 "MB04PY.f"
	t3 = *tau * v3;
#line 493 "MB04PY.f"
	v4 = v[4];
#line 494 "MB04PY.f"
	t4 = *tau * v4;
#line 495 "MB04PY.f"
	v5 = v[5];
#line 496 "MB04PY.f"
	t5 = *tau * v5;
#line 497 "MB04PY.f"
	i__1 = *m;
#line 497 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 498 "MB04PY.f"
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6];
#line 500 "MB04PY.f"
	    c__[j + c_dim1] -= sum * *tau;
#line 501 "MB04PY.f"
	    c__[j + (c_dim1 << 1)] -= sum * t1;
#line 502 "MB04PY.f"
	    c__[j + c_dim1 * 3] -= sum * t2;
#line 503 "MB04PY.f"
	    c__[j + (c_dim1 << 2)] -= sum * t3;
#line 504 "MB04PY.f"
	    c__[j + c_dim1 * 5] -= sum * t4;
#line 505 "MB04PY.f"
	    c__[j + c_dim1 * 6] -= sum * t5;
#line 506 "MB04PY.f"
/* L320: */
#line 506 "MB04PY.f"
	}
#line 507 "MB04PY.f"
	goto L410;
#line 508 "MB04PY.f"
L330:

/*        Special code for 7 x 7 Householder. */

#line 512 "MB04PY.f"
	v1 = v[1];
#line 513 "MB04PY.f"
	t1 = *tau * v1;
#line 514 "MB04PY.f"
	v2 = v[2];
#line 515 "MB04PY.f"
	t2 = *tau * v2;
#line 516 "MB04PY.f"
	v3 = v[3];
#line 517 "MB04PY.f"
	t3 = *tau * v3;
#line 518 "MB04PY.f"
	v4 = v[4];
#line 519 "MB04PY.f"
	t4 = *tau * v4;
#line 520 "MB04PY.f"
	v5 = v[5];
#line 521 "MB04PY.f"
	t5 = *tau * v5;
#line 522 "MB04PY.f"
	v6 = v[6];
#line 523 "MB04PY.f"
	t6 = *tau * v6;
#line 524 "MB04PY.f"
	i__1 = *m;
#line 524 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 525 "MB04PY.f"
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6] + v6 * c__[j + 
		    c_dim1 * 7];
#line 528 "MB04PY.f"
	    c__[j + c_dim1] -= sum * *tau;
#line 529 "MB04PY.f"
	    c__[j + (c_dim1 << 1)] -= sum * t1;
#line 530 "MB04PY.f"
	    c__[j + c_dim1 * 3] -= sum * t2;
#line 531 "MB04PY.f"
	    c__[j + (c_dim1 << 2)] -= sum * t3;
#line 532 "MB04PY.f"
	    c__[j + c_dim1 * 5] -= sum * t4;
#line 533 "MB04PY.f"
	    c__[j + c_dim1 * 6] -= sum * t5;
#line 534 "MB04PY.f"
	    c__[j + c_dim1 * 7] -= sum * t6;
#line 535 "MB04PY.f"
/* L340: */
#line 535 "MB04PY.f"
	}
#line 536 "MB04PY.f"
	goto L410;
#line 537 "MB04PY.f"
L350:

/*        Special code for 8 x 8 Householder. */

#line 541 "MB04PY.f"
	v1 = v[1];
#line 542 "MB04PY.f"
	t1 = *tau * v1;
#line 543 "MB04PY.f"
	v2 = v[2];
#line 544 "MB04PY.f"
	t2 = *tau * v2;
#line 545 "MB04PY.f"
	v3 = v[3];
#line 546 "MB04PY.f"
	t3 = *tau * v3;
#line 547 "MB04PY.f"
	v4 = v[4];
#line 548 "MB04PY.f"
	t4 = *tau * v4;
#line 549 "MB04PY.f"
	v5 = v[5];
#line 550 "MB04PY.f"
	t5 = *tau * v5;
#line 551 "MB04PY.f"
	v6 = v[6];
#line 552 "MB04PY.f"
	t6 = *tau * v6;
#line 553 "MB04PY.f"
	v7 = v[7];
#line 554 "MB04PY.f"
	t7 = *tau * v7;
#line 555 "MB04PY.f"
	i__1 = *m;
#line 555 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 556 "MB04PY.f"
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6] + v6 * c__[j + 
		    c_dim1 * 7] + v7 * c__[j + (c_dim1 << 3)];
#line 559 "MB04PY.f"
	    c__[j + c_dim1] -= sum * *tau;
#line 560 "MB04PY.f"
	    c__[j + (c_dim1 << 1)] -= sum * t1;
#line 561 "MB04PY.f"
	    c__[j + c_dim1 * 3] -= sum * t2;
#line 562 "MB04PY.f"
	    c__[j + (c_dim1 << 2)] -= sum * t3;
#line 563 "MB04PY.f"
	    c__[j + c_dim1 * 5] -= sum * t4;
#line 564 "MB04PY.f"
	    c__[j + c_dim1 * 6] -= sum * t5;
#line 565 "MB04PY.f"
	    c__[j + c_dim1 * 7] -= sum * t6;
#line 566 "MB04PY.f"
	    c__[j + (c_dim1 << 3)] -= sum * t7;
#line 567 "MB04PY.f"
/* L360: */
#line 567 "MB04PY.f"
	}
#line 568 "MB04PY.f"
	goto L410;
#line 569 "MB04PY.f"
L370:

/*        Special code for 9 x 9 Householder. */

#line 573 "MB04PY.f"
	v1 = v[1];
#line 574 "MB04PY.f"
	t1 = *tau * v1;
#line 575 "MB04PY.f"
	v2 = v[2];
#line 576 "MB04PY.f"
	t2 = *tau * v2;
#line 577 "MB04PY.f"
	v3 = v[3];
#line 578 "MB04PY.f"
	t3 = *tau * v3;
#line 579 "MB04PY.f"
	v4 = v[4];
#line 580 "MB04PY.f"
	t4 = *tau * v4;
#line 581 "MB04PY.f"
	v5 = v[5];
#line 582 "MB04PY.f"
	t5 = *tau * v5;
#line 583 "MB04PY.f"
	v6 = v[6];
#line 584 "MB04PY.f"
	t6 = *tau * v6;
#line 585 "MB04PY.f"
	v7 = v[7];
#line 586 "MB04PY.f"
	t7 = *tau * v7;
#line 587 "MB04PY.f"
	v8 = v[8];
#line 588 "MB04PY.f"
	t8 = *tau * v8;
#line 589 "MB04PY.f"
	i__1 = *m;
#line 589 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 590 "MB04PY.f"
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6] + v6 * c__[j + 
		    c_dim1 * 7] + v7 * c__[j + (c_dim1 << 3)] + v8 * c__[j + 
		    c_dim1 * 9];
#line 593 "MB04PY.f"
	    c__[j + c_dim1] -= sum * *tau;
#line 594 "MB04PY.f"
	    c__[j + (c_dim1 << 1)] -= sum * t1;
#line 595 "MB04PY.f"
	    c__[j + c_dim1 * 3] -= sum * t2;
#line 596 "MB04PY.f"
	    c__[j + (c_dim1 << 2)] -= sum * t3;
#line 597 "MB04PY.f"
	    c__[j + c_dim1 * 5] -= sum * t4;
#line 598 "MB04PY.f"
	    c__[j + c_dim1 * 6] -= sum * t5;
#line 599 "MB04PY.f"
	    c__[j + c_dim1 * 7] -= sum * t6;
#line 600 "MB04PY.f"
	    c__[j + (c_dim1 << 3)] -= sum * t7;
#line 601 "MB04PY.f"
	    c__[j + c_dim1 * 9] -= sum * t8;
#line 602 "MB04PY.f"
/* L380: */
#line 602 "MB04PY.f"
	}
#line 603 "MB04PY.f"
	goto L410;
#line 604 "MB04PY.f"
L390:

/*        Special code for 10 x 10 Householder. */

#line 608 "MB04PY.f"
	v1 = v[1];
#line 609 "MB04PY.f"
	t1 = *tau * v1;
#line 610 "MB04PY.f"
	v2 = v[2];
#line 611 "MB04PY.f"
	t2 = *tau * v2;
#line 612 "MB04PY.f"
	v3 = v[3];
#line 613 "MB04PY.f"
	t3 = *tau * v3;
#line 614 "MB04PY.f"
	v4 = v[4];
#line 615 "MB04PY.f"
	t4 = *tau * v4;
#line 616 "MB04PY.f"
	v5 = v[5];
#line 617 "MB04PY.f"
	t5 = *tau * v5;
#line 618 "MB04PY.f"
	v6 = v[6];
#line 619 "MB04PY.f"
	t6 = *tau * v6;
#line 620 "MB04PY.f"
	v7 = v[7];
#line 621 "MB04PY.f"
	t7 = *tau * v7;
#line 622 "MB04PY.f"
	v8 = v[8];
#line 623 "MB04PY.f"
	t8 = *tau * v8;
#line 624 "MB04PY.f"
	v9 = v[9];
#line 625 "MB04PY.f"
	t9 = *tau * v9;
#line 626 "MB04PY.f"
	i__1 = *m;
#line 626 "MB04PY.f"
	for (j = 1; j <= i__1; ++j) {
#line 627 "MB04PY.f"
	    sum = c__[j + c_dim1] + v1 * c__[j + (c_dim1 << 1)] + v2 * c__[j 
		    + c_dim1 * 3] + v3 * c__[j + (c_dim1 << 2)] + v4 * c__[j 
		    + c_dim1 * 5] + v5 * c__[j + c_dim1 * 6] + v6 * c__[j + 
		    c_dim1 * 7] + v7 * c__[j + (c_dim1 << 3)] + v8 * c__[j + 
		    c_dim1 * 9] + v9 * c__[j + c_dim1 * 10];
#line 631 "MB04PY.f"
	    c__[j + c_dim1] -= sum * *tau;
#line 632 "MB04PY.f"
	    c__[j + (c_dim1 << 1)] -= sum * t1;
#line 633 "MB04PY.f"
	    c__[j + c_dim1 * 3] -= sum * t2;
#line 634 "MB04PY.f"
	    c__[j + (c_dim1 << 2)] -= sum * t3;
#line 635 "MB04PY.f"
	    c__[j + c_dim1 * 5] -= sum * t4;
#line 636 "MB04PY.f"
	    c__[j + c_dim1 * 6] -= sum * t5;
#line 637 "MB04PY.f"
	    c__[j + c_dim1 * 7] -= sum * t6;
#line 638 "MB04PY.f"
	    c__[j + (c_dim1 << 3)] -= sum * t7;
#line 639 "MB04PY.f"
	    c__[j + c_dim1 * 9] -= sum * t8;
#line 640 "MB04PY.f"
	    c__[j + c_dim1 * 10] -= sum * t9;
#line 641 "MB04PY.f"
/* L400: */
#line 641 "MB04PY.f"
	}
#line 642 "MB04PY.f"
	goto L410;
#line 643 "MB04PY.f"
    }
#line 644 "MB04PY.f"
L410:
#line 645 "MB04PY.f"
    return 0;

/* *** Last line of MB04PY *** */
} /* mb04py_ */

