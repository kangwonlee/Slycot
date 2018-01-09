#line 1 "MB04OY.f"
/* MB04OY.f -- translated by f2c (version 20100827).
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

#line 1 "MB04OY.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b14 = 1.;

/* Subroutine */ int mb04oy_(integer *m, integer *n, doublereal *v, 
	doublereal *tau, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *dwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal t1, t2, t3, t4, t5, t6, t7, t8, t9, v1, v2, v3, v4, v5, 
	    v6, v7, v8, v9, sum;
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

/*     To apply a real elementary reflector H to a real (m+1)-by-n */
/*     matrix C = [ A ], from the left, where A has one row. H is */
/*                [ B ] */
/*     represented in the form */
/*                                        ( 1 ) */
/*           H = I - tau * u *u',    u  = (   ), */
/*                                        ( v ) */
/*     where tau is a real scalar and v is a real m-vector. */

/*     If tau = 0, then H is taken to be the unit matrix. */

/*     In-line code is used if H has order < 11. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix B.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrices A and B.  N >= 0. */

/*     V       (input) DOUBLE PRECISION array, dimension (M) */
/*             The vector v in the representation of H. */

/*     TAU     (input) DOUBLE PRECISION */
/*             The scalar factor of the elementary reflector H. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading 1-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading 1-by-N part of this array contains */
/*             the updated matrix A (the first row of H * C). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= 1. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the updated matrix B (the last m rows of H * C). */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,M). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */
/*             DWORK is not referenced if H has order less than 11. */

/*     METHOD */

/*     The routine applies the elementary reflector H, taking the special */
/*     structure of C into account. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */
/*     Based on LAPACK routines DLARFX and DLATZM. */

/*     REVISIONS */

/*     Dec. 1997. */

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

#line 117 "MB04OY.f"
    /* Parameter adjustments */
#line 117 "MB04OY.f"
    --v;
#line 117 "MB04OY.f"
    a_dim1 = *lda;
#line 117 "MB04OY.f"
    a_offset = 1 + a_dim1;
#line 117 "MB04OY.f"
    a -= a_offset;
#line 117 "MB04OY.f"
    b_dim1 = *ldb;
#line 117 "MB04OY.f"
    b_offset = 1 + b_dim1;
#line 117 "MB04OY.f"
    b -= b_offset;
#line 117 "MB04OY.f"
    --dwork;
#line 117 "MB04OY.f"

#line 117 "MB04OY.f"
    /* Function Body */
#line 117 "MB04OY.f"
    if (*tau == 0.) {
#line 117 "MB04OY.f"
	return 0;
#line 117 "MB04OY.f"
    }

/*     Form  H * C, where H has order m+1. */

#line 122 "MB04OY.f"
    switch (*m + 1) {
#line 122 "MB04OY.f"
	case 1:  goto L10;
#line 122 "MB04OY.f"
	case 2:  goto L30;
#line 122 "MB04OY.f"
	case 3:  goto L50;
#line 122 "MB04OY.f"
	case 4:  goto L70;
#line 122 "MB04OY.f"
	case 5:  goto L90;
#line 122 "MB04OY.f"
	case 6:  goto L110;
#line 122 "MB04OY.f"
	case 7:  goto L130;
#line 122 "MB04OY.f"
	case 8:  goto L150;
#line 122 "MB04OY.f"
	case 9:  goto L170;
#line 122 "MB04OY.f"
	case 10:  goto L190;
#line 122 "MB04OY.f"
    }

/*     Code for general M. Compute */

/*     w := C'*u,  C := C - tau * u * w'. */

#line 129 "MB04OY.f"
    dcopy_(n, &a[a_offset], lda, &dwork[1], &c__1);
#line 130 "MB04OY.f"
    dgemv_("Transpose", m, n, &c_b14, &b[b_offset], ldb, &v[1], &c__1, &c_b14,
	     &dwork[1], &c__1, (ftnlen)9);
#line 131 "MB04OY.f"
    d__1 = -(*tau);
#line 131 "MB04OY.f"
    daxpy_(n, &d__1, &dwork[1], &c__1, &a[a_offset], lda);
#line 132 "MB04OY.f"
    d__1 = -(*tau);
#line 132 "MB04OY.f"
    dger_(m, n, &d__1, &v[1], &c__1, &dwork[1], &c__1, &b[b_offset], ldb);
#line 133 "MB04OY.f"
    goto L210;
#line 134 "MB04OY.f"
L10:

/*     Special code for 1 x 1 Householder */

#line 138 "MB04OY.f"
    t1 = 1. - *tau;
#line 139 "MB04OY.f"
    i__1 = *n;
#line 139 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 140 "MB04OY.f"
	a[j * a_dim1 + 1] = t1 * a[j * a_dim1 + 1];
#line 141 "MB04OY.f"
/* L20: */
#line 141 "MB04OY.f"
    }
#line 142 "MB04OY.f"
    goto L210;
#line 143 "MB04OY.f"
L30:

/*     Special code for 2 x 2 Householder */

#line 147 "MB04OY.f"
    v1 = v[1];
#line 148 "MB04OY.f"
    t1 = *tau * v1;
#line 149 "MB04OY.f"
    i__1 = *n;
#line 149 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 150 "MB04OY.f"
	sum = a[j * a_dim1 + 1] + v1 * b[j * b_dim1 + 1];
#line 151 "MB04OY.f"
	a[j * a_dim1 + 1] -= sum * *tau;
#line 152 "MB04OY.f"
	b[j * b_dim1 + 1] -= sum * t1;
#line 153 "MB04OY.f"
/* L40: */
#line 153 "MB04OY.f"
    }
#line 154 "MB04OY.f"
    goto L210;
#line 155 "MB04OY.f"
L50:

/*     Special code for 3 x 3 Householder */

#line 159 "MB04OY.f"
    v1 = v[1];
#line 160 "MB04OY.f"
    t1 = *tau * v1;
#line 161 "MB04OY.f"
    v2 = v[2];
#line 162 "MB04OY.f"
    t2 = *tau * v2;
#line 163 "MB04OY.f"
    i__1 = *n;
#line 163 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 164 "MB04OY.f"
	sum = a[j * a_dim1 + 1] + v1 * b[j * b_dim1 + 1] + v2 * b[j * b_dim1 
		+ 2];
#line 165 "MB04OY.f"
	a[j * a_dim1 + 1] -= sum * *tau;
#line 166 "MB04OY.f"
	b[j * b_dim1 + 1] -= sum * t1;
#line 167 "MB04OY.f"
	b[j * b_dim1 + 2] -= sum * t2;
#line 168 "MB04OY.f"
/* L60: */
#line 168 "MB04OY.f"
    }
#line 169 "MB04OY.f"
    goto L210;
#line 170 "MB04OY.f"
L70:

/*     Special code for 4 x 4 Householder */

#line 174 "MB04OY.f"
    v1 = v[1];
#line 175 "MB04OY.f"
    t1 = *tau * v1;
#line 176 "MB04OY.f"
    v2 = v[2];
#line 177 "MB04OY.f"
    t2 = *tau * v2;
#line 178 "MB04OY.f"
    v3 = v[3];
#line 179 "MB04OY.f"
    t3 = *tau * v3;
#line 180 "MB04OY.f"
    i__1 = *n;
#line 180 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 181 "MB04OY.f"
	sum = a[j * a_dim1 + 1] + v1 * b[j * b_dim1 + 1] + v2 * b[j * b_dim1 
		+ 2] + v3 * b[j * b_dim1 + 3];
#line 182 "MB04OY.f"
	a[j * a_dim1 + 1] -= sum * *tau;
#line 183 "MB04OY.f"
	b[j * b_dim1 + 1] -= sum * t1;
#line 184 "MB04OY.f"
	b[j * b_dim1 + 2] -= sum * t2;
#line 185 "MB04OY.f"
	b[j * b_dim1 + 3] -= sum * t3;
#line 186 "MB04OY.f"
/* L80: */
#line 186 "MB04OY.f"
    }
#line 187 "MB04OY.f"
    goto L210;
#line 188 "MB04OY.f"
L90:

/*     Special code for 5 x 5 Householder */

#line 192 "MB04OY.f"
    v1 = v[1];
#line 193 "MB04OY.f"
    t1 = *tau * v1;
#line 194 "MB04OY.f"
    v2 = v[2];
#line 195 "MB04OY.f"
    t2 = *tau * v2;
#line 196 "MB04OY.f"
    v3 = v[3];
#line 197 "MB04OY.f"
    t3 = *tau * v3;
#line 198 "MB04OY.f"
    v4 = v[4];
#line 199 "MB04OY.f"
    t4 = *tau * v4;
#line 200 "MB04OY.f"
    i__1 = *n;
#line 200 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 201 "MB04OY.f"
	sum = a[j * a_dim1 + 1] + v1 * b[j * b_dim1 + 1] + v2 * b[j * b_dim1 
		+ 2] + v3 * b[j * b_dim1 + 3] + v4 * b[j * b_dim1 + 4];
#line 203 "MB04OY.f"
	a[j * a_dim1 + 1] -= sum * *tau;
#line 204 "MB04OY.f"
	b[j * b_dim1 + 1] -= sum * t1;
#line 205 "MB04OY.f"
	b[j * b_dim1 + 2] -= sum * t2;
#line 206 "MB04OY.f"
	b[j * b_dim1 + 3] -= sum * t3;
#line 207 "MB04OY.f"
	b[j * b_dim1 + 4] -= sum * t4;
#line 208 "MB04OY.f"
/* L100: */
#line 208 "MB04OY.f"
    }
#line 209 "MB04OY.f"
    goto L210;
#line 210 "MB04OY.f"
L110:

/*     Special code for 6 x 6 Householder */

#line 214 "MB04OY.f"
    v1 = v[1];
#line 215 "MB04OY.f"
    t1 = *tau * v1;
#line 216 "MB04OY.f"
    v2 = v[2];
#line 217 "MB04OY.f"
    t2 = *tau * v2;
#line 218 "MB04OY.f"
    v3 = v[3];
#line 219 "MB04OY.f"
    t3 = *tau * v3;
#line 220 "MB04OY.f"
    v4 = v[4];
#line 221 "MB04OY.f"
    t4 = *tau * v4;
#line 222 "MB04OY.f"
    v5 = v[5];
#line 223 "MB04OY.f"
    t5 = *tau * v5;
#line 224 "MB04OY.f"
    i__1 = *n;
#line 224 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 225 "MB04OY.f"
	sum = a[j * a_dim1 + 1] + v1 * b[j * b_dim1 + 1] + v2 * b[j * b_dim1 
		+ 2] + v3 * b[j * b_dim1 + 3] + v4 * b[j * b_dim1 + 4] + v5 * 
		b[j * b_dim1 + 5];
#line 227 "MB04OY.f"
	a[j * a_dim1 + 1] -= sum * *tau;
#line 228 "MB04OY.f"
	b[j * b_dim1 + 1] -= sum * t1;
#line 229 "MB04OY.f"
	b[j * b_dim1 + 2] -= sum * t2;
#line 230 "MB04OY.f"
	b[j * b_dim1 + 3] -= sum * t3;
#line 231 "MB04OY.f"
	b[j * b_dim1 + 4] -= sum * t4;
#line 232 "MB04OY.f"
	b[j * b_dim1 + 5] -= sum * t5;
#line 233 "MB04OY.f"
/* L120: */
#line 233 "MB04OY.f"
    }
#line 234 "MB04OY.f"
    goto L210;
#line 235 "MB04OY.f"
L130:

/*     Special code for 7 x 7 Householder */

#line 239 "MB04OY.f"
    v1 = v[1];
#line 240 "MB04OY.f"
    t1 = *tau * v1;
#line 241 "MB04OY.f"
    v2 = v[2];
#line 242 "MB04OY.f"
    t2 = *tau * v2;
#line 243 "MB04OY.f"
    v3 = v[3];
#line 244 "MB04OY.f"
    t3 = *tau * v3;
#line 245 "MB04OY.f"
    v4 = v[4];
#line 246 "MB04OY.f"
    t4 = *tau * v4;
#line 247 "MB04OY.f"
    v5 = v[5];
#line 248 "MB04OY.f"
    t5 = *tau * v5;
#line 249 "MB04OY.f"
    v6 = v[6];
#line 250 "MB04OY.f"
    t6 = *tau * v6;
#line 251 "MB04OY.f"
    i__1 = *n;
#line 251 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 252 "MB04OY.f"
	sum = a[j * a_dim1 + 1] + v1 * b[j * b_dim1 + 1] + v2 * b[j * b_dim1 
		+ 2] + v3 * b[j * b_dim1 + 3] + v4 * b[j * b_dim1 + 4] + v5 * 
		b[j * b_dim1 + 5] + v6 * b[j * b_dim1 + 6];
#line 254 "MB04OY.f"
	a[j * a_dim1 + 1] -= sum * *tau;
#line 255 "MB04OY.f"
	b[j * b_dim1 + 1] -= sum * t1;
#line 256 "MB04OY.f"
	b[j * b_dim1 + 2] -= sum * t2;
#line 257 "MB04OY.f"
	b[j * b_dim1 + 3] -= sum * t3;
#line 258 "MB04OY.f"
	b[j * b_dim1 + 4] -= sum * t4;
#line 259 "MB04OY.f"
	b[j * b_dim1 + 5] -= sum * t5;
#line 260 "MB04OY.f"
	b[j * b_dim1 + 6] -= sum * t6;
#line 261 "MB04OY.f"
/* L140: */
#line 261 "MB04OY.f"
    }
#line 262 "MB04OY.f"
    goto L210;
#line 263 "MB04OY.f"
L150:

/*     Special code for 8 x 8 Householder */

#line 267 "MB04OY.f"
    v1 = v[1];
#line 268 "MB04OY.f"
    t1 = *tau * v1;
#line 269 "MB04OY.f"
    v2 = v[2];
#line 270 "MB04OY.f"
    t2 = *tau * v2;
#line 271 "MB04OY.f"
    v3 = v[3];
#line 272 "MB04OY.f"
    t3 = *tau * v3;
#line 273 "MB04OY.f"
    v4 = v[4];
#line 274 "MB04OY.f"
    t4 = *tau * v4;
#line 275 "MB04OY.f"
    v5 = v[5];
#line 276 "MB04OY.f"
    t5 = *tau * v5;
#line 277 "MB04OY.f"
    v6 = v[6];
#line 278 "MB04OY.f"
    t6 = *tau * v6;
#line 279 "MB04OY.f"
    v7 = v[7];
#line 280 "MB04OY.f"
    t7 = *tau * v7;
#line 281 "MB04OY.f"
    i__1 = *n;
#line 281 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 282 "MB04OY.f"
	sum = a[j * a_dim1 + 1] + v1 * b[j * b_dim1 + 1] + v2 * b[j * b_dim1 
		+ 2] + v3 * b[j * b_dim1 + 3] + v4 * b[j * b_dim1 + 4] + v5 * 
		b[j * b_dim1 + 5] + v6 * b[j * b_dim1 + 6] + v7 * b[j * 
		b_dim1 + 7];
#line 285 "MB04OY.f"
	a[j * a_dim1 + 1] -= sum * *tau;
#line 286 "MB04OY.f"
	b[j * b_dim1 + 1] -= sum * t1;
#line 287 "MB04OY.f"
	b[j * b_dim1 + 2] -= sum * t2;
#line 288 "MB04OY.f"
	b[j * b_dim1 + 3] -= sum * t3;
#line 289 "MB04OY.f"
	b[j * b_dim1 + 4] -= sum * t4;
#line 290 "MB04OY.f"
	b[j * b_dim1 + 5] -= sum * t5;
#line 291 "MB04OY.f"
	b[j * b_dim1 + 6] -= sum * t6;
#line 292 "MB04OY.f"
	b[j * b_dim1 + 7] -= sum * t7;
#line 293 "MB04OY.f"
/* L160: */
#line 293 "MB04OY.f"
    }
#line 294 "MB04OY.f"
    goto L210;
#line 295 "MB04OY.f"
L170:

/*     Special code for 9 x 9 Householder */

#line 299 "MB04OY.f"
    v1 = v[1];
#line 300 "MB04OY.f"
    t1 = *tau * v1;
#line 301 "MB04OY.f"
    v2 = v[2];
#line 302 "MB04OY.f"
    t2 = *tau * v2;
#line 303 "MB04OY.f"
    v3 = v[3];
#line 304 "MB04OY.f"
    t3 = *tau * v3;
#line 305 "MB04OY.f"
    v4 = v[4];
#line 306 "MB04OY.f"
    t4 = *tau * v4;
#line 307 "MB04OY.f"
    v5 = v[5];
#line 308 "MB04OY.f"
    t5 = *tau * v5;
#line 309 "MB04OY.f"
    v6 = v[6];
#line 310 "MB04OY.f"
    t6 = *tau * v6;
#line 311 "MB04OY.f"
    v7 = v[7];
#line 312 "MB04OY.f"
    t7 = *tau * v7;
#line 313 "MB04OY.f"
    v8 = v[8];
#line 314 "MB04OY.f"
    t8 = *tau * v8;
#line 315 "MB04OY.f"
    i__1 = *n;
#line 315 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 316 "MB04OY.f"
	sum = a[j * a_dim1 + 1] + v1 * b[j * b_dim1 + 1] + v2 * b[j * b_dim1 
		+ 2] + v3 * b[j * b_dim1 + 3] + v4 * b[j * b_dim1 + 4] + v5 * 
		b[j * b_dim1 + 5] + v6 * b[j * b_dim1 + 6] + v7 * b[j * 
		b_dim1 + 7] + v8 * b[j * b_dim1 + 8];
#line 319 "MB04OY.f"
	a[j * a_dim1 + 1] -= sum * *tau;
#line 320 "MB04OY.f"
	b[j * b_dim1 + 1] -= sum * t1;
#line 321 "MB04OY.f"
	b[j * b_dim1 + 2] -= sum * t2;
#line 322 "MB04OY.f"
	b[j * b_dim1 + 3] -= sum * t3;
#line 323 "MB04OY.f"
	b[j * b_dim1 + 4] -= sum * t4;
#line 324 "MB04OY.f"
	b[j * b_dim1 + 5] -= sum * t5;
#line 325 "MB04OY.f"
	b[j * b_dim1 + 6] -= sum * t6;
#line 326 "MB04OY.f"
	b[j * b_dim1 + 7] -= sum * t7;
#line 327 "MB04OY.f"
	b[j * b_dim1 + 8] -= sum * t8;
#line 328 "MB04OY.f"
/* L180: */
#line 328 "MB04OY.f"
    }
#line 329 "MB04OY.f"
    goto L210;
#line 330 "MB04OY.f"
L190:

/*     Special code for 10 x 10 Householder */

#line 334 "MB04OY.f"
    v1 = v[1];
#line 335 "MB04OY.f"
    t1 = *tau * v1;
#line 336 "MB04OY.f"
    v2 = v[2];
#line 337 "MB04OY.f"
    t2 = *tau * v2;
#line 338 "MB04OY.f"
    v3 = v[3];
#line 339 "MB04OY.f"
    t3 = *tau * v3;
#line 340 "MB04OY.f"
    v4 = v[4];
#line 341 "MB04OY.f"
    t4 = *tau * v4;
#line 342 "MB04OY.f"
    v5 = v[5];
#line 343 "MB04OY.f"
    t5 = *tau * v5;
#line 344 "MB04OY.f"
    v6 = v[6];
#line 345 "MB04OY.f"
    t6 = *tau * v6;
#line 346 "MB04OY.f"
    v7 = v[7];
#line 347 "MB04OY.f"
    t7 = *tau * v7;
#line 348 "MB04OY.f"
    v8 = v[8];
#line 349 "MB04OY.f"
    t8 = *tau * v8;
#line 350 "MB04OY.f"
    v9 = v[9];
#line 351 "MB04OY.f"
    t9 = *tau * v9;
#line 352 "MB04OY.f"
    i__1 = *n;
#line 352 "MB04OY.f"
    for (j = 1; j <= i__1; ++j) {
#line 353 "MB04OY.f"
	sum = a[j * a_dim1 + 1] + v1 * b[j * b_dim1 + 1] + v2 * b[j * b_dim1 
		+ 2] + v3 * b[j * b_dim1 + 3] + v4 * b[j * b_dim1 + 4] + v5 * 
		b[j * b_dim1 + 5] + v6 * b[j * b_dim1 + 6] + v7 * b[j * 
		b_dim1 + 7] + v8 * b[j * b_dim1 + 8] + v9 * b[j * b_dim1 + 9];
#line 356 "MB04OY.f"
	a[j * a_dim1 + 1] -= sum * *tau;
#line 357 "MB04OY.f"
	b[j * b_dim1 + 1] -= sum * t1;
#line 358 "MB04OY.f"
	b[j * b_dim1 + 2] -= sum * t2;
#line 359 "MB04OY.f"
	b[j * b_dim1 + 3] -= sum * t3;
#line 360 "MB04OY.f"
	b[j * b_dim1 + 4] -= sum * t4;
#line 361 "MB04OY.f"
	b[j * b_dim1 + 5] -= sum * t5;
#line 362 "MB04OY.f"
	b[j * b_dim1 + 6] -= sum * t6;
#line 363 "MB04OY.f"
	b[j * b_dim1 + 7] -= sum * t7;
#line 364 "MB04OY.f"
	b[j * b_dim1 + 8] -= sum * t8;
#line 365 "MB04OY.f"
	b[j * b_dim1 + 9] -= sum * t9;
#line 366 "MB04OY.f"
/* L200: */
#line 366 "MB04OY.f"
    }
#line 367 "MB04OY.f"
L210:
#line 368 "MB04OY.f"
    return 0;
/* *** Last line of MB04OY *** */
} /* mb04oy_ */

