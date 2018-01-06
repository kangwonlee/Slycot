#line 1 "SB01MD.f"
/* SB01MD.f -- translated by f2c (version 20100827).
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

#line 1 "SB01MD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b23 = 1.;
static doublereal c_b25 = 0.;

/* Subroutine */ int sb01md_(integer *ncont, integer *n, doublereal *a, 
	integer *lda, doublereal *b, doublereal *wr, doublereal *wi, 
	doublereal *z__, integer *ldz, doublereal *g, doublereal *dwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, k, l;
    static doublereal p, q, r__, s, t, b1;
    static integer ni, ll, nj, nl, im1, lp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dscal_(
	    integer *, doublereal *, doublereal *, integer *), dgemv_(char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    static logical compl;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static integer ncont2;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xerbla_(char *, integer *, ftnlen);


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

/*     To determine the one-dimensional state feedback matrix G of the */
/*     linear time-invariant single-input system */

/*           dX/dt = A * X + B * U, */

/*     where A is an NCONT-by-NCONT matrix and B is an NCONT element */
/*     vector such that the closed-loop system */

/*           dX/dt = (A - B * G) * X */

/*     has desired poles. The system must be preliminarily reduced */
/*     to orthogonal canonical form using the SLICOT Library routine */
/*     AB01MD. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NCONT   (input) INTEGER */
/*             The order of the matrix A as produced by SLICOT Library */
/*             routine AB01MD.  NCONT >= 0. */

/*     N       (input) INTEGER */
/*             The order of the matrix Z.  N >= NCONT. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDA,NCONT) */
/*             On entry, the leading NCONT-by-NCONT part of this array */
/*             must contain the canonical form of the state dynamics */
/*             matrix A as produced by SLICOT Library routine AB01MD. */
/*             On exit, the leading NCONT-by-NCONT part of this array */
/*             contains the upper quasi-triangular form S of the closed- */
/*             loop system matrix (A - B * G), that is triangular except */
/*             for possible 2-by-2 diagonal blocks. */
/*             (To reconstruct the closed-loop system matrix see */
/*             FURTHER COMMENTS below.) */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,NCONT). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (NCONT) */
/*             On entry, this array must contain the canonical form of */
/*             the input/state vector B as produced by SLICOT Library */
/*             routine AB01MD. */
/*             On exit, this array contains the transformed vector Z * B */
/*             of the closed-loop system. */

/*     WR      (input) DOUBLE PRECISION array, dimension (NCONT) */
/*     WI      (input) DOUBLE PRECISION array, dimension (NCONT) */
/*             These arrays must contain the real and imaginary parts, */
/*             respectively, of the desired poles of the closed-loop */
/*             system. The poles can be unordered, except that complex */
/*             conjugate pairs of poles must appear consecutively. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the orthogonal transformation matrix as produced */
/*             by SLICOT Library routine AB01MD, which reduces the system */
/*             to canonical form. */
/*             On exit, the leading NCONT-by-NCONT part of this array */
/*             contains the orthogonal matrix Z which reduces the closed- */
/*             loop system matrix (A - B * G) to upper quasi-triangular */
/*             form. */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z.  LDZ >= MAX(1,N). */

/*     G       (output) DOUBLE PRECISION array, dimension (NCONT) */
/*             This array contains the one-dimensional state feedback */
/*             matrix G of the original system. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (3*NCONT) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The method is based on the orthogonal reduction of the closed-loop */
/*     system matrix (A - B * G) to upper quasi-triangular form S whose */
/*     1-by-1 and 2-by-2 diagonal blocks correspond to the desired poles. */
/*     That is, S = Z'*(A - B * G)*Z, where Z is an orthogonal matrix. */

/*     REFERENCES */

/*     [1] Petkov, P. Hr. */
/*         A Computational Algorithm for Pole Assignment of Linear */
/*         Single Input Systems. */
/*         Internal Report 81/2, Control Systems Research Group, School */
/*         of Electronic Engineering and Computer Science, Kingston */
/*         Polytechnic, 1981. */

/*     NUMERICAL ASPECTS */
/*                                   3 */
/*     The algorithm requires 0(NCONT ) operations and is backward */
/*     stable. */

/*     FURTHER COMMENTS */

/*     If required, the closed-loop system matrix (A - B * G) can be */
/*     formed from the matrix product Z * S * Z' (where S and Z are the */
/*     matrices output in arrays A and Z respectively). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB01AD by Control Systems Research */
/*     Group, Kingston Polytechnic, United Kingdom, May 1981. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Closed loop spectrum, closed loop systems, eigenvalue assignment, */
/*     orthogonal canonical form, orthogonal transformation, pole */
/*     placement, Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 173 "SB01MD.f"
    /* Parameter adjustments */
#line 173 "SB01MD.f"
    a_dim1 = *lda;
#line 173 "SB01MD.f"
    a_offset = 1 + a_dim1;
#line 173 "SB01MD.f"
    a -= a_offset;
#line 173 "SB01MD.f"
    --b;
#line 173 "SB01MD.f"
    --wr;
#line 173 "SB01MD.f"
    --wi;
#line 173 "SB01MD.f"
    z_dim1 = *ldz;
#line 173 "SB01MD.f"
    z_offset = 1 + z_dim1;
#line 173 "SB01MD.f"
    z__ -= z_offset;
#line 173 "SB01MD.f"
    --g;
#line 173 "SB01MD.f"
    --dwork;
#line 173 "SB01MD.f"

#line 173 "SB01MD.f"
    /* Function Body */
#line 173 "SB01MD.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 177 "SB01MD.f"
    if (*ncont < 0) {
#line 178 "SB01MD.f"
	*info = -1;
#line 179 "SB01MD.f"
    } else if (*n < *ncont) {
#line 180 "SB01MD.f"
	*info = -2;
#line 181 "SB01MD.f"
    } else if (*lda < max(1,*ncont)) {
#line 182 "SB01MD.f"
	*info = -4;
#line 183 "SB01MD.f"
    } else if (*ldz < max(1,*n)) {
#line 184 "SB01MD.f"
	*info = -9;
#line 185 "SB01MD.f"
    }

#line 187 "SB01MD.f"
    if (*info != 0) {

/*        Error return */

#line 191 "SB01MD.f"
	i__1 = -(*info);
#line 191 "SB01MD.f"
	xerbla_("SB01MD", &i__1, (ftnlen)6);
#line 192 "SB01MD.f"
	return 0;
#line 193 "SB01MD.f"
    }

/*     Quick return if possible. */

#line 197 "SB01MD.f"
    if (*ncont == 0 || *n == 0) {
#line 197 "SB01MD.f"
	return 0;
#line 197 "SB01MD.f"
    }

/*     Return if the system is not complete controllable. */

#line 202 "SB01MD.f"
    if (b[1] == 0.) {
#line 202 "SB01MD.f"
	return 0;
#line 202 "SB01MD.f"
    }

#line 205 "SB01MD.f"
    if (*ncont == 1) {

/*        1-by-1 case. */

#line 209 "SB01MD.f"
	p = a[a_dim1 + 1] - wr[1];
#line 210 "SB01MD.f"
	a[a_dim1 + 1] = wr[1];
#line 211 "SB01MD.f"
	g[1] = p / b[1];
#line 212 "SB01MD.f"
	z__[z_dim1 + 1] = 1.;
#line 213 "SB01MD.f"
	return 0;
#line 214 "SB01MD.f"
    }

/*     General case.  Save the contents of WI in DWORK. */

#line 218 "SB01MD.f"
    ncont2 = *ncont << 1;
#line 219 "SB01MD.f"
    dcopy_(ncont, &wi[1], &c__1, &dwork[ncont2 + 1], &c__1);

#line 221 "SB01MD.f"
    b1 = b[1];
#line 222 "SB01MD.f"
    b[1] = 1.;
#line 223 "SB01MD.f"
    l = 0;
#line 224 "SB01MD.f"
    ll = 0;
#line 225 "SB01MD.f"
L20:
#line 226 "SB01MD.f"
    ++l;
#line 227 "SB01MD.f"
    ++ll;
#line 228 "SB01MD.f"
    compl = dwork[ncont2 + l] != 0.;
#line 229 "SB01MD.f"
    if (l != *ncont) {
#line 230 "SB01MD.f"
	lp1 = l + 1;
#line 231 "SB01MD.f"
	nl = *ncont - l;
#line 232 "SB01MD.f"
	if (ll != 2) {
#line 233 "SB01MD.f"
	    if (compl) {

/*              Compute complex eigenvector. */

#line 237 "SB01MD.f"
		dwork[*ncont] = 1.;
#line 238 "SB01MD.f"
		dwork[ncont2] = 1.;
#line 239 "SB01MD.f"
		p = wr[l];
#line 240 "SB01MD.f"
		t = dwork[ncont2 + l];
#line 241 "SB01MD.f"
		q = t * dwork[ncont2 + lp1];
#line 242 "SB01MD.f"
		dwork[ncont2 + l] = 1.;
#line 243 "SB01MD.f"
		dwork[ncont2 + lp1] = q;

#line 245 "SB01MD.f"
		i__1 = lp1;
#line 245 "SB01MD.f"
		for (i__ = *ncont; i__ >= i__1; --i__) {
#line 246 "SB01MD.f"
		    im1 = i__ - 1;
#line 247 "SB01MD.f"
		    i__2 = *ncont - im1;
#line 247 "SB01MD.f"
		    dwork[im1] = (p * dwork[i__] + q * dwork[*ncont + i__] - 
			    ddot_(&i__2, &a[i__ + i__ * a_dim1], lda, &dwork[
			    i__], &c__1)) / a[i__ + im1 * a_dim1];
#line 250 "SB01MD.f"
		    i__2 = *ncont - im1;
#line 250 "SB01MD.f"
		    dwork[*ncont + im1] = (p * dwork[*ncont + i__] + dwork[
			    i__] - ddot_(&i__2, &a[i__ + i__ * a_dim1], lda, &
			    dwork[*ncont + i__], &c__1)) / a[i__ + im1 * 
			    a_dim1];
#line 253 "SB01MD.f"
/* L40: */
#line 253 "SB01MD.f"
		}

#line 255 "SB01MD.f"
	    } else {

/*              Compute real eigenvector. */

#line 259 "SB01MD.f"
		dwork[*ncont] = 1.;
#line 260 "SB01MD.f"
		p = wr[l];

#line 262 "SB01MD.f"
		i__1 = lp1;
#line 262 "SB01MD.f"
		for (i__ = *ncont; i__ >= i__1; --i__) {
#line 263 "SB01MD.f"
		    im1 = i__ - 1;
#line 264 "SB01MD.f"
		    i__2 = *ncont - im1;
#line 264 "SB01MD.f"
		    dwork[im1] = (p * dwork[i__] - ddot_(&i__2, &a[i__ + i__ *
			     a_dim1], lda, &dwork[i__], &c__1)) / a[i__ + im1 
			    * a_dim1];
#line 267 "SB01MD.f"
/* L60: */
#line 267 "SB01MD.f"
		}

#line 269 "SB01MD.f"
	    }
#line 270 "SB01MD.f"
	}

/*        Transform eigenvector. */

#line 274 "SB01MD.f"
	i__1 = l;
#line 274 "SB01MD.f"
	for (k = *ncont - 1; k >= i__1; --k) {
#line 275 "SB01MD.f"
	    if (ll != 2) {
#line 276 "SB01MD.f"
		r__ = dwork[k];
#line 277 "SB01MD.f"
		s = dwork[k + 1];
#line 278 "SB01MD.f"
	    } else {
#line 279 "SB01MD.f"
		r__ = dwork[*ncont + k];
#line 280 "SB01MD.f"
		s = dwork[*ncont + k + 1];
#line 281 "SB01MD.f"
	    }
#line 282 "SB01MD.f"
	    dlartg_(&r__, &s, &p, &q, &t);
#line 283 "SB01MD.f"
	    dwork[k] = t;
#line 284 "SB01MD.f"
	    if (ll != 2) {
/* Computing MAX */
#line 285 "SB01MD.f"
		i__2 = k - 1;
#line 285 "SB01MD.f"
		nj = max(i__2,l);
#line 286 "SB01MD.f"
	    } else {
#line 287 "SB01MD.f"
		dwork[*ncont + k] = t;
#line 288 "SB01MD.f"
		nj = l - 1;
#line 289 "SB01MD.f"
	    }

/*           Transform  A. */

#line 293 "SB01MD.f"
	    i__2 = *ncont - nj + 1;
#line 293 "SB01MD.f"
	    drot_(&i__2, &a[k + nj * a_dim1], lda, &a[k + 1 + nj * a_dim1], 
		    lda, &p, &q);

#line 295 "SB01MD.f"
	    if (compl && ll == 1) {
#line 296 "SB01MD.f"
		ni = *ncont;
#line 297 "SB01MD.f"
	    } else {
/* Computing MIN */
#line 298 "SB01MD.f"
		i__2 = k + 2;
#line 298 "SB01MD.f"
		ni = min(i__2,*ncont);
#line 299 "SB01MD.f"
	    }
#line 300 "SB01MD.f"
	    drot_(&ni, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1], &
		    c__1, &p, &q);

#line 302 "SB01MD.f"
	    if (k == l) {

/*              Transform  B. */

#line 306 "SB01MD.f"
		t = b[k];
#line 307 "SB01MD.f"
		b[k] = p * t;
#line 308 "SB01MD.f"
		b[k + 1] = -q * t;
#line 309 "SB01MD.f"
	    }

/*           Accumulate transformations. */

#line 313 "SB01MD.f"
	    drot_(ncont, &z__[k * z_dim1 + 1], &c__1, &z__[(k + 1) * z_dim1 + 
		    1], &c__1, &p, &q);

#line 315 "SB01MD.f"
	    if (compl && ll != 2) {
#line 316 "SB01MD.f"
		t = dwork[*ncont + k];
#line 317 "SB01MD.f"
		dwork[*ncont + k] = p * t + q * dwork[*ncont + k + 1];
#line 318 "SB01MD.f"
		dwork[*ncont + k + 1] = p * dwork[*ncont + k + 1] - q * t;
#line 319 "SB01MD.f"
	    }
#line 320 "SB01MD.f"
/* L80: */
#line 320 "SB01MD.f"
	}

#line 322 "SB01MD.f"
    }

#line 324 "SB01MD.f"
    if (! compl) {

/*        Find one element of  G. */

#line 328 "SB01MD.f"
	k = l;
#line 329 "SB01MD.f"
	r__ = b[l];
#line 330 "SB01MD.f"
	if (l != *ncont) {
#line 331 "SB01MD.f"
	    if ((d__1 = b[lp1], abs(d__1)) > (d__2 = b[l], abs(d__2))) {
#line 332 "SB01MD.f"
		k = lp1;
#line 333 "SB01MD.f"
		r__ = b[lp1];
#line 334 "SB01MD.f"
	    }
#line 335 "SB01MD.f"
	}
#line 336 "SB01MD.f"
	p = a[k + l * a_dim1];
#line 337 "SB01MD.f"
	if (k == l) {
#line 337 "SB01MD.f"
	    p -= wr[l];
#line 337 "SB01MD.f"
	}
#line 338 "SB01MD.f"
	p /= r__;

#line 340 "SB01MD.f"
	d__1 = -p;
#line 340 "SB01MD.f"
	daxpy_(&lp1, &d__1, &b[1], &c__1, &a[l * a_dim1 + 1], &c__1);

#line 342 "SB01MD.f"
	g[l] = p / b1;
#line 343 "SB01MD.f"
	if (l != *ncont) {
#line 344 "SB01MD.f"
	    ll = 0;
#line 345 "SB01MD.f"
	    goto L20;
#line 346 "SB01MD.f"
	}
#line 347 "SB01MD.f"
    } else if (ll == 1) {
#line 348 "SB01MD.f"
	goto L20;
#line 349 "SB01MD.f"
    } else {

/*        Find two elements of  G. */

#line 353 "SB01MD.f"
	k = l;
#line 354 "SB01MD.f"
	r__ = b[l];
#line 355 "SB01MD.f"
	if (l != *ncont) {
#line 356 "SB01MD.f"
	    if ((d__1 = b[lp1], abs(d__1)) > (d__2 = b[l], abs(d__2))) {
#line 357 "SB01MD.f"
		k = lp1;
#line 358 "SB01MD.f"
		r__ = b[lp1];
#line 359 "SB01MD.f"
	    }
#line 360 "SB01MD.f"
	}
#line 361 "SB01MD.f"
	p = a[k + (l - 1) * a_dim1];
#line 362 "SB01MD.f"
	q = a[k + l * a_dim1];
#line 363 "SB01MD.f"
	if (k == l) {
#line 364 "SB01MD.f"
	    p -= dwork[*ncont + l] / dwork[l - 1] * dwork[ncont2 + l];
#line 365 "SB01MD.f"
	    q = q - wr[l] + dwork[*ncont + l - 1] / dwork[l - 1] * dwork[
		    ncont2 + l];
#line 367 "SB01MD.f"
	}
#line 368 "SB01MD.f"
	p /= r__;
#line 369 "SB01MD.f"
	q /= r__;

#line 371 "SB01MD.f"
	d__1 = -p;
#line 371 "SB01MD.f"
	daxpy_(&lp1, &d__1, &b[1], &c__1, &a[(l - 1) * a_dim1 + 1], &c__1);
#line 372 "SB01MD.f"
	d__1 = -q;
#line 372 "SB01MD.f"
	daxpy_(&lp1, &d__1, &b[1], &c__1, &a[l * a_dim1 + 1], &c__1);

#line 374 "SB01MD.f"
	g[l - 1] = p / b1;
#line 375 "SB01MD.f"
	g[l] = q / b1;
#line 376 "SB01MD.f"
	if (l != *ncont) {
#line 377 "SB01MD.f"
	    ll = 0;
#line 378 "SB01MD.f"
	    goto L20;
#line 379 "SB01MD.f"
	}
#line 380 "SB01MD.f"
    }

/*     Transform  G. */

#line 384 "SB01MD.f"
    dgemv_("No transpose", ncont, ncont, &c_b23, &z__[z_offset], ldz, &g[1], &
	    c__1, &c_b25, &dwork[1], &c__1, (ftnlen)12);
#line 386 "SB01MD.f"
    dcopy_(ncont, &dwork[1], &c__1, &g[1], &c__1);
#line 387 "SB01MD.f"
    dscal_(ncont, &b1, &b[1], &c__1);

/*     Annihilate A after the first subdiagonal. */

#line 391 "SB01MD.f"
    if (*ncont > 2) {
#line 391 "SB01MD.f"
	i__1 = *ncont - 2;
#line 391 "SB01MD.f"
	i__2 = *ncont - 2;
#line 391 "SB01MD.f"
	dlaset_("Lower", &i__1, &i__2, &c_b25, &c_b25, &a[a_dim1 + 3], lda, (
		ftnlen)5);
#line 391 "SB01MD.f"
    }

#line 395 "SB01MD.f"
    return 0;
/* *** Last line of SB01MD *** */
} /* sb01md_ */

