#line 1 "MB03YT.f"
/* MB03YT.f -- translated by f2c (version 20100827).
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

#line 1 "MB03YT.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int mb03yt_(doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *csl, doublereal *snl, doublereal *csr, doublereal *
	snr)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static doublereal r__, t, h1, h2, h3, wi, qq, rr, wr1, wr2, ulp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dlag2_(
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal anorm, bnorm, scale1, scale2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlasv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal safmin;
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

/*     To compute the periodic Schur factorization of a real 2-by-2 */
/*     matrix pair (A,B) where B is upper triangular. This routine */
/*     computes orthogonal (rotation) matrices given by CSL, SNL and CSR, */
/*     SNR such that */

/*     1) if the pair (A,B) has two real eigenvalues, then */

/*        [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ] */
/*        [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ] */

/*        [ b11 b12 ] := [  CSR  SNR ] [ b11 b12 ] [  CSL -SNL ] */
/*        [  0  b22 ]    [ -SNR  CSR ] [  0  b22 ] [  SNL  CSL ], */

/*     2) if the pair (A,B) has a pair of complex conjugate eigenvalues, */
/*        then */

/*        [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ] */
/*        [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ] */

/*        [ b11  0  ] := [  CSR  SNR ] [ b11 b12 ] [  CSL -SNL ] */
/*        [  0  b22 ]    [ -SNR  CSR ] [  0  b22 ] [  SNL  CSL ]. */

/*     This is a modified version of the LAPACK routine DLAGV2 for */
/*     computing the real, generalized Schur decomposition of a */
/*     two-by-two matrix pencil. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,2) */
/*             On entry, the leading 2-by-2 part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading 2-by-2 part of this array contains */
/*             the matrix A of the pair in periodic Schur form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= 2. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,2) */
/*             On entry, the leading 2-by-2 part of this array must */
/*             contain the upper triangular matrix B. */
/*             On exit, the leading 2-by-2 part of this array contains */
/*             the matrix B of the pair in periodic Schur form. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= 2. */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (2) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (2) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (2) */
/*             (ALPHAR(k)+i*ALPHAI(k))*BETA(k) are the eigenvalues of the */
/*             pair (A,B), k=1,2, i = sqrt(-1). ALPHAI(1) >= 0. */

/*     CSL     (output) DOUBLE PRECISION */
/*             The cosine of the first rotation matrix. */

/*     SNL     (output) DOUBLE PRECISION */
/*             The sine of the first rotation matrix. */

/*     CSR     (output) DOUBLE PRECISION */
/*             The cosine of the second rotation matrix. */

/*     SNR     (output) DOUBLE PRECISION */
/*             The sine of the second rotation matrix. */

/*     REFERENCES */

/*     [1] Van Loan, C. */
/*         Generalized Singular Values with Algorithms and Applications. */
/*         Ph. D. Thesis, University of Michigan, 1973. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLAPV2). */
/*     V. Sima, July 2008, May 2009. */

/*     KEYWORDS */

/*     Eigenvalue, periodic Schur form */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

#line 134 "MB03YT.f"
    /* Parameter adjustments */
#line 134 "MB03YT.f"
    a_dim1 = *lda;
#line 134 "MB03YT.f"
    a_offset = 1 + a_dim1;
#line 134 "MB03YT.f"
    a -= a_offset;
#line 134 "MB03YT.f"
    b_dim1 = *ldb;
#line 134 "MB03YT.f"
    b_offset = 1 + b_dim1;
#line 134 "MB03YT.f"
    b -= b_offset;
#line 134 "MB03YT.f"
    --alphar;
#line 134 "MB03YT.f"
    --alphai;
#line 134 "MB03YT.f"
    --beta;
#line 134 "MB03YT.f"

#line 134 "MB03YT.f"
    /* Function Body */
#line 134 "MB03YT.f"
    safmin = dlamch_("S", (ftnlen)1);
#line 135 "MB03YT.f"
    ulp = dlamch_("P", (ftnlen)1);

/*     Scale A. */

/* Computing MAX */
#line 139 "MB03YT.f"
    d__5 = (d__1 = a[a_dim1 + 1], abs(d__1)) + (d__2 = a[a_dim1 + 2], abs(
	    d__2)), d__6 = (d__3 = a[(a_dim1 << 1) + 1], abs(d__3)) + (d__4 = 
	    a[(a_dim1 << 1) + 2], abs(d__4)), d__5 = max(d__5,d__6);
#line 139 "MB03YT.f"
    anorm = max(d__5,safmin);
#line 141 "MB03YT.f"
    a[a_dim1 + 1] /= anorm;
#line 142 "MB03YT.f"
    a[(a_dim1 << 1) + 1] /= anorm;
#line 143 "MB03YT.f"
    a[a_dim1 + 2] /= anorm;
#line 144 "MB03YT.f"
    a[(a_dim1 << 1) + 2] /= anorm;

/*     Scale B. */

/* Computing MAX */
#line 148 "MB03YT.f"
    d__4 = (d__3 = b[b_dim1 + 1], abs(d__3)), d__5 = (d__1 = b[(b_dim1 << 1) 
	    + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 1) + 2], abs(d__2)), d__4 
	    = max(d__4,d__5);
#line 148 "MB03YT.f"
    bnorm = max(d__4,safmin);
#line 149 "MB03YT.f"
    b[b_dim1 + 1] /= bnorm;
#line 150 "MB03YT.f"
    b[(b_dim1 << 1) + 1] /= bnorm;
#line 151 "MB03YT.f"
    b[(b_dim1 << 1) + 2] /= bnorm;

/*     Check if A can be deflated. */

#line 155 "MB03YT.f"
    if ((d__1 = a[a_dim1 + 2], abs(d__1)) <= ulp) {
#line 156 "MB03YT.f"
	*csl = 1.;
#line 157 "MB03YT.f"
	*snl = 0.;
#line 158 "MB03YT.f"
	*csr = 1.;
#line 159 "MB03YT.f"
	*snr = 0.;
#line 160 "MB03YT.f"
	wi = 0.;
#line 161 "MB03YT.f"
	a[a_dim1 + 2] = 0.;
#line 162 "MB03YT.f"
	b[b_dim1 + 2] = 0.;

/*     Check if B is singular. */

#line 166 "MB03YT.f"
    } else if ((d__1 = b[b_dim1 + 1], abs(d__1)) <= ulp) {
#line 167 "MB03YT.f"
	dlartg_(&a[(a_dim1 << 1) + 2], &a[a_dim1 + 2], csr, snr, &t);
#line 168 "MB03YT.f"
	*snr = -(*snr);
#line 169 "MB03YT.f"
	drot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, csr,
		 snr);
#line 170 "MB03YT.f"
	drot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csr, snr);
#line 171 "MB03YT.f"
	*csl = 1.;
#line 172 "MB03YT.f"
	*snl = 0.;
#line 173 "MB03YT.f"
	wi = 0.;
#line 174 "MB03YT.f"
	a[a_dim1 + 2] = 0.;
#line 175 "MB03YT.f"
	b[b_dim1 + 1] = 0.;
#line 176 "MB03YT.f"
	b[b_dim1 + 2] = 0.;
#line 177 "MB03YT.f"
    } else if ((d__1 = b[(b_dim1 << 1) + 2], abs(d__1)) <= ulp) {
#line 178 "MB03YT.f"
	dlartg_(&a[a_dim1 + 1], &a[a_dim1 + 2], csl, snl, &r__);
#line 179 "MB03YT.f"
	*csr = 1.;
#line 180 "MB03YT.f"
	*snr = 0.;
#line 181 "MB03YT.f"
	wi = 0.;
#line 182 "MB03YT.f"
	drot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
#line 183 "MB03YT.f"
	drot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, csl,
		 snl);
#line 184 "MB03YT.f"
	a[a_dim1 + 2] = 0.;
#line 185 "MB03YT.f"
	b[b_dim1 + 2] = 0.;
#line 186 "MB03YT.f"
	b[(b_dim1 << 1) + 2] = 0.;
#line 187 "MB03YT.f"
    } else {

/*        B is nonsingular, first compute the eigenvalues of A / adj(B). */

#line 191 "MB03YT.f"
	r__ = b[b_dim1 + 1];
#line 192 "MB03YT.f"
	b[b_dim1 + 1] = b[(b_dim1 << 1) + 2];
#line 193 "MB03YT.f"
	b[(b_dim1 << 1) + 2] = r__;
#line 194 "MB03YT.f"
	b[(b_dim1 << 1) + 1] = -b[(b_dim1 << 1) + 1];
#line 195 "MB03YT.f"
	dlag2_(&a[a_offset], lda, &b[b_offset], ldb, &safmin, &scale1, &
		scale2, &wr1, &wr2, &wi);

#line 198 "MB03YT.f"
	if (wi == 0.) {

/*           Two real eigenvalues, compute s*A-w*B. */

#line 202 "MB03YT.f"
	    h1 = scale1 * a[a_dim1 + 1] - wr1 * b[b_dim1 + 1];
#line 203 "MB03YT.f"
	    h2 = scale1 * a[(a_dim1 << 1) + 1] - wr1 * b[(b_dim1 << 1) + 1];
#line 204 "MB03YT.f"
	    h3 = scale1 * a[(a_dim1 << 1) + 2] - wr1 * b[(b_dim1 << 1) + 2];

#line 206 "MB03YT.f"
	    rr = dlapy2_(&h1, &h2);
#line 207 "MB03YT.f"
	    d__1 = scale1 * a[a_dim1 + 2];
#line 207 "MB03YT.f"
	    qq = dlapy2_(&d__1, &h3);

#line 209 "MB03YT.f"
	    if (rr > qq) {

/*              Find right rotation matrix to zero 1,1 element of */
/*              (sA - wB). */

#line 214 "MB03YT.f"
		dlartg_(&h2, &h1, csr, snr, &t);

#line 216 "MB03YT.f"
	    } else {

/*              Find right rotation matrix to zero 2,1 element of */
/*              (sA - wB). */

#line 221 "MB03YT.f"
		d__1 = scale1 * a[a_dim1 + 2];
#line 221 "MB03YT.f"
		dlartg_(&h3, &d__1, csr, snr, &t);

#line 223 "MB03YT.f"
	    }

#line 225 "MB03YT.f"
	    *snr = -(*snr);
#line 226 "MB03YT.f"
	    drot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, 
		    csr, snr);
#line 227 "MB03YT.f"
	    drot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, 
		    csr, snr);

/*           Compute inf norms of A and B. */

/* Computing MAX */
#line 231 "MB03YT.f"
	    d__5 = (d__1 = a[a_dim1 + 1], abs(d__1)) + (d__2 = a[(a_dim1 << 1)
		     + 1], abs(d__2)), d__6 = (d__3 = a[a_dim1 + 2], abs(d__3)
		    ) + (d__4 = a[(a_dim1 << 1) + 2], abs(d__4));
#line 231 "MB03YT.f"
	    h1 = max(d__5,d__6);
/* Computing MAX */
#line 233 "MB03YT.f"
	    d__5 = (d__1 = b[b_dim1 + 1], abs(d__1)) + (d__2 = b[(b_dim1 << 1)
		     + 1], abs(d__2)), d__6 = (d__3 = b[b_dim1 + 2], abs(d__3)
		    ) + (d__4 = b[(b_dim1 << 1) + 2], abs(d__4));
#line 233 "MB03YT.f"
	    h2 = max(d__5,d__6);

#line 236 "MB03YT.f"
	    if (scale1 * h1 >= abs(wr1) * h2) {

/*              Find left rotation matrix Q to zero out B(2,1). */

#line 240 "MB03YT.f"
		dlartg_(&b[b_dim1 + 1], &b[b_dim1 + 2], csl, snl, &r__);

#line 242 "MB03YT.f"
	    } else {

/*              Find left rotation matrix Q to zero out A(2,1). */

#line 246 "MB03YT.f"
		dlartg_(&a[a_dim1 + 1], &a[a_dim1 + 2], csl, snl, &r__);

#line 248 "MB03YT.f"
	    }

#line 250 "MB03YT.f"
	    drot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
#line 251 "MB03YT.f"
	    drot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csl, snl);

#line 253 "MB03YT.f"
	    a[a_dim1 + 2] = 0.;
#line 254 "MB03YT.f"
	    b[b_dim1 + 2] = 0.;

/*           Re-adjoint B. */

#line 258 "MB03YT.f"
	    r__ = b[b_dim1 + 1];
#line 259 "MB03YT.f"
	    b[b_dim1 + 1] = b[(b_dim1 << 1) + 2];
#line 260 "MB03YT.f"
	    b[(b_dim1 << 1) + 2] = r__;
#line 261 "MB03YT.f"
	    b[(b_dim1 << 1) + 1] = -b[(b_dim1 << 1) + 1];

#line 263 "MB03YT.f"
	} else {

/*           A pair of complex conjugate eigenvalues: */
/*           first compute the SVD of the matrix adj(B). */

#line 268 "MB03YT.f"
	    r__ = b[b_dim1 + 1];
#line 269 "MB03YT.f"
	    b[b_dim1 + 1] = b[(b_dim1 << 1) + 2];
#line 270 "MB03YT.f"
	    b[(b_dim1 << 1) + 2] = r__;
#line 271 "MB03YT.f"
	    b[(b_dim1 << 1) + 1] = -b[(b_dim1 << 1) + 1];
#line 272 "MB03YT.f"
	    dlasv2_(&b[b_dim1 + 1], &b[(b_dim1 << 1) + 1], &b[(b_dim1 << 1) + 
		    2], &r__, &t, snl, csl, snr, csr);

/*           Form (A,B) := Q(A,adj(B))Z' where Q is left rotation matrix */
/*           and Z is right rotation matrix computed from DLASV2. */

#line 278 "MB03YT.f"
	    drot_(&c__2, &a[a_dim1 + 1], lda, &a[a_dim1 + 2], lda, csl, snl);
#line 279 "MB03YT.f"
	    drot_(&c__2, &b[b_dim1 + 1], ldb, &b[b_dim1 + 2], ldb, csr, snr);
#line 280 "MB03YT.f"
	    drot_(&c__2, &a[a_dim1 + 1], &c__1, &a[(a_dim1 << 1) + 1], &c__1, 
		    csr, snr);
#line 281 "MB03YT.f"
	    drot_(&c__2, &b[b_dim1 + 1], &c__1, &b[(b_dim1 << 1) + 1], &c__1, 
		    csl, snl);

#line 283 "MB03YT.f"
	    b[b_dim1 + 2] = 0.;
#line 284 "MB03YT.f"
	    b[(b_dim1 << 1) + 1] = 0.;
#line 285 "MB03YT.f"
	}

#line 287 "MB03YT.f"
    }

/*     Unscaling */

#line 291 "MB03YT.f"
    r__ = b[b_dim1 + 1];
#line 292 "MB03YT.f"
    t = b[(b_dim1 << 1) + 2];
#line 293 "MB03YT.f"
    a[a_dim1 + 1] = anorm * a[a_dim1 + 1];
#line 294 "MB03YT.f"
    a[a_dim1 + 2] = anorm * a[a_dim1 + 2];
#line 295 "MB03YT.f"
    a[(a_dim1 << 1) + 1] = anorm * a[(a_dim1 << 1) + 1];
#line 296 "MB03YT.f"
    a[(a_dim1 << 1) + 2] = anorm * a[(a_dim1 << 1) + 2];
#line 297 "MB03YT.f"
    b[b_dim1 + 1] = bnorm * b[b_dim1 + 1];
#line 298 "MB03YT.f"
    b[b_dim1 + 2] = bnorm * b[b_dim1 + 2];
#line 299 "MB03YT.f"
    b[(b_dim1 << 1) + 1] = bnorm * b[(b_dim1 << 1) + 1];
#line 300 "MB03YT.f"
    b[(b_dim1 << 1) + 2] = bnorm * b[(b_dim1 << 1) + 2];

#line 302 "MB03YT.f"
    if (wi == 0.) {
#line 303 "MB03YT.f"
	alphar[1] = a[a_dim1 + 1];
#line 304 "MB03YT.f"
	alphar[2] = a[(a_dim1 << 1) + 2];
#line 305 "MB03YT.f"
	alphai[1] = 0.;
#line 306 "MB03YT.f"
	alphai[2] = 0.;
#line 307 "MB03YT.f"
	beta[1] = b[b_dim1 + 1];
#line 308 "MB03YT.f"
	beta[2] = b[(b_dim1 << 1) + 2];
#line 309 "MB03YT.f"
    } else {
#line 310 "MB03YT.f"
	wr1 = anorm * wr1;
#line 311 "MB03YT.f"
	wi = anorm * wi;
#line 312 "MB03YT.f"
	if (abs(wr1) > 1. || wi > 1.) {
#line 313 "MB03YT.f"
	    wr1 *= r__;
#line 314 "MB03YT.f"
	    wi *= r__;
#line 315 "MB03YT.f"
	    r__ = 1.;
#line 316 "MB03YT.f"
	}
#line 317 "MB03YT.f"
	if (abs(wr1) > 1. || abs(wi) > 1.) {
#line 318 "MB03YT.f"
	    wr1 *= t;
#line 319 "MB03YT.f"
	    wi *= t;
#line 320 "MB03YT.f"
	    t = 1.;
#line 321 "MB03YT.f"
	}
#line 322 "MB03YT.f"
	alphar[1] = wr1 / scale1 * r__ * t;
#line 323 "MB03YT.f"
	alphai[1] = (d__1 = wi / scale1 * r__ * t, abs(d__1));
#line 324 "MB03YT.f"
	alphar[2] = alphar[1];
#line 325 "MB03YT.f"
	alphai[2] = -alphai[1];
#line 326 "MB03YT.f"
	beta[1] = bnorm;
#line 327 "MB03YT.f"
	beta[2] = bnorm;
#line 328 "MB03YT.f"
    }
#line 329 "MB03YT.f"
    return 0;
/* *** Last line of MB03YT *** */
} /* mb03yt_ */

