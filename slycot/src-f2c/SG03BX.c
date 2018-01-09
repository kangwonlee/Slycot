#line 1 "SG03BX.f"
/* SG03BX.f -- translated by f2c (version 20100827).
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

#line 1 "SG03BX.f"
/* Table of constant values */

static integer c__2 = 2;
static doublereal c_b8 = 0.;
static doublereal c_b14 = 1.;
static doublereal c_b60 = -1.;
static integer c__4 = 4;
static integer c__1 = 1;

/* Subroutine */ int sg03bx_(char *dico, char *trans, doublereal *a, integer *
	lda, doublereal *e, integer *lde, doublereal *b, integer *ldb, 
	doublereal *u, integer *ldu, doublereal *scale, doublereal *m1, 
	integer *ldm1, doublereal *m2, integer *ldm2, integer *info, ftnlen 
	dico_len, ftnlen trans_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, e_dim1, e_offset, m1_dim1, 
	    m1_offset, m2_dim1, m2_offset, u_dim1, u_offset;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal l, t, v, w, aa[4]	/* was [2][2] */, b11, bb[4]	/* 
	    was [2][2] */, b22, ai[4]	/* was [2][2] */, bi[4]	/* was [2][2] 
	    */, ci, ee[4]	/* was [2][2] */, ei[4]	/* was [2][2] */, ar[
	    4]	/* was [2][2] */, br[4]	/* was [2][2] */, cr, er[4]	/* 
	    was [2][2] */, qi[4]	/* was [2][2] */, si, ti[4]	/* 
	    was [2][2] */, ui[4]	/* was [2][2] */, xi, yi, qr[4]	/* 
	    was [2][2] */, zi[4]	/* was [2][2] */, sr, tr[4]	/* 
	    was [2][2] */, ur[4]	/* was [2][2] */, xr, yr, zr[4]	/* 
	    was [2][2] */, m1i[4]	/* was [2][2] */, m2i[4]	/* 
	    was [2][2] */, m1r[4]	/* was [2][2] */, m2r[4]	/* 
	    was [2][2] */, b12i, b12r, qbi[4]	/* was [2][2] */, qbr[4]	
	    /* was [2][2] */, eps, qui[4]	/* was [2][2] */, qur[4]	
	    /* was [2][2] */, lami, lamr;
    extern /* Subroutine */ int dlag2_(doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal betai, alpha;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal betar;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), sg03by_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal scale1, scale2;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dladiv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal bignum;
    static logical iscont;
    static doublereal smlnum;
    static logical istrns;


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

/*     To solve for X = op(U)**T * op(U) either the generalized c-stable */
/*     continuous-time Lyapunov equation */

/*             T                    T */
/*        op(A)  * X * op(E) + op(E)  * X * op(A) */

/*                 2        T */
/*        = - SCALE  * op(B)  * op(B),                                (1) */

/*     or the generalized d-stable discrete-time Lyapunov equation */

/*             T                    T */
/*        op(A)  * X * op(A) - op(E)  * X * op(E) */

/*                 2        T */
/*        = - SCALE  * op(B)  * op(B),                                (2) */

/*     where op(K) is either K or K**T for K = A, B, E, U. The Cholesky */
/*     factor U of the solution is computed without first finding X. */

/*     Furthermore, the auxiliary matrices */

/*                                   -1        -1 */
/*        M1 := op(U) * op(A) * op(E)   * op(U) */

/*                           -1        -1 */
/*        M2 := op(B) * op(E)   * op(U) */

/*     are computed in a numerically reliable way. */

/*     The matrices A, B, E, M1, M2, and U are real 2-by-2 matrices. The */
/*     pencil A - lambda * E must have a pair of complex conjugate */
/*     eigenvalues. The eigenvalues must be in the open right half plane */
/*     (in the continuous-time case) or inside the unit circle (in the */
/*     discrete-time case). */

/*     The resulting matrix U is upper triangular. The entries on its */
/*     main diagonal are non-negative. SCALE is an output scale factor */
/*     set to avoid overflow in U. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies whether the continuous-time or the discrete-time */
/*             equation is to be solved: */
/*             = 'C':  Solve continuous-time equation (1); */
/*             = 'D':  Solve discrete-time equation (2). */

/*     TRANS   CHARACTER*1 */
/*             Specifies whether the transposed equation is to be solved */
/*             or not: */
/*             = 'N':  op(K) = K,     K = A, B, E, U; */
/*             = 'T':  op(K) = K**T,  K = A, B, E, U. */

/*     Input/Output Parameters */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,2) */
/*             The leading 2-by-2 part of this array must contain the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= 2. */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,2) */
/*             The leading 2-by-2 upper triangular part of this array */
/*             must contain the matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E.  LDE >= 2. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,2) */
/*             The leading 2-by-2 upper triangular part of this array */
/*             must contain the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= 2. */

/*     U       (output) DOUBLE PRECISION array, dimension (LDU,2) */
/*             The leading 2-by-2 part of this array contains the upper */
/*             triangular matrix U. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U.  LDU >= 2. */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor set to avoid overflow in U. */
/*             0 < SCALE <= 1. */

/*     M1      (output) DOUBLE PRECISION array, dimension (LDM1,2) */
/*             The leading 2-by-2 part of this array contains the */
/*             matrix M1. */

/*     LDM1    INTEGER */
/*             The leading dimension of the array M1.  LDM1 >= 2. */

/*     M2      (output) DOUBLE PRECISION array, dimension (LDM2,2) */
/*             The leading 2-by-2 part of this array contains the */
/*             matrix M2. */

/*     LDM2    INTEGER */
/*             The leading dimension of the array M2.  LDM2 >= 2. */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 2:  the eigenvalues of the pencil A - lambda * E are not */
/*                   a pair of complex conjugate numbers; */
/*             = 3:  the eigenvalues of the pencil A - lambda * E are */
/*                   not in the open right half plane (in the continuous- */
/*                   time case) or inside the unit circle (in the */
/*                   discrete-time case). */

/*     METHOD */

/*     The method used by the routine is based on a generalization of the */
/*     method due to Hammarling ([1], section 6) for Lyapunov equations */
/*     of order 2. A more detailed description is given in [2]. */

/*     REFERENCES */

/*     [1] Hammarling, S.J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-323, 1982. */

/*     [2] Penzl, T. */
/*         Numerical solution of generalized Lyapunov equations. */
/*         Advances in Comp. Math., vol. 8, pp. 33-48, 1998. */

/*     FURTHER COMMENTS */

/*     If the solution matrix U is singular, the matrices M1 and M2 are */
/*     properly set (see [1], equation (6.21)). */

/*     CONTRIBUTOR */

/*     T. Penzl, Technical University Chemnitz, Germany, Aug. 1998. */

/*     REVISIONS */

/*     Sep. 1998 (V. Sima). */
/*     Dec. 1998 (V. Sima). */
/*     July 2003 (V. Sima; suggested by Klaus Schnepper). */
/*     Oct. 2003 (A. Varga). */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     Decode input parameters. */

#line 209 "SG03BX.f"
    /* Parameter adjustments */
#line 209 "SG03BX.f"
    a_dim1 = *lda;
#line 209 "SG03BX.f"
    a_offset = 1 + a_dim1;
#line 209 "SG03BX.f"
    a -= a_offset;
#line 209 "SG03BX.f"
    e_dim1 = *lde;
#line 209 "SG03BX.f"
    e_offset = 1 + e_dim1;
#line 209 "SG03BX.f"
    e -= e_offset;
#line 209 "SG03BX.f"
    b_dim1 = *ldb;
#line 209 "SG03BX.f"
    b_offset = 1 + b_dim1;
#line 209 "SG03BX.f"
    b -= b_offset;
#line 209 "SG03BX.f"
    u_dim1 = *ldu;
#line 209 "SG03BX.f"
    u_offset = 1 + u_dim1;
#line 209 "SG03BX.f"
    u -= u_offset;
#line 209 "SG03BX.f"
    m1_dim1 = *ldm1;
#line 209 "SG03BX.f"
    m1_offset = 1 + m1_dim1;
#line 209 "SG03BX.f"
    m1 -= m1_offset;
#line 209 "SG03BX.f"
    m2_dim1 = *ldm2;
#line 209 "SG03BX.f"
    m2_offset = 1 + m2_dim1;
#line 209 "SG03BX.f"
    m2 -= m2_offset;
#line 209 "SG03BX.f"

#line 209 "SG03BX.f"
    /* Function Body */
#line 209 "SG03BX.f"
    istrns = lsame_(trans, "T", (ftnlen)1, (ftnlen)1);
#line 210 "SG03BX.f"
    iscont = lsame_(dico, "C", (ftnlen)1, (ftnlen)1);

/*     Do not check input parameters for errors. */

/*     Set constants to control overflow. */

#line 216 "SG03BX.f"
    eps = dlamch_("P", (ftnlen)1);
#line 217 "SG03BX.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 218 "SG03BX.f"
    bignum = 1. / smlnum;
#line 219 "SG03BX.f"
    dlabad_(&smlnum, &bignum);

#line 221 "SG03BX.f"
    *info = 0;
#line 222 "SG03BX.f"
    *scale = 1.;

/*     Make copies of A, E, and B. */

#line 226 "SG03BX.f"
    aa[0] = a[a_dim1 + 1];
#line 227 "SG03BX.f"
    aa[1] = a[a_dim1 + 2];
#line 228 "SG03BX.f"
    aa[2] = a[(a_dim1 << 1) + 1];
#line 229 "SG03BX.f"
    aa[3] = a[(a_dim1 << 1) + 2];
#line 230 "SG03BX.f"
    ee[0] = e[e_dim1 + 1];
#line 231 "SG03BX.f"
    ee[1] = 0.;
#line 232 "SG03BX.f"
    ee[2] = e[(e_dim1 << 1) + 1];
#line 233 "SG03BX.f"
    ee[3] = e[(e_dim1 << 1) + 2];
#line 234 "SG03BX.f"
    bb[0] = b[b_dim1 + 1];
#line 235 "SG03BX.f"
    bb[1] = 0.;
#line 236 "SG03BX.f"
    bb[2] = b[(b_dim1 << 1) + 1];
#line 237 "SG03BX.f"
    bb[3] = b[(b_dim1 << 1) + 2];

/*     If the transposed equation (op(K)=K**T, K=A,B,E,U) is to be */
/*     solved, transpose the matrices A, E, B with respect to the */
/*     anti-diagonal. This results in a non-transposed equation. */

#line 243 "SG03BX.f"
    if (istrns) {
#line 244 "SG03BX.f"
	v = aa[0];
#line 245 "SG03BX.f"
	aa[0] = aa[3];
#line 246 "SG03BX.f"
	aa[3] = v;
#line 247 "SG03BX.f"
	v = ee[0];
#line 248 "SG03BX.f"
	ee[0] = ee[3];
#line 249 "SG03BX.f"
	ee[3] = v;
#line 250 "SG03BX.f"
	v = bb[0];
#line 251 "SG03BX.f"
	bb[0] = bb[3];
#line 252 "SG03BX.f"
	bb[3] = v;
#line 253 "SG03BX.f"
    }

/*     Perform QZ-step to transform the pencil A - lambda * E to */
/*     generalized Schur form. The main diagonal of the Schur factor of E */
/*     is real and positive. */

/*     Compute eigenvalues (LAMR + LAMI * I, LAMR - LAMI * I). */

/* Computing MAX */
/* Computing MAX */
#line 261 "SG03BX.f"
    d__2 = abs(ee[0]), d__3 = abs(ee[2]), d__2 = max(d__2,d__3), d__3 = abs(
	    ee[3]);
#line 261 "SG03BX.f"
    d__1 = eps * max(d__2,d__3);
#line 261 "SG03BX.f"
    t = max(d__1,smlnum);
/* Computing MIN */
#line 263 "SG03BX.f"
    d__1 = abs(ee[0]), d__2 = abs(ee[3]);
#line 263 "SG03BX.f"
    if (min(d__1,d__2) < t) {
#line 264 "SG03BX.f"
	*info = 3;
#line 265 "SG03BX.f"
	return 0;
#line 266 "SG03BX.f"
    }
#line 267 "SG03BX.f"
    d__1 = smlnum * eps;
#line 267 "SG03BX.f"
    dlag2_(aa, &c__2, ee, &c__2, &d__1, &scale1, &scale2, &lamr, &w, &lami);
#line 269 "SG03BX.f"
    if (lami <= 0.) {
#line 270 "SG03BX.f"
	*info = 2;
#line 271 "SG03BX.f"
	return 0;
#line 272 "SG03BX.f"
    }

/*     Compute right orthogonal transformation matrix Q. */

#line 276 "SG03BX.f"
    d__1 = scale1 * aa[0] - ee[0] * lamr;
#line 276 "SG03BX.f"
    d__2 = -ee[0] * lami;
#line 276 "SG03BX.f"
    d__3 = scale1 * aa[1];
#line 276 "SG03BX.f"
    sg03by_(&d__1, &d__2, &d__3, &c_b8, &cr, &ci, &sr, &si, &l);
#line 278 "SG03BX.f"
    qr[0] = cr;
#line 279 "SG03BX.f"
    qr[2] = sr;
#line 280 "SG03BX.f"
    qr[1] = -sr;
#line 281 "SG03BX.f"
    qr[3] = cr;
#line 282 "SG03BX.f"
    qi[0] = -ci;
#line 283 "SG03BX.f"
    qi[2] = -si;
#line 284 "SG03BX.f"
    qi[1] = -si;
#line 285 "SG03BX.f"
    qi[3] = ci;

/*     A := Q * A */

#line 289 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qr, &c__2, aa, &c__2, &c_b8,
	     ar, &c__2, (ftnlen)1, (ftnlen)1);
#line 290 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qi, &c__2, aa, &c__2, &c_b8,
	     ai, &c__2, (ftnlen)1, (ftnlen)1);

/*     E := Q * E */

#line 294 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qr, &c__2, ee, &c__2, &c_b8,
	     er, &c__2, (ftnlen)1, (ftnlen)1);
#line 295 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qi, &c__2, ee, &c__2, &c_b8,
	     ei, &c__2, (ftnlen)1, (ftnlen)1);

/*     Compute left orthogonal transformation matrix Z. */

#line 299 "SG03BX.f"
    sg03by_(&er[3], &ei[3], &er[1], &ei[1], &cr, &ci, &sr, &si, &l);
#line 301 "SG03BX.f"
    zr[0] = cr;
#line 302 "SG03BX.f"
    zr[2] = sr;
#line 303 "SG03BX.f"
    zr[1] = -sr;
#line 304 "SG03BX.f"
    zr[3] = cr;
#line 305 "SG03BX.f"
    zi[0] = ci;
#line 306 "SG03BX.f"
    zi[2] = -si;
#line 307 "SG03BX.f"
    zi[1] = -si;
#line 308 "SG03BX.f"
    zi[3] = -ci;

/*     E := E * Z */

#line 312 "SG03BX.f"
    dgemv_("T", &c__2, &c__2, &c_b14, zr, &c__2, er, &c__2, &c_b8, tr, &c__2, 
	    (ftnlen)1);
#line 313 "SG03BX.f"
    dgemv_("T", &c__2, &c__2, &c_b60, zi, &c__2, ei, &c__2, &c_b14, tr, &c__2,
	     (ftnlen)1);
#line 314 "SG03BX.f"
    dgemv_("T", &c__2, &c__2, &c_b14, zi, &c__2, er, &c__2, &c_b8, ti, &c__2, 
	    (ftnlen)1);
#line 315 "SG03BX.f"
    dgemv_("T", &c__2, &c__2, &c_b14, zr, &c__2, ei, &c__2, &c_b14, ti, &c__2,
	     (ftnlen)1);
#line 316 "SG03BX.f"
    dcopy_(&c__2, tr, &c__2, er, &c__2);
#line 317 "SG03BX.f"
    dcopy_(&c__2, ti, &c__2, ei, &c__2);
#line 318 "SG03BX.f"
    er[1] = 0.;
#line 319 "SG03BX.f"
    er[3] = l;
#line 320 "SG03BX.f"
    ei[1] = 0.;
#line 321 "SG03BX.f"
    ei[3] = 0.;

/*     Make main diagonal entries of E real and positive. */
/*     (Note:  Z and E are altered.) */

#line 326 "SG03BX.f"
    v = dlapy2_(er, ei);
#line 327 "SG03BX.f"
    dladiv_(&v, &c_b8, er, ei, &xr, &xi);
#line 328 "SG03BX.f"
    er[0] = v;
#line 329 "SG03BX.f"
    ei[0] = 0.;
#line 330 "SG03BX.f"
    yr = zr[0];
#line 331 "SG03BX.f"
    yi = zi[0];
#line 332 "SG03BX.f"
    zr[0] = xr * yr - xi * yi;
#line 333 "SG03BX.f"
    zi[0] = xr * yi + xi * yr;
#line 334 "SG03BX.f"
    yr = zr[1];
#line 335 "SG03BX.f"
    yi = zi[1];
#line 336 "SG03BX.f"
    zr[1] = xr * yr - xi * yi;
#line 337 "SG03BX.f"
    zi[1] = xr * yi + xi * yr;

/*     A := A * Z */

#line 341 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ar, &c__2, zr, &c__2, &c_b8,
	     tr, &c__2, (ftnlen)1, (ftnlen)1);
#line 342 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b60, ai, &c__2, zi, &c__2, &
	    c_b14, tr, &c__2, (ftnlen)1, (ftnlen)1);
#line 343 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ar, &c__2, zi, &c__2, &c_b8,
	     ti, &c__2, (ftnlen)1, (ftnlen)1);
#line 344 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ai, &c__2, zr, &c__2, &
	    c_b14, ti, &c__2, (ftnlen)1, (ftnlen)1);
#line 345 "SG03BX.f"
    dcopy_(&c__4, tr, &c__1, ar, &c__1);
#line 346 "SG03BX.f"
    dcopy_(&c__4, ti, &c__1, ai, &c__1);

/*     End of QZ-step. */

/*     B := B * Z */

#line 352 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, bb, &c__2, zr, &c__2, &c_b8,
	     br, &c__2, (ftnlen)1, (ftnlen)1);
#line 353 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, bb, &c__2, zi, &c__2, &c_b8,
	     bi, &c__2, (ftnlen)1, (ftnlen)1);

/*     Overwrite B with the upper triangular matrix of its */
/*     QR-factorization. The elements on the main diagonal are real */
/*     and non-negative. */

#line 359 "SG03BX.f"
    sg03by_(br, bi, &br[1], &bi[1], &cr, &ci, &sr, &si, &l);
#line 361 "SG03BX.f"
    qbr[0] = cr;
#line 362 "SG03BX.f"
    qbr[2] = sr;
#line 363 "SG03BX.f"
    qbr[1] = -sr;
#line 364 "SG03BX.f"
    qbr[3] = cr;
#line 365 "SG03BX.f"
    qbi[0] = -ci;
#line 366 "SG03BX.f"
    qbi[2] = -si;
#line 367 "SG03BX.f"
    qbi[1] = -si;
#line 368 "SG03BX.f"
    qbi[3] = ci;
#line 369 "SG03BX.f"
    dgemv_("N", &c__2, &c__2, &c_b14, qbr, &c__2, &br[2], &c__1, &c_b8, tr, &
	    c__1, (ftnlen)1);
#line 370 "SG03BX.f"
    dgemv_("N", &c__2, &c__2, &c_b60, qbi, &c__2, &bi[2], &c__1, &c_b14, tr, &
	    c__1, (ftnlen)1);
#line 371 "SG03BX.f"
    dgemv_("N", &c__2, &c__2, &c_b14, qbi, &c__2, &br[2], &c__1, &c_b8, ti, &
	    c__1, (ftnlen)1);
#line 372 "SG03BX.f"
    dgemv_("N", &c__2, &c__2, &c_b14, qbr, &c__2, &bi[2], &c__1, &c_b14, ti, &
	    c__1, (ftnlen)1);
#line 373 "SG03BX.f"
    dcopy_(&c__2, tr, &c__1, &br[2], &c__1);
#line 374 "SG03BX.f"
    dcopy_(&c__2, ti, &c__1, &bi[2], &c__1);
#line 375 "SG03BX.f"
    br[0] = l;
#line 376 "SG03BX.f"
    br[1] = 0.;
#line 377 "SG03BX.f"
    bi[0] = 0.;
#line 378 "SG03BX.f"
    bi[1] = 0.;
#line 379 "SG03BX.f"
    v = dlapy2_(&br[3], &bi[3]);
/* Computing MAX */
/* Computing MAX */
#line 380 "SG03BX.f"
    d__2 = br[0], d__3 = dlapy2_(&br[2], &bi[2]);
#line 380 "SG03BX.f"
    d__1 = eps * max(d__2,d__3);
#line 380 "SG03BX.f"
    if (v >= max(d__1,smlnum)) {
#line 382 "SG03BX.f"
	dladiv_(&v, &c_b8, &br[3], &bi[3], &xr, &xi);
#line 383 "SG03BX.f"
	br[3] = v;
#line 384 "SG03BX.f"
	yr = qbr[1];
#line 385 "SG03BX.f"
	yi = qbi[1];
#line 386 "SG03BX.f"
	qbr[1] = xr * yr - xi * yi;
#line 387 "SG03BX.f"
	qbi[1] = xr * yi + xi * yr;
#line 388 "SG03BX.f"
	yr = qbr[3];
#line 389 "SG03BX.f"
	yi = qbi[3];
#line 390 "SG03BX.f"
	qbr[3] = xr * yr - xi * yi;
#line 391 "SG03BX.f"
	qbi[3] = xr * yi + xi * yr;
#line 392 "SG03BX.f"
    } else {
#line 393 "SG03BX.f"
	br[3] = 0.;
#line 394 "SG03BX.f"
    }
#line 395 "SG03BX.f"
    bi[3] = 0.;

/*     Compute the Cholesky factor of the solution of the reduced */
/*     equation. The solution may be scaled to avoid overflow. */

#line 400 "SG03BX.f"
    if (iscont) {

/*        Continuous-time equation. */

/*        Step I:  Compute U(1,1). Set U(2,1) = 0. */

#line 406 "SG03BX.f"
	v = (ar[0] * er[0] + ai[0] * ei[0]) * -2.;
#line 407 "SG03BX.f"
	if (v <= 0.) {
#line 408 "SG03BX.f"
	    *info = 3;
#line 409 "SG03BX.f"
	    return 0;
#line 410 "SG03BX.f"
	}
#line 411 "SG03BX.f"
	v = sqrt(v);
#line 412 "SG03BX.f"
	t = abs(br[0]) * 2. * smlnum;
#line 413 "SG03BX.f"
	if (t > v) {
#line 414 "SG03BX.f"
	    scale1 = v / t;
#line 415 "SG03BX.f"
	    *scale = scale1 * *scale;
#line 416 "SG03BX.f"
	    br[0] = scale1 * br[0];
#line 417 "SG03BX.f"
	    br[2] = scale1 * br[2];
#line 418 "SG03BX.f"
	    bi[2] = scale1 * bi[2];
#line 419 "SG03BX.f"
	    br[3] = scale1 * br[3];
#line 420 "SG03BX.f"
	}
#line 421 "SG03BX.f"
	ur[0] = br[0] / v;
#line 422 "SG03BX.f"
	ui[0] = 0.;
#line 423 "SG03BX.f"
	ur[1] = 0.;
#line 424 "SG03BX.f"
	ui[1] = 0.;

/*        Step II:  Compute U(1,2). */

/* Computing MAX */
/* Computing MAX */
#line 428 "SG03BX.f"
	d__2 = br[3], d__3 = dlapy2_(&br[2], &bi[2]);
#line 428 "SG03BX.f"
	d__1 = eps * max(d__2,d__3);
#line 428 "SG03BX.f"
	t = max(d__1,smlnum);
#line 430 "SG03BX.f"
	if (abs(br[0]) < t) {
#line 431 "SG03BX.f"
	    ur[2] = 0.;
#line 432 "SG03BX.f"
	    ui[2] = 0.;
#line 433 "SG03BX.f"
	} else {
#line 434 "SG03BX.f"
	    xr = ar[0] * er[2] + ai[0] * ei[2];
#line 435 "SG03BX.f"
	    xi = ai[0] * er[2] - ar[0] * ei[2];
#line 436 "SG03BX.f"
	    xr = xr + ar[2] * er[0] + ai[2] * ei[0];
#line 437 "SG03BX.f"
	    xi = xi - ai[2] * er[0] + ar[2] * ei[0];
#line 438 "SG03BX.f"
	    xr = -br[2] * v - xr * ur[0];
#line 439 "SG03BX.f"
	    xi = bi[2] * v - xi * ur[0];
#line 440 "SG03BX.f"
	    yr = ar[3] * er[0] + ai[3] * ei[0];
#line 441 "SG03BX.f"
	    yi = -ai[3] * er[0] + ar[3] * ei[0];
#line 442 "SG03BX.f"
	    yr = yr + er[3] * ar[0] + ei[3] * ai[0];
#line 443 "SG03BX.f"
	    yi = yi - ei[3] * ar[0] + er[3] * ai[0];
#line 444 "SG03BX.f"
	    t = dlapy2_(&xr, &xi) * 2. * smlnum;
#line 445 "SG03BX.f"
	    if (t > dlapy2_(&yr, &yi)) {
#line 446 "SG03BX.f"
		scale1 = dlapy2_(&yr, &yi) / t;
#line 447 "SG03BX.f"
		*scale = scale1 * *scale;
#line 448 "SG03BX.f"
		br[0] = scale1 * br[0];
#line 449 "SG03BX.f"
		br[2] = scale1 * br[2];
#line 450 "SG03BX.f"
		bi[2] = scale1 * bi[2];
#line 451 "SG03BX.f"
		br[3] = scale1 * br[3];
#line 452 "SG03BX.f"
		ur[0] = scale1 * ur[0];
#line 453 "SG03BX.f"
		xr = scale1 * xr;
#line 454 "SG03BX.f"
		xi = scale1 * xi;
#line 455 "SG03BX.f"
	    }
#line 456 "SG03BX.f"
	    dladiv_(&xr, &xi, &yr, &yi, &ur[2], &ui[2]);
#line 457 "SG03BX.f"
	    ui[2] = -ui[2];
#line 458 "SG03BX.f"
	}

/*        Step III:  Compute U(2,2). */

#line 462 "SG03BX.f"
	xr = (er[2] * ur[0] + er[3] * ur[2] - ei[3] * ui[2]) * v;
#line 463 "SG03BX.f"
	xi = (-ei[2] * ur[0] - er[3] * ui[2] - ei[3] * ur[2]) * v;
#line 464 "SG03BX.f"
	t = dlapy2_(&xr, &xi) * 2. * smlnum;
#line 465 "SG03BX.f"
	if (t > dlapy2_(er, ei)) {
#line 466 "SG03BX.f"
	    scale1 = dlapy2_(er, ei) / t;
#line 467 "SG03BX.f"
	    *scale = scale1 * *scale;
#line 468 "SG03BX.f"
	    ur[0] = scale1 * ur[0];
#line 469 "SG03BX.f"
	    ur[2] = scale1 * ur[2];
#line 470 "SG03BX.f"
	    ui[2] = scale1 * ui[2];
#line 471 "SG03BX.f"
	    br[0] = scale1 * br[0];
#line 472 "SG03BX.f"
	    br[2] = scale1 * br[2];
#line 473 "SG03BX.f"
	    bi[2] = scale1 * bi[2];
#line 474 "SG03BX.f"
	    br[3] = scale1 * br[3];
#line 475 "SG03BX.f"
	    xr = scale1 * xr;
#line 476 "SG03BX.f"
	    xi = scale1 * xi;
#line 477 "SG03BX.f"
	}
#line 478 "SG03BX.f"
	d__1 = -ei[0];
#line 478 "SG03BX.f"
	dladiv_(&xr, &xi, er, &d__1, &yr, &yi);
#line 479 "SG03BX.f"
	yr = br[2] - yr;
#line 480 "SG03BX.f"
	yi = -bi[2] - yi;
#line 481 "SG03BX.f"
	v = (ar[3] * er[3] + ai[3] * ei[3]) * -2.;
#line 482 "SG03BX.f"
	if (v <= 0.) {
#line 483 "SG03BX.f"
	    *info = 3;
#line 484 "SG03BX.f"
	    return 0;
#line 485 "SG03BX.f"
	}
#line 486 "SG03BX.f"
	v = sqrt(v);
#line 487 "SG03BX.f"
	d__1 = dlapy2_(&br[3], &bi[3]);
#line 487 "SG03BX.f"
	d__2 = dlapy2_(&yr, &yi);
#line 487 "SG03BX.f"
	w = dlapy2_(&d__1, &d__2);
#line 488 "SG03BX.f"
	t = w * 2. * smlnum;
#line 489 "SG03BX.f"
	if (t > v) {
#line 490 "SG03BX.f"
	    scale1 = v / t;
#line 491 "SG03BX.f"
	    *scale = scale1 * *scale;
#line 492 "SG03BX.f"
	    ur[0] = scale1 * ur[0];
#line 493 "SG03BX.f"
	    ur[2] = scale1 * ur[2];
#line 494 "SG03BX.f"
	    ui[2] = scale1 * ui[2];
#line 495 "SG03BX.f"
	    br[0] = scale1 * br[0];
#line 496 "SG03BX.f"
	    br[2] = scale1 * br[2];
#line 497 "SG03BX.f"
	    bi[2] = scale1 * bi[2];
#line 498 "SG03BX.f"
	    br[3] = scale1 * br[3];
#line 499 "SG03BX.f"
	    w = scale1 * w;
#line 500 "SG03BX.f"
	}
#line 501 "SG03BX.f"
	ur[3] = w / v;
#line 502 "SG03BX.f"
	ui[3] = 0.;

/*        Compute matrices M1 and M2 for the reduced equation. */

#line 506 "SG03BX.f"
	m1r[1] = 0.;
#line 507 "SG03BX.f"
	m1i[1] = 0.;
#line 508 "SG03BX.f"
	m2r[1] = 0.;
#line 509 "SG03BX.f"
	m2i[1] = 0.;
#line 510 "SG03BX.f"
	dladiv_(ar, ai, er, ei, &betar, &betai);
#line 511 "SG03BX.f"
	m1r[0] = betar;
#line 512 "SG03BX.f"
	m1i[0] = betai;
#line 513 "SG03BX.f"
	m1r[3] = betar;
#line 514 "SG03BX.f"
	m1i[3] = -betai;
#line 515 "SG03BX.f"
	alpha = sqrt(betar * -2.);
#line 516 "SG03BX.f"
	m2r[0] = alpha;
#line 517 "SG03BX.f"
	m2i[0] = 0.;
#line 518 "SG03BX.f"
	v = er[0] * er[3];
#line 519 "SG03BX.f"
	xr = (-br[0] * er[2] + er[0] * br[2]) / v;
#line 520 "SG03BX.f"
	xi = (-br[0] * ei[2] + er[0] * bi[2]) / v;
#line 521 "SG03BX.f"
	yr = xr - alpha * ur[2];
#line 522 "SG03BX.f"
	yi = -xi + alpha * ui[2];
#line 523 "SG03BX.f"
	if (yr != 0. || yi != 0.) {
#line 524 "SG03BX.f"
	    m2r[2] = yr / ur[3];
#line 525 "SG03BX.f"
	    m2i[2] = -yi / ur[3];
#line 526 "SG03BX.f"
	    m2r[3] = br[3] / (er[3] * ur[3]);
#line 527 "SG03BX.f"
	    m2i[3] = 0.;
#line 528 "SG03BX.f"
	    m1r[2] = -alpha * m2r[2];
#line 529 "SG03BX.f"
	    m1i[2] = -alpha * m2i[2];
#line 530 "SG03BX.f"
	} else {
#line 531 "SG03BX.f"
	    m2r[2] = 0.;
#line 532 "SG03BX.f"
	    m2i[2] = 0.;
#line 533 "SG03BX.f"
	    m2r[3] = alpha;
#line 534 "SG03BX.f"
	    m2i[3] = 0.;
#line 535 "SG03BX.f"
	    m1r[2] = 0.;
#line 536 "SG03BX.f"
	    m1i[2] = 0.;
#line 537 "SG03BX.f"
	}
#line 538 "SG03BX.f"
    } else {

/*        Discrete-time equation. */

/*        Step I:  Compute U(1,1). Set U(2,1) = 0. */

/* Computing 2nd power */
#line 544 "SG03BX.f"
	d__1 = er[0];
/* Computing 2nd power */
#line 544 "SG03BX.f"
	d__2 = ei[0];
/* Computing 2nd power */
#line 544 "SG03BX.f"
	d__3 = ar[0];
/* Computing 2nd power */
#line 544 "SG03BX.f"
	d__4 = ai[0];
#line 544 "SG03BX.f"
	v = d__1 * d__1 + d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
#line 545 "SG03BX.f"
	if (v <= 0.) {
#line 546 "SG03BX.f"
	    *info = 3;
#line 547 "SG03BX.f"
	    return 0;
#line 548 "SG03BX.f"
	}
#line 549 "SG03BX.f"
	v = sqrt(v);
#line 550 "SG03BX.f"
	t = abs(br[0]) * 2. * smlnum;
#line 551 "SG03BX.f"
	if (t > v) {
#line 552 "SG03BX.f"
	    scale1 = v / t;
#line 553 "SG03BX.f"
	    *scale = scale1 * *scale;
#line 554 "SG03BX.f"
	    br[0] = scale1 * br[0];
#line 555 "SG03BX.f"
	    br[2] = scale1 * br[2];
#line 556 "SG03BX.f"
	    bi[2] = scale1 * bi[2];
#line 557 "SG03BX.f"
	    br[3] = scale1 * br[3];
#line 558 "SG03BX.f"
	}
#line 559 "SG03BX.f"
	ur[0] = br[0] / v;
#line 560 "SG03BX.f"
	ui[0] = 0.;
#line 561 "SG03BX.f"
	ur[1] = 0.;
#line 562 "SG03BX.f"
	ui[1] = 0.;

/*        Step II:  Compute U(1,2). */

/* Computing MAX */
/* Computing MAX */
#line 566 "SG03BX.f"
	d__2 = br[3], d__3 = dlapy2_(&br[2], &bi[2]);
#line 566 "SG03BX.f"
	d__1 = eps * max(d__2,d__3);
#line 566 "SG03BX.f"
	t = max(d__1,smlnum);
#line 568 "SG03BX.f"
	if (abs(br[0]) < t) {
#line 569 "SG03BX.f"
	    ur[2] = 0.;
#line 570 "SG03BX.f"
	    ui[2] = 0.;
#line 571 "SG03BX.f"
	} else {
#line 572 "SG03BX.f"
	    xr = ar[0] * ar[2] + ai[0] * ai[2];
#line 573 "SG03BX.f"
	    xi = ai[0] * ar[2] - ar[0] * ai[2];
#line 574 "SG03BX.f"
	    xr = xr - er[2] * er[0] - ei[2] * ei[0];
#line 575 "SG03BX.f"
	    xi = xi + ei[2] * er[0] - er[2] * ei[0];
#line 576 "SG03BX.f"
	    xr = -br[2] * v - xr * ur[0];
#line 577 "SG03BX.f"
	    xi = bi[2] * v - xi * ur[0];
#line 578 "SG03BX.f"
	    yr = ar[3] * ar[0] + ai[3] * ai[0];
#line 579 "SG03BX.f"
	    yi = -ai[3] * ar[0] + ar[3] * ai[0];
#line 580 "SG03BX.f"
	    yr = yr - er[3] * er[0] - ei[3] * ei[0];
#line 581 "SG03BX.f"
	    yi = yi + ei[3] * er[0] - er[3] * ei[0];
#line 582 "SG03BX.f"
	    t = dlapy2_(&xr, &xi) * 2. * smlnum;
#line 583 "SG03BX.f"
	    if (t > dlapy2_(&yr, &yi)) {
#line 584 "SG03BX.f"
		scale1 = dlapy2_(&yr, &yi) / t;
#line 585 "SG03BX.f"
		*scale = scale1 * *scale;
#line 586 "SG03BX.f"
		br[0] = scale1 * br[0];
#line 587 "SG03BX.f"
		br[2] = scale1 * br[2];
#line 588 "SG03BX.f"
		bi[2] = scale1 * bi[2];
#line 589 "SG03BX.f"
		br[3] = scale1 * br[3];
#line 590 "SG03BX.f"
		ur[0] = scale1 * ur[0];
#line 591 "SG03BX.f"
		xr = scale1 * xr;
#line 592 "SG03BX.f"
		xi = scale1 * xi;
#line 593 "SG03BX.f"
	    }
#line 594 "SG03BX.f"
	    dladiv_(&xr, &xi, &yr, &yi, &ur[2], &ui[2]);
#line 595 "SG03BX.f"
	    ui[2] = -ui[2];
#line 596 "SG03BX.f"
	}

/*        Step III:  Compute U(2,2). */

#line 600 "SG03BX.f"
	xr = er[2] * ur[0] + er[3] * ur[2] - ei[3] * ui[2];
#line 601 "SG03BX.f"
	xi = -ei[2] * ur[0] - er[3] * ui[2] - ei[3] * ur[2];
#line 602 "SG03BX.f"
	yr = ar[2] * ur[0] + ar[3] * ur[2] - ai[3] * ui[2];
#line 603 "SG03BX.f"
	yi = -ai[2] * ur[0] - ar[3] * ui[2] - ai[3] * ur[2];
/* Computing 2nd power */
#line 604 "SG03BX.f"
	d__1 = er[3];
/* Computing 2nd power */
#line 604 "SG03BX.f"
	d__2 = ei[3];
/* Computing 2nd power */
#line 604 "SG03BX.f"
	d__3 = ar[3];
/* Computing 2nd power */
#line 604 "SG03BX.f"
	d__4 = ai[3];
#line 604 "SG03BX.f"
	v = d__1 * d__1 + d__2 * d__2 - d__3 * d__3 - d__4 * d__4;
#line 605 "SG03BX.f"
	if (v <= 0.) {
#line 606 "SG03BX.f"
	    *info = 3;
#line 607 "SG03BX.f"
	    return 0;
#line 608 "SG03BX.f"
	}
#line 609 "SG03BX.f"
	v = sqrt(v);
/* Computing MAX */
#line 610 "SG03BX.f"
	d__1 = abs(br[3]), d__2 = abs(br[2]), d__1 = max(d__1,d__2), d__2 = 
		abs(bi[2]), d__1 = max(d__1,d__2), d__2 = abs(xr), d__1 = max(
		d__1,d__2), d__2 = abs(xi), d__1 = max(d__1,d__2), d__2 = abs(
		yr), d__1 = max(d__1,d__2), d__2 = abs(yi);
#line 610 "SG03BX.f"
	t = max(d__1,d__2);
#line 612 "SG03BX.f"
	if (t <= smlnum) {
#line 612 "SG03BX.f"
	    t = 1.;
#line 612 "SG03BX.f"
	}
/* Computing 2nd power */
#line 613 "SG03BX.f"
	d__1 = br[3] / t;
/* Computing 2nd power */
#line 613 "SG03BX.f"
	d__2 = br[2] / t;
/* Computing 2nd power */
#line 613 "SG03BX.f"
	d__3 = bi[2] / t;
/* Computing 2nd power */
#line 613 "SG03BX.f"
	d__4 = xr / t;
/* Computing 2nd power */
#line 613 "SG03BX.f"
	d__5 = xi / t;
/* Computing 2nd power */
#line 613 "SG03BX.f"
	d__6 = yr / t;
/* Computing 2nd power */
#line 613 "SG03BX.f"
	d__7 = yi / t;
#line 613 "SG03BX.f"
	w = d__1 * d__1 + d__2 * d__2 + d__3 * d__3 - d__4 * d__4 - d__5 * 
		d__5 + d__6 * d__6 + d__7 * d__7;
#line 615 "SG03BX.f"
	if (w < 0.) {
#line 616 "SG03BX.f"
	    *info = 3;
#line 617 "SG03BX.f"
	    return 0;
#line 618 "SG03BX.f"
	}
#line 619 "SG03BX.f"
	w = t * sqrt(w);
#line 620 "SG03BX.f"
	t = w * 2. * smlnum;
#line 621 "SG03BX.f"
	if (t > v) {
#line 622 "SG03BX.f"
	    scale1 = v / t;
#line 623 "SG03BX.f"
	    *scale = scale1 * *scale;
#line 624 "SG03BX.f"
	    ur[0] = scale1 * ur[0];
#line 625 "SG03BX.f"
	    ur[2] = scale1 * ur[2];
#line 626 "SG03BX.f"
	    ui[2] = scale1 * ui[2];
#line 627 "SG03BX.f"
	    br[0] = scale1 * br[0];
#line 628 "SG03BX.f"
	    br[2] = scale1 * br[2];
#line 629 "SG03BX.f"
	    bi[2] = scale1 * bi[2];
#line 630 "SG03BX.f"
	    br[3] = scale1 * br[3];
#line 631 "SG03BX.f"
	    w = scale1 * w;
#line 632 "SG03BX.f"
	}
#line 633 "SG03BX.f"
	ur[3] = w / v;
#line 634 "SG03BX.f"
	ui[3] = 0.;

/*        Compute matrices M1 and M2 for the reduced equation. */

#line 638 "SG03BX.f"
	b11 = br[0] / er[0];
#line 639 "SG03BX.f"
	t = er[0] * er[3];
#line 640 "SG03BX.f"
	b12r = (er[0] * br[2] - br[0] * er[2]) / t;
#line 641 "SG03BX.f"
	b12i = (er[0] * bi[2] - br[0] * ei[2]) / t;
#line 642 "SG03BX.f"
	b22 = br[3] / er[3];
#line 643 "SG03BX.f"
	m1r[1] = 0.;
#line 644 "SG03BX.f"
	m1i[1] = 0.;
#line 645 "SG03BX.f"
	m2r[1] = 0.;
#line 646 "SG03BX.f"
	m2i[1] = 0.;
#line 647 "SG03BX.f"
	dladiv_(ar, ai, er, ei, &betar, &betai);
#line 648 "SG03BX.f"
	m1r[0] = betar;
#line 649 "SG03BX.f"
	m1i[0] = betai;
#line 650 "SG03BX.f"
	m1r[3] = betar;
#line 651 "SG03BX.f"
	m1i[3] = -betai;
#line 652 "SG03BX.f"
	v = dlapy2_(&betar, &betai);
#line 653 "SG03BX.f"
	alpha = sqrt((1. - v) * (v + 1.));
#line 654 "SG03BX.f"
	m2r[0] = alpha;
#line 655 "SG03BX.f"
	m2i[0] = 0.;
#line 656 "SG03BX.f"
	xr = (ai[0] * ei[2] - ar[0] * er[2]) / t + ar[2] / er[3];
#line 657 "SG03BX.f"
	xi = (ar[0] * ei[2] + ai[0] * er[2]) / t - ai[2] / er[3];
#line 658 "SG03BX.f"
	xr = betai * -2. * b12i - b11 * xr;
#line 659 "SG03BX.f"
	xi = betai * -2. * b12r - b11 * xi;
#line 660 "SG03BX.f"
	v = (betai - betar) * (betai + betar) + 1.;
#line 661 "SG03BX.f"
	w = betai * -2. * betar;
#line 662 "SG03BX.f"
	dladiv_(&xr, &xi, &v, &w, &yr, &yi);
#line 663 "SG03BX.f"
	if (yr != 0. || yi != 0.) {
#line 664 "SG03BX.f"
	    m2r[2] = (yr * betar - yi * betai) / ur[3];
#line 665 "SG03BX.f"
	    m2i[2] = -(yi * betar + yr * betai) / ur[3];
#line 666 "SG03BX.f"
	    m2r[3] = b22 / ur[3];
#line 667 "SG03BX.f"
	    m2i[3] = 0.;
#line 668 "SG03BX.f"
	    m1r[2] = -alpha * yr / ur[3];
#line 669 "SG03BX.f"
	    m1i[2] = alpha * yi / ur[3];
#line 670 "SG03BX.f"
	} else {
#line 671 "SG03BX.f"
	    m2r[2] = 0.;
#line 672 "SG03BX.f"
	    m2i[2] = 0.;
#line 673 "SG03BX.f"
	    m2r[3] = alpha;
#line 674 "SG03BX.f"
	    m2i[3] = 0.;
#line 675 "SG03BX.f"
	    m1r[2] = 0.;
#line 676 "SG03BX.f"
	    m1i[2] = 0.;
#line 677 "SG03BX.f"
	}
#line 678 "SG03BX.f"
    }

/*     Transform U back:  U := U * Q. */
/*     (Note:  Z is used as workspace.) */

#line 683 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ur, &c__2, qr, &c__2, &c_b8,
	     zr, &c__2, (ftnlen)1, (ftnlen)1);
#line 684 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b60, ui, &c__2, qi, &c__2, &
	    c_b14, zr, &c__2, (ftnlen)1, (ftnlen)1);
#line 685 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ur, &c__2, qi, &c__2, &c_b8,
	     zi, &c__2, (ftnlen)1, (ftnlen)1);
#line 686 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, ui, &c__2, qr, &c__2, &
	    c_b14, zi, &c__2, (ftnlen)1, (ftnlen)1);

/*     Overwrite U with the upper triangular matrix of its */
/*     QR-factorization. The elements on the main diagonal are real */
/*     and non-negative. */

#line 692 "SG03BX.f"
    sg03by_(zr, zi, &zr[1], &zi[1], &cr, &ci, &sr, &si, &l);
#line 694 "SG03BX.f"
    qur[0] = cr;
#line 695 "SG03BX.f"
    qur[2] = sr;
#line 696 "SG03BX.f"
    qur[1] = -sr;
#line 697 "SG03BX.f"
    qur[3] = cr;
#line 698 "SG03BX.f"
    qui[0] = -ci;
#line 699 "SG03BX.f"
    qui[2] = -si;
#line 700 "SG03BX.f"
    qui[1] = -si;
#line 701 "SG03BX.f"
    qui[3] = ci;
#line 702 "SG03BX.f"
    dgemv_("N", &c__2, &c__2, &c_b14, qur, &c__2, &zr[2], &c__1, &c_b8, &u[(
	    u_dim1 << 1) + 1], &c__1, (ftnlen)1);
#line 703 "SG03BX.f"
    dgemv_("N", &c__2, &c__2, &c_b60, qui, &c__2, &zi[2], &c__1, &c_b14, &u[(
	    u_dim1 << 1) + 1], &c__1, (ftnlen)1);
#line 704 "SG03BX.f"
    dgemv_("N", &c__2, &c__2, &c_b14, qui, &c__2, &zr[2], &c__1, &c_b8, &ui[2]
	    , &c__1, (ftnlen)1);
#line 705 "SG03BX.f"
    dgemv_("N", &c__2, &c__2, &c_b14, qur, &c__2, &zi[2], &c__1, &c_b14, &ui[
	    2], &c__1, (ftnlen)1);
#line 706 "SG03BX.f"
    u[u_dim1 + 1] = l;
#line 707 "SG03BX.f"
    u[u_dim1 + 2] = 0.;
#line 708 "SG03BX.f"
    v = dlapy2_(&u[(u_dim1 << 1) + 2], &ui[3]);
#line 709 "SG03BX.f"
    if (v != 0.) {
#line 710 "SG03BX.f"
	dladiv_(&v, &c_b8, &u[(u_dim1 << 1) + 2], &ui[3], &xr, &xi);
#line 711 "SG03BX.f"
	yr = qur[1];
#line 712 "SG03BX.f"
	yi = qui[1];
#line 713 "SG03BX.f"
	qur[1] = xr * yr - xi * yi;
#line 714 "SG03BX.f"
	qui[1] = xr * yi + xi * yr;
#line 715 "SG03BX.f"
	yr = qur[3];
#line 716 "SG03BX.f"
	yi = qui[3];
#line 717 "SG03BX.f"
	qur[3] = xr * yr - xi * yi;
#line 718 "SG03BX.f"
	qui[3] = xr * yi + xi * yr;
#line 719 "SG03BX.f"
    }
#line 720 "SG03BX.f"
    u[(u_dim1 << 1) + 2] = v;

/*     Transform the matrices M1 and M2 back. */

/*        M1 := QU * M1 * QU**H */
/*        M2 := QB**H * M2 * QU**H */

#line 727 "SG03BX.f"
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m1r, &c__2, qur, &c__2, &
	    c_b8, tr, &c__2, (ftnlen)1, (ftnlen)1);
#line 728 "SG03BX.f"
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m1i, &c__2, qui, &c__2, &
	    c_b14, tr, &c__2, (ftnlen)1, (ftnlen)1);
#line 729 "SG03BX.f"
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b60, m1r, &c__2, qui, &c__2, &
	    c_b8, ti, &c__2, (ftnlen)1, (ftnlen)1);
#line 730 "SG03BX.f"
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m1i, &c__2, qur, &c__2, &
	    c_b14, ti, &c__2, (ftnlen)1, (ftnlen)1);
#line 731 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b14, qur, &c__2, tr, &c__2, &
	    c_b8, &m1[m1_offset], ldm1, (ftnlen)1, (ftnlen)1);
#line 733 "SG03BX.f"
    dgemm_("N", "N", &c__2, &c__2, &c__2, &c_b60, qui, &c__2, ti, &c__2, &
	    c_b14, &m1[m1_offset], ldm1, (ftnlen)1, (ftnlen)1);

#line 736 "SG03BX.f"
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m2r, &c__2, qur, &c__2, &
	    c_b8, tr, &c__2, (ftnlen)1, (ftnlen)1);
#line 737 "SG03BX.f"
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m2i, &c__2, qui, &c__2, &
	    c_b14, tr, &c__2, (ftnlen)1, (ftnlen)1);
#line 738 "SG03BX.f"
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b60, m2r, &c__2, qui, &c__2, &
	    c_b8, ti, &c__2, (ftnlen)1, (ftnlen)1);
#line 739 "SG03BX.f"
    dgemm_("N", "T", &c__2, &c__2, &c__2, &c_b14, m2i, &c__2, qur, &c__2, &
	    c_b14, ti, &c__2, (ftnlen)1, (ftnlen)1);
#line 740 "SG03BX.f"
    dgemm_("T", "N", &c__2, &c__2, &c__2, &c_b14, qbr, &c__2, tr, &c__2, &
	    c_b8, &m2[m2_offset], ldm2, (ftnlen)1, (ftnlen)1);
#line 742 "SG03BX.f"
    dgemm_("T", "N", &c__2, &c__2, &c__2, &c_b14, qbi, &c__2, ti, &c__2, &
	    c_b14, &m2[m2_offset], ldm2, (ftnlen)1, (ftnlen)1);

/*     If the transposed equation (op(K)=K**T, K=A,B,E,U) is to be */
/*     solved, transpose the matrix U with respect to the */
/*     anti-diagonal and the matrices M1, M2 with respect to the diagonal */
/*     and the anti-diagonal. */

#line 750 "SG03BX.f"
    if (istrns) {
#line 751 "SG03BX.f"
	v = u[u_dim1 + 1];
#line 752 "SG03BX.f"
	u[u_dim1 + 1] = u[(u_dim1 << 1) + 2];
#line 753 "SG03BX.f"
	u[(u_dim1 << 1) + 2] = v;
#line 754 "SG03BX.f"
	v = m1[m1_dim1 + 1];
#line 755 "SG03BX.f"
	m1[m1_dim1 + 1] = m1[(m1_dim1 << 1) + 2];
#line 756 "SG03BX.f"
	m1[(m1_dim1 << 1) + 2] = v;
#line 757 "SG03BX.f"
	v = m2[m2_dim1 + 1];
#line 758 "SG03BX.f"
	m2[m2_dim1 + 1] = m2[(m2_dim1 << 1) + 2];
#line 759 "SG03BX.f"
	m2[(m2_dim1 << 1) + 2] = v;
#line 760 "SG03BX.f"
    }

#line 762 "SG03BX.f"
    return 0;
/* *** Last line of SG03BX *** */
} /* sg03bx_ */

