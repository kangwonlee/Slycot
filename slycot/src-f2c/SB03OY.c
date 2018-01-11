#line 1 "SB03OY.f"
/* SB03OY.f -- translated by f2c (version 20100827).
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

#line 1 "SB03OY.f"
/* Table of constant values */

static doublereal c_b4 = 1.;

/* Subroutine */ int sb03oy_(logical *discr, logical *ltrans, integer *isgn, 
	doublereal *s, integer *lds, doublereal *r__, integer *ldr, 
	doublereal *a, integer *lda, doublereal *scale, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, r_dim1, r_offset, s_dim1, s_offset;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), sqrt(doublereal);

    /* Local variables */
    static doublereal g[2], t[2], y[2], e1, e2, p1, p2[2], p3, v1, v2[2], v3, 
	    dp[2], s11, s12, s21, s22, dt[2], x11[2], x12[2], x21[2], x22[2], 
	    p3i, p3r, eta, csp[2], csq[2], eps, sgn, cst[2], snp, snq, snt, 
	    absb, absg, abst, temp[2], smin, gamma[2], alpha, delta[2];
    extern /* Subroutine */ int sb03ov_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal tempi, tempr;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    extern doublereal dlapy2_(doublereal *, doublereal *), dlapy3_(doublereal 
	    *, doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scaloc, bignum, smlnum;


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

/*     To solve for the Cholesky factor  U  of  X, */

/*        op(U)'*op(U) = X, */

/*     where  U  is a two-by-two upper triangular matrix, either the */
/*     continuous-time two-by-two Lyapunov equation */
/*                                         2 */
/*         op(S)'*X + X*op(S) = -ISGN*scale *op(R)'*op(R), */

/*     when DISCR = .FALSE., or the discrete-time two-by-two Lyapunov */
/*     equation */
/*                                         2 */
/*         op(S)'*X*op(S) - X = -ISGN*scale *op(R)'*op(R), */

/*     when DISCR = .TRUE., where op(K) = K or K' (i.e., the transpose of */
/*     the matrix K),  S  is a two-by-two matrix with complex conjugate */
/*     eigenvalues,  R  is a two-by-two upper triangular matrix, */
/*     ISGN = -1 or 1,  and  scale  is an output scale factor, set less */
/*     than or equal to 1 to avoid overflow in  X.  The routine also */
/*     computes two matrices, B and A, so that */
/*                                   2 */
/*        B*U = U*S  and  A*U = scale *R,  if  LTRANS = .FALSE.,  or */
/*                                   2 */
/*        U*B = S*U  and  U*A = scale *R,  if  LTRANS = .TRUE., */
/*     which are used by the general Lyapunov solver. */
/*     In the continuous-time case  ISGN*S  must be stable, so that its */
/*     eigenvalues must have strictly negative real parts. */
/*     In the discrete-time case  S  must be convergent if ISGN = 1, that */
/*     is, its eigenvalues must have moduli less than unity, or  S  must */
/*     be completely divergent if ISGN = -1, that is, its eigenvalues */
/*     must have moduli greater than unity. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DISCR   LOGICAL */
/*             Specifies the equation to be solved:       2 */
/*             = .FALSE.: op(S)'*X + X*op(S) = -ISGN*scale *op(R)'*op(R); */
/*                                                        2 */
/*             = .TRUE. : op(S)'*X*op(S) - X = -ISGN*scale *op(R)'*op(R). */

/*     LTRANS  LOGICAL */
/*             Specifies the form of op(K) to be used, as follows: */
/*             = .FALSE.:  op(K) = K    (No transpose); */
/*             = .TRUE. :  op(K) = K**T (Transpose). */

/*     ISGN    INTEGER */
/*             Specifies the sign of the equation as described before. */
/*             ISGN may only be 1 or -1. */

/*     Input/Output Parameters */

/*     S       (input/output) DOUBLE PRECISION array, dimension (LDS,2) */
/*             On entry, S must contain a 2-by-2 matrix. */
/*             On exit, S contains a 2-by-2 matrix B such that B*U = U*S, */
/*             if LTRANS = .FALSE., or U*B = S*U, if LTRANS = .TRUE.. */
/*             Notice that if U is nonsingular then */
/*               B = U*S*inv( U ),  if LTRANS = .FALSE. */
/*               B = inv( U )*S*U,  if LTRANS = .TRUE.. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= 2. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,2) */
/*             On entry, R must contain a 2-by-2 upper triangular matrix. */
/*             The element R( 2, 1 ) is not referenced. */
/*             On exit, R contains U, the 2-by-2 upper triangular */
/*             Cholesky factor of the solution X, X = op(U)'*op(U). */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= 2. */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,2) */
/*             A contains a 2-by-2 upper triangular matrix A satisfying */
/*             A*U/scale = scale*R, if LTRANS = .FALSE., or */
/*             U*A/scale = scale*R, if LTRANS = .TRUE.. */
/*             Notice that if U is nonsingular then */
/*               A = scale*scale*R*inv( U ),  if LTRANS = .FALSE. */
/*               A = scale*scale*inv( U )*R,  if LTRANS = .TRUE.. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= 2. */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if the Lyapunov equation is (nearly) singular */
/*                   (warning indicator); */
/*                   if DISCR = .FALSE., this means that while the */
/*                   matrix S has computed eigenvalues with negative real */
/*                   parts, it is only just stable in the sense that */
/*                   small perturbations in S can make one or more of the */
/*                   eigenvalues have a non-negative real part; */
/*                   if DISCR = .TRUE., this means that while the */
/*                   matrix S has computed eigenvalues inside the unit */
/*                   circle, it is nevertheless only just convergent, in */
/*                   the sense that small perturbations in S can make one */
/*                   or more of the eigenvalues lie outside the unit */
/*                   circle; */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrix S is unchanged); */
/*             = 2:  if DISCR = .FALSE., and ISGN*S is not stable or */
/*                   if DISCR = .TRUE., ISGN = 1 and S is not convergent */
/*                   or if DISCR = .TRUE., ISGN = -1 and S is not */
/*                   completely divergent; */
/*             = 4:  if S has real eigenvalues. */

/*     NOTE: In the interests of speed, this routine does not check all */
/*           inputs for errors. */

/*     METHOD */

/*     The LAPACK scheme for solving 2-by-2 Sylvester equations is */
/*     adapted for 2-by-2 Lyapunov equations, but directly computing the */
/*     Cholesky factor of the solution. */

/*     REFERENCES */

/*     [1] Hammarling S. J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-325, 1982. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB03CY by Sven Hammarling, */
/*     NAG Ltd., United Kingdom, November 1986. */
/*     Partly based on SB03CY and PLYAP2 by A. Varga, University of */
/*     Bochum, May 1992. */

/*     REVISIONS */

/*     Dec. 1997, April 1998. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

/*     ***************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     The comments in this routine refer to notation and equation */
/*     numbers in sections 6 and 10 of [1]. */

/*     Find the eigenvalue  lambda = E1 - i*E2  of s11. */

#line 208 "SB03OY.f"
    /* Parameter adjustments */
#line 208 "SB03OY.f"
    s_dim1 = *lds;
#line 208 "SB03OY.f"
    s_offset = 1 + s_dim1;
#line 208 "SB03OY.f"
    s -= s_offset;
#line 208 "SB03OY.f"
    r_dim1 = *ldr;
#line 208 "SB03OY.f"
    r_offset = 1 + r_dim1;
#line 208 "SB03OY.f"
    r__ -= r_offset;
#line 208 "SB03OY.f"
    a_dim1 = *lda;
#line 208 "SB03OY.f"
    a_offset = 1 + a_dim1;
#line 208 "SB03OY.f"
    a -= a_offset;
#line 208 "SB03OY.f"

#line 208 "SB03OY.f"
    /* Function Body */
#line 208 "SB03OY.f"
    *info = 0;
#line 209 "SB03OY.f"
    sgn = (doublereal) (*isgn);
#line 210 "SB03OY.f"
    s11 = s[s_dim1 + 1];
#line 211 "SB03OY.f"
    s12 = s[(s_dim1 << 1) + 1];
#line 212 "SB03OY.f"
    s21 = s[s_dim1 + 2];
#line 213 "SB03OY.f"
    s22 = s[(s_dim1 << 1) + 2];

/*     Set constants to control overflow. */

#line 217 "SB03OY.f"
    eps = dlamch_("P", (ftnlen)1);
#line 218 "SB03OY.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 219 "SB03OY.f"
    bignum = 1. / smlnum;
#line 220 "SB03OY.f"
    dlabad_(&smlnum, &bignum);
#line 221 "SB03OY.f"
    smlnum = smlnum * 4. / eps;
#line 222 "SB03OY.f"
    bignum = 1. / smlnum;

/* Computing MAX */
/* Computing MAX */
#line 224 "SB03OY.f"
    d__3 = abs(s11), d__4 = abs(s12), d__3 = max(d__3,d__4), d__4 = abs(s21), 
	    d__3 = max(d__3,d__4), d__4 = abs(s22);
#line 224 "SB03OY.f"
    d__1 = smlnum, d__2 = eps * max(d__3,d__4);
#line 224 "SB03OY.f"
    smin = max(d__1,d__2);
#line 226 "SB03OY.f"
    *scale = 1.;

#line 228 "SB03OY.f"
    dlanv2_(&s11, &s12, &s21, &s22, &tempr, &tempi, &e1, &e2, csp, csq);
#line 229 "SB03OY.f"
    if (tempi == 0.) {
#line 230 "SB03OY.f"
	*info = 4;
#line 231 "SB03OY.f"
	return 0;
#line 232 "SB03OY.f"
    }
#line 233 "SB03OY.f"
    absb = dlapy2_(&e1, &e2);
#line 234 "SB03OY.f"
    if (*discr) {
#line 235 "SB03OY.f"
	if (sgn * (absb - 1.) >= 0.) {
#line 236 "SB03OY.f"
	    *info = 2;
#line 237 "SB03OY.f"
	    return 0;
#line 238 "SB03OY.f"
	}
#line 239 "SB03OY.f"
    } else {
#line 240 "SB03OY.f"
	if (sgn * e1 >= 0.) {
#line 241 "SB03OY.f"
	    *info = 2;
#line 242 "SB03OY.f"
	    return 0;
#line 243 "SB03OY.f"
	}
#line 244 "SB03OY.f"
    }

/*     Compute the cos and sine that define  Qhat.  The sine is real. */

#line 248 "SB03OY.f"
    temp[0] = s[s_dim1 + 1] - e1;
#line 249 "SB03OY.f"
    temp[1] = e2;
#line 250 "SB03OY.f"
    if (*ltrans) {
#line 250 "SB03OY.f"
	temp[1] = -e2;
#line 250 "SB03OY.f"
    }
#line 251 "SB03OY.f"
    sb03ov_(temp, &s[s_dim1 + 2], csq, &snq);

/*     beta in (6.9) is given by  beta = E1 + i*E2,  compute  t. */

#line 255 "SB03OY.f"
    temp[0] = csq[0] * s[(s_dim1 << 1) + 1] - snq * s[s_dim1 + 1];
#line 256 "SB03OY.f"
    temp[1] = csq[1] * s[(s_dim1 << 1) + 1];
#line 257 "SB03OY.f"
    tempr = csq[0] * s[(s_dim1 << 1) + 2] - snq * s[s_dim1 + 2];
#line 258 "SB03OY.f"
    tempi = csq[1] * s[(s_dim1 << 1) + 2];
#line 259 "SB03OY.f"
    t[0] = csq[0] * temp[0] - csq[1] * temp[1] + snq * tempr;
#line 260 "SB03OY.f"
    t[1] = csq[0] * temp[1] + csq[1] * temp[0] + snq * tempi;

#line 262 "SB03OY.f"
    if (*ltrans) {
/*                                                         (     -- ) */
/*        Case op(M) = M'.  Note that the modified  R  is  ( p3  p2 ). */
/*                                                         ( 0   p1 ) */

/*        Compute the cos and sine that define  Phat. */

#line 269 "SB03OY.f"
	temp[0] = csq[0] * r__[(r_dim1 << 1) + 2] - snq * r__[(r_dim1 << 1) + 
		1];
#line 270 "SB03OY.f"
	temp[1] = -csq[1] * r__[(r_dim1 << 1) + 2];
#line 271 "SB03OY.f"
	d__1 = -snq * r__[r_dim1 + 1];
#line 271 "SB03OY.f"
	sb03ov_(temp, &d__1, csp, &snp);

/*        Compute p1, p2 and p3 of the relation corresponding to (6.11). */

#line 275 "SB03OY.f"
	p1 = temp[0];
#line 276 "SB03OY.f"
	temp[0] = csq[0] * r__[(r_dim1 << 1) + 1] + snq * r__[(r_dim1 << 1) + 
		2];
#line 277 "SB03OY.f"
	temp[1] = -csq[1] * r__[(r_dim1 << 1) + 1];
#line 278 "SB03OY.f"
	tempr = csq[0] * r__[r_dim1 + 1];
#line 279 "SB03OY.f"
	tempi = -csq[1] * r__[r_dim1 + 1];
#line 280 "SB03OY.f"
	p2[0] = csp[0] * temp[0] - csp[1] * temp[1] + snp * tempr;
#line 281 "SB03OY.f"
	p2[1] = -csp[0] * temp[1] - csp[1] * temp[0] - snp * tempi;
#line 282 "SB03OY.f"
	p3r = csp[0] * tempr + csp[1] * tempi - snp * temp[0];
#line 283 "SB03OY.f"
	p3i = csp[0] * tempi - csp[1] * tempr - snp * temp[1];
#line 284 "SB03OY.f"
    } else {

/*        Case op(M) = M. */

/*        Compute the cos and sine that define  Phat. */

#line 290 "SB03OY.f"
	temp[0] = csq[0] * r__[r_dim1 + 1] + snq * r__[(r_dim1 << 1) + 1];
#line 291 "SB03OY.f"
	temp[1] = csq[1] * r__[r_dim1 + 1];
#line 292 "SB03OY.f"
	d__1 = snq * r__[(r_dim1 << 1) + 2];
#line 292 "SB03OY.f"
	sb03ov_(temp, &d__1, csp, &snp);

/*        Compute p1, p2 and p3 of (6.11). */

#line 296 "SB03OY.f"
	p1 = temp[0];
#line 297 "SB03OY.f"
	temp[0] = csq[0] * r__[(r_dim1 << 1) + 1] - snq * r__[r_dim1 + 1];
#line 298 "SB03OY.f"
	temp[1] = csq[1] * r__[(r_dim1 << 1) + 1];
#line 299 "SB03OY.f"
	tempr = csq[0] * r__[(r_dim1 << 1) + 2];
#line 300 "SB03OY.f"
	tempi = csq[1] * r__[(r_dim1 << 1) + 2];
#line 301 "SB03OY.f"
	p2[0] = csp[0] * temp[0] - csp[1] * temp[1] + snp * tempr;
#line 302 "SB03OY.f"
	p2[1] = csp[0] * temp[1] + csp[1] * temp[0] + snp * tempi;
#line 303 "SB03OY.f"
	p3r = csp[0] * tempr + csp[1] * tempi - snp * temp[0];
#line 304 "SB03OY.f"
	p3i = csp[1] * tempr - csp[0] * tempi + snp * temp[1];
#line 305 "SB03OY.f"
    }

/*     Make  p3  real by multiplying by  conjg ( p3 )/abs( p3 )  to give */

/*     p3 := abs( p3 ). */

#line 311 "SB03OY.f"
    if (p3i == 0.) {
#line 312 "SB03OY.f"
	p3 = abs(p3r);
#line 313 "SB03OY.f"
	dp[0] = d_sign(&c_b4, &p3r);
#line 314 "SB03OY.f"
	dp[1] = 0.;
#line 315 "SB03OY.f"
    } else {
#line 316 "SB03OY.f"
	p3 = dlapy2_(&p3r, &p3i);
#line 317 "SB03OY.f"
	dp[0] = p3r / p3;
#line 318 "SB03OY.f"
	dp[1] = -p3i / p3;
#line 319 "SB03OY.f"
    }

/*     Now compute the quantities v1, v2, v3 and y in (6.13) - (6.15), */
/*     or (10.23) - (10.25). Care is taken to avoid overflows. */

#line 324 "SB03OY.f"
    if (*discr) {
#line 325 "SB03OY.f"
	alpha = sqrt((d__1 = 1. - absb, abs(d__1)) * (absb + 1.));
#line 326 "SB03OY.f"
    } else {
#line 327 "SB03OY.f"
	alpha = sqrt((d__1 = e1 * 2., abs(d__1)));
#line 328 "SB03OY.f"
    }

#line 330 "SB03OY.f"
    scaloc = 1.;
#line 331 "SB03OY.f"
    if (alpha < smin) {
#line 332 "SB03OY.f"
	alpha = smin;
#line 333 "SB03OY.f"
	*info = 1;
#line 334 "SB03OY.f"
    }
#line 335 "SB03OY.f"
    abst = abs(p1);
#line 336 "SB03OY.f"
    if (alpha < 1. && abst > 1.) {
#line 337 "SB03OY.f"
	if (abst > bignum * alpha) {
#line 337 "SB03OY.f"
	    scaloc = 1. / abst;
#line 337 "SB03OY.f"
	}
#line 339 "SB03OY.f"
    }
#line 340 "SB03OY.f"
    if (scaloc != 1.) {
#line 341 "SB03OY.f"
	p1 = scaloc * p1;
#line 342 "SB03OY.f"
	p2[0] = scaloc * p2[0];
#line 343 "SB03OY.f"
	p2[1] = scaloc * p2[1];
#line 344 "SB03OY.f"
	p3 = scaloc * p3;
#line 345 "SB03OY.f"
	*scale = scaloc * *scale;
#line 346 "SB03OY.f"
    }
#line 347 "SB03OY.f"
    v1 = p1 / alpha;

#line 349 "SB03OY.f"
    if (*discr) {
/* Computing 2nd power */
#line 350 "SB03OY.f"
	d__1 = e2;
#line 350 "SB03OY.f"
	g[0] = (1. - e1) * (e1 + 1.) + d__1 * d__1;
#line 351 "SB03OY.f"
	g[1] = e1 * -2. * e2;
#line 352 "SB03OY.f"
	absg = dlapy2_(g, &g[1]);
#line 353 "SB03OY.f"
	scaloc = 1.;
#line 354 "SB03OY.f"
	if (absg < smin) {
#line 355 "SB03OY.f"
	    absg = smin;
#line 356 "SB03OY.f"
	    *info = 1;
#line 357 "SB03OY.f"
	}
#line 358 "SB03OY.f"
	temp[0] = sgn * alpha * p2[0] + v1 * (e1 * t[0] - e2 * t[1]);
#line 359 "SB03OY.f"
	temp[1] = sgn * alpha * p2[1] + v1 * (e1 * t[1] + e2 * t[0]);
/* Computing MAX */
#line 360 "SB03OY.f"
	d__1 = abs(temp[0]), d__2 = abs(temp[1]);
#line 360 "SB03OY.f"
	abst = max(d__1,d__2);
#line 361 "SB03OY.f"
	if (absg < 1. && abst > 1.) {
#line 362 "SB03OY.f"
	    if (abst > bignum * absg) {
#line 362 "SB03OY.f"
		scaloc = 1. / abst;
#line 362 "SB03OY.f"
	    }
#line 364 "SB03OY.f"
	}
#line 365 "SB03OY.f"
	if (scaloc != 1.) {
#line 366 "SB03OY.f"
	    v1 = scaloc * v1;
#line 367 "SB03OY.f"
	    temp[0] = scaloc * temp[0];
#line 368 "SB03OY.f"
	    temp[1] = scaloc * temp[1];
#line 369 "SB03OY.f"
	    p1 = scaloc * p1;
#line 370 "SB03OY.f"
	    p2[0] = scaloc * p2[0];
#line 371 "SB03OY.f"
	    p2[1] = scaloc * p2[1];
#line 372 "SB03OY.f"
	    p3 = scaloc * p3;
#line 373 "SB03OY.f"
	    *scale = scaloc * *scale;
#line 374 "SB03OY.f"
	}
#line 375 "SB03OY.f"
	temp[0] /= absg;
#line 376 "SB03OY.f"
	temp[1] /= absg;

#line 378 "SB03OY.f"
	scaloc = 1.;
#line 379 "SB03OY.f"
	v2[0] = g[0] * temp[0] + g[1] * temp[1];
#line 380 "SB03OY.f"
	v2[1] = g[0] * temp[1] - g[1] * temp[0];
/* Computing MAX */
#line 381 "SB03OY.f"
	d__1 = abs(v2[0]), d__2 = abs(v2[1]);
#line 381 "SB03OY.f"
	abst = max(d__1,d__2);
#line 382 "SB03OY.f"
	if (absg < 1. && abst > 1.) {
#line 383 "SB03OY.f"
	    if (abst > bignum * absg) {
#line 383 "SB03OY.f"
		scaloc = 1. / abst;
#line 383 "SB03OY.f"
	    }
#line 385 "SB03OY.f"
	}
#line 386 "SB03OY.f"
	if (scaloc != 1.) {
#line 387 "SB03OY.f"
	    v1 = scaloc * v1;
#line 388 "SB03OY.f"
	    v2[0] = scaloc * v2[0];
#line 389 "SB03OY.f"
	    v2[1] = scaloc * v2[1];
#line 390 "SB03OY.f"
	    p1 = scaloc * p1;
#line 391 "SB03OY.f"
	    p2[0] = scaloc * p2[0];
#line 392 "SB03OY.f"
	    p2[1] = scaloc * p2[1];
#line 393 "SB03OY.f"
	    p3 = scaloc * p3;
#line 394 "SB03OY.f"
	    *scale = scaloc * *scale;
#line 395 "SB03OY.f"
	}
#line 396 "SB03OY.f"
	v2[0] /= absg;
#line 397 "SB03OY.f"
	v2[1] /= absg;

#line 399 "SB03OY.f"
	scaloc = 1.;
#line 400 "SB03OY.f"
	temp[0] = p1 * t[0] - e2 * 2. * p2[1];
#line 401 "SB03OY.f"
	temp[1] = p1 * t[1] + e2 * 2. * p2[0];
/* Computing MAX */
#line 402 "SB03OY.f"
	d__1 = abs(temp[0]), d__2 = abs(temp[1]);
#line 402 "SB03OY.f"
	abst = max(d__1,d__2);
#line 403 "SB03OY.f"
	if (absg < 1. && abst > 1.) {
#line 404 "SB03OY.f"
	    if (abst > bignum * absg) {
#line 404 "SB03OY.f"
		scaloc = 1. / abst;
#line 404 "SB03OY.f"
	    }
#line 406 "SB03OY.f"
	}
#line 407 "SB03OY.f"
	if (scaloc != 1.) {
#line 408 "SB03OY.f"
	    temp[0] = scaloc * temp[0];
#line 409 "SB03OY.f"
	    temp[1] = scaloc * temp[1];
#line 410 "SB03OY.f"
	    v1 = scaloc * v1;
#line 411 "SB03OY.f"
	    v2[0] = scaloc * v2[0];
#line 412 "SB03OY.f"
	    v2[1] = scaloc * v2[1];
#line 413 "SB03OY.f"
	    p3 = scaloc * p3;
#line 414 "SB03OY.f"
	    *scale = scaloc * *scale;
#line 415 "SB03OY.f"
	}
#line 416 "SB03OY.f"
	temp[0] /= absg;
#line 417 "SB03OY.f"
	temp[1] /= absg;

#line 419 "SB03OY.f"
	scaloc = 1.;
#line 420 "SB03OY.f"
	y[0] = -(g[0] * temp[0] + g[1] * temp[1]);
#line 421 "SB03OY.f"
	y[1] = -(g[0] * temp[1] - g[1] * temp[0]);
/* Computing MAX */
#line 422 "SB03OY.f"
	d__1 = abs(y[0]), d__2 = abs(y[1]);
#line 422 "SB03OY.f"
	abst = max(d__1,d__2);
#line 423 "SB03OY.f"
	if (absg < 1. && abst > 1.) {
#line 424 "SB03OY.f"
	    if (abst > bignum * absg) {
#line 424 "SB03OY.f"
		scaloc = 1. / abst;
#line 424 "SB03OY.f"
	    }
#line 426 "SB03OY.f"
	}
#line 427 "SB03OY.f"
	if (scaloc != 1.) {
#line 428 "SB03OY.f"
	    y[0] = scaloc * y[0];
#line 429 "SB03OY.f"
	    y[1] = scaloc * y[1];
#line 430 "SB03OY.f"
	    v1 = scaloc * v1;
#line 431 "SB03OY.f"
	    v2[0] = scaloc * v2[0];
#line 432 "SB03OY.f"
	    v2[1] = scaloc * v2[1];
#line 433 "SB03OY.f"
	    p3 = scaloc * p3;
#line 434 "SB03OY.f"
	    *scale = scaloc * *scale;
#line 435 "SB03OY.f"
	}
#line 436 "SB03OY.f"
	y[0] /= absg;
#line 437 "SB03OY.f"
	y[1] /= absg;
#line 438 "SB03OY.f"
    } else {

#line 440 "SB03OY.f"
	scaloc = 1.;
#line 441 "SB03OY.f"
	if (absb < smin) {
#line 442 "SB03OY.f"
	    absb = smin;
#line 443 "SB03OY.f"
	    *info = 1;
#line 444 "SB03OY.f"
	}
#line 445 "SB03OY.f"
	temp[0] = sgn * alpha * p2[0] + v1 * t[0];
#line 446 "SB03OY.f"
	temp[1] = sgn * alpha * p2[1] + v1 * t[1];
/* Computing MAX */
#line 447 "SB03OY.f"
	d__1 = abs(temp[0]), d__2 = abs(temp[1]);
#line 447 "SB03OY.f"
	abst = max(d__1,d__2);
#line 448 "SB03OY.f"
	if (absb < 1. && abst > 1.) {
#line 449 "SB03OY.f"
	    if (abst > bignum * absb) {
#line 449 "SB03OY.f"
		scaloc = 1. / abst;
#line 449 "SB03OY.f"
	    }
#line 451 "SB03OY.f"
	}
#line 452 "SB03OY.f"
	if (scaloc != 1.) {
#line 453 "SB03OY.f"
	    v1 = scaloc * v1;
#line 454 "SB03OY.f"
	    temp[0] = scaloc * temp[0];
#line 455 "SB03OY.f"
	    temp[1] = scaloc * temp[1];
#line 456 "SB03OY.f"
	    p2[0] = scaloc * p2[0];
#line 457 "SB03OY.f"
	    p2[1] = scaloc * p2[1];
#line 458 "SB03OY.f"
	    p3 = scaloc * p3;
#line 459 "SB03OY.f"
	    *scale = scaloc * *scale;
#line 460 "SB03OY.f"
	}
#line 461 "SB03OY.f"
	temp[0] /= absb * 2.;
#line 462 "SB03OY.f"
	temp[1] /= absb * 2.;
#line 463 "SB03OY.f"
	scaloc = 1.;
#line 464 "SB03OY.f"
	v2[0] = -(e1 * temp[0] + e2 * temp[1]);
#line 465 "SB03OY.f"
	v2[1] = -(e1 * temp[1] - e2 * temp[0]);
/* Computing MAX */
#line 466 "SB03OY.f"
	d__1 = abs(v2[0]), d__2 = abs(v2[1]);
#line 466 "SB03OY.f"
	abst = max(d__1,d__2);
#line 467 "SB03OY.f"
	if (absb < 1. && abst > 1.) {
#line 468 "SB03OY.f"
	    if (abst > bignum * absb) {
#line 468 "SB03OY.f"
		scaloc = 1. / abst;
#line 468 "SB03OY.f"
	    }
#line 470 "SB03OY.f"
	}
#line 471 "SB03OY.f"
	if (scaloc != 1.) {
#line 472 "SB03OY.f"
	    v1 = scaloc * v1;
#line 473 "SB03OY.f"
	    v2[0] = scaloc * v2[0];
#line 474 "SB03OY.f"
	    v2[1] = scaloc * v2[1];
#line 475 "SB03OY.f"
	    p2[0] = scaloc * p2[0];
#line 476 "SB03OY.f"
	    p2[1] = scaloc * p2[1];
#line 477 "SB03OY.f"
	    p3 = scaloc * p3;
#line 478 "SB03OY.f"
	    *scale = scaloc * *scale;
#line 479 "SB03OY.f"
	}
#line 480 "SB03OY.f"
	v2[0] /= absb;
#line 481 "SB03OY.f"
	v2[1] /= absb;
#line 482 "SB03OY.f"
	y[0] = p2[0] - alpha * v2[0];
#line 483 "SB03OY.f"
	y[1] = p2[1] - alpha * v2[1];
#line 484 "SB03OY.f"
    }

#line 486 "SB03OY.f"
    scaloc = 1.;
#line 487 "SB03OY.f"
    v3 = dlapy3_(&p3, y, &y[1]);
#line 488 "SB03OY.f"
    if (alpha < 1. && v3 > 1.) {
#line 489 "SB03OY.f"
	if (v3 > bignum * alpha) {
#line 489 "SB03OY.f"
	    scaloc = 1. / v3;
#line 489 "SB03OY.f"
	}
#line 491 "SB03OY.f"
    }
#line 492 "SB03OY.f"
    if (scaloc != 1.) {
#line 493 "SB03OY.f"
	v1 = scaloc * v1;
#line 494 "SB03OY.f"
	v2[0] = scaloc * v2[0];
#line 495 "SB03OY.f"
	v2[1] = scaloc * v2[1];
#line 496 "SB03OY.f"
	v3 = scaloc * v3;
#line 497 "SB03OY.f"
	p3 = scaloc * p3;
#line 498 "SB03OY.f"
	*scale = scaloc * *scale;
#line 499 "SB03OY.f"
    }
#line 500 "SB03OY.f"
    v3 /= alpha;

#line 502 "SB03OY.f"
    if (*ltrans) {

/*        Case op(M) = M'. */

/*        Form  X = conjg( Qhat' )*v11. */

#line 508 "SB03OY.f"
	x11[0] = csq[0] * v3;
#line 509 "SB03OY.f"
	x11[1] = csq[1] * v3;
#line 510 "SB03OY.f"
	x21[0] = snq * v3;
#line 511 "SB03OY.f"
	x12[0] = csq[0] * v2[0] + csq[1] * v2[1] - snq * v1;
#line 512 "SB03OY.f"
	x12[1] = -csq[0] * v2[1] + csq[1] * v2[0];
#line 513 "SB03OY.f"
	x22[0] = csq[0] * v1 + snq * v2[0];
#line 514 "SB03OY.f"
	x22[1] = -csq[1] * v1 - snq * v2[1];

/*        Obtain u11 from the RQ-factorization of X. The conjugate of */
/*        X22 should be taken. */

#line 519 "SB03OY.f"
	x22[1] = -x22[1];
#line 520 "SB03OY.f"
	sb03ov_(x22, x21, cst, &snt);
#line 521 "SB03OY.f"
	r__[(r_dim1 << 1) + 2] = x22[0];
#line 522 "SB03OY.f"
	r__[(r_dim1 << 1) + 1] = cst[0] * x12[0] - cst[1] * x12[1] + snt * 
		x11[0];
#line 523 "SB03OY.f"
	tempr = cst[0] * x11[0] + cst[1] * x11[1] - snt * x12[0];
#line 524 "SB03OY.f"
	tempi = cst[0] * x11[1] - cst[1] * x11[0] - snt * x12[1];
#line 525 "SB03OY.f"
	if (tempi == 0.) {
#line 526 "SB03OY.f"
	    r__[r_dim1 + 1] = abs(tempr);
#line 527 "SB03OY.f"
	    dt[0] = d_sign(&c_b4, &tempr);
#line 528 "SB03OY.f"
	    dt[1] = 0.;
#line 529 "SB03OY.f"
	} else {
#line 530 "SB03OY.f"
	    r__[r_dim1 + 1] = dlapy2_(&tempr, &tempi);
#line 531 "SB03OY.f"
	    dt[0] = tempr / r__[r_dim1 + 1];
#line 532 "SB03OY.f"
	    dt[1] = -tempi / r__[r_dim1 + 1];
#line 533 "SB03OY.f"
	}
#line 534 "SB03OY.f"
    } else {

/*        Case op(M) = M. */

/*        Now form  X = v11*conjg( Qhat' ). */

#line 540 "SB03OY.f"
	x11[0] = csq[0] * v1 - snq * v2[0];
#line 541 "SB03OY.f"
	x11[1] = -csq[1] * v1 + snq * v2[1];
#line 542 "SB03OY.f"
	x21[0] = -snq * v3;
#line 543 "SB03OY.f"
	x12[0] = csq[0] * v2[0] + csq[1] * v2[1] + snq * v1;
#line 544 "SB03OY.f"
	x12[1] = -csq[0] * v2[1] + csq[1] * v2[0];
#line 545 "SB03OY.f"
	x22[0] = csq[0] * v3;
#line 546 "SB03OY.f"
	x22[1] = csq[1] * v3;

/*        Obtain u11 from the QR-factorization of X. */

#line 550 "SB03OY.f"
	sb03ov_(x11, x21, cst, &snt);
#line 551 "SB03OY.f"
	r__[r_dim1 + 1] = x11[0];
#line 552 "SB03OY.f"
	r__[(r_dim1 << 1) + 1] = cst[0] * x12[0] + cst[1] * x12[1] + snt * 
		x22[0];
#line 553 "SB03OY.f"
	tempr = cst[0] * x22[0] - cst[1] * x22[1] - snt * x12[0];
#line 554 "SB03OY.f"
	tempi = cst[0] * x22[1] + cst[1] * x22[0] - snt * x12[1];
#line 555 "SB03OY.f"
	if (tempi == 0.) {
#line 556 "SB03OY.f"
	    r__[(r_dim1 << 1) + 2] = abs(tempr);
#line 557 "SB03OY.f"
	    dt[0] = d_sign(&c_b4, &tempr);
#line 558 "SB03OY.f"
	    dt[1] = 0.;
#line 559 "SB03OY.f"
	} else {
#line 560 "SB03OY.f"
	    r__[(r_dim1 << 1) + 2] = dlapy2_(&tempr, &tempi);
#line 561 "SB03OY.f"
	    dt[0] = tempr / r__[(r_dim1 << 1) + 2];
#line 562 "SB03OY.f"
	    dt[1] = -tempi / r__[(r_dim1 << 1) + 2];
#line 563 "SB03OY.f"
	}
#line 564 "SB03OY.f"
    }

/*     The computations below are not needed when B and A are not */
/*     useful. Compute delta, eta and gamma as in (6.21) or (10.26). */

#line 569 "SB03OY.f"
    if (y[0] == 0. && y[1] == 0.) {
#line 570 "SB03OY.f"
	delta[0] = 0.;
#line 571 "SB03OY.f"
	delta[1] = 0.;
#line 572 "SB03OY.f"
	gamma[0] = 0.;
#line 573 "SB03OY.f"
	gamma[1] = 0.;
#line 574 "SB03OY.f"
	eta = alpha;
#line 575 "SB03OY.f"
    } else {
#line 576 "SB03OY.f"
	delta[0] = y[0] / v3;
#line 577 "SB03OY.f"
	delta[1] = y[1] / v3;
#line 578 "SB03OY.f"
	gamma[0] = -alpha * delta[0];
#line 579 "SB03OY.f"
	gamma[1] = -alpha * delta[1];
#line 580 "SB03OY.f"
	eta = p3 / v3;
#line 581 "SB03OY.f"
	if (*discr) {
#line 582 "SB03OY.f"
	    tempr = e1 * delta[0] - e2 * delta[1];
#line 583 "SB03OY.f"
	    delta[1] = e1 * delta[1] + e2 * delta[0];
#line 584 "SB03OY.f"
	    delta[0] = tempr;
#line 585 "SB03OY.f"
	}
#line 586 "SB03OY.f"
    }

#line 588 "SB03OY.f"
    if (*ltrans) {

/*        Case op(M) = M'. */

/*        Find  X = conjg( That' )*( inv( v11 )*s11hat*v11 ). */
/*        ( Defer the scaling.) */

#line 595 "SB03OY.f"
	x11[0] = cst[0] * e1 + cst[1] * e2;
#line 596 "SB03OY.f"
	x11[1] = -cst[0] * e2 + cst[1] * e1;
#line 597 "SB03OY.f"
	x21[0] = snt * e1;
#line 598 "SB03OY.f"
	x21[1] = -snt * e2;
#line 599 "SB03OY.f"
	x12[0] = sgn * (cst[0] * gamma[0] + cst[1] * gamma[1]) - snt * e1;
#line 600 "SB03OY.f"
	x12[1] = sgn * (-cst[0] * gamma[1] + cst[1] * gamma[0]) - snt * e2;
#line 601 "SB03OY.f"
	x22[0] = cst[0] * e1 + cst[1] * e2 + sgn * snt * gamma[0];
#line 602 "SB03OY.f"
	x22[1] = cst[0] * e2 - cst[1] * e1 - sgn * snt * gamma[1];

/*        Now find  B = X*That. ( Include the scaling here.) */

#line 606 "SB03OY.f"
	s[s_dim1 + 1] = cst[0] * x11[0] + cst[1] * x11[1] - snt * x12[0];
#line 607 "SB03OY.f"
	tempr = cst[0] * x21[0] + cst[1] * x21[1] - snt * x22[0];
#line 608 "SB03OY.f"
	tempi = cst[0] * x21[1] - cst[1] * x21[0] - snt * x22[1];
#line 609 "SB03OY.f"
	s[s_dim1 + 2] = dt[0] * tempr - dt[1] * tempi;
#line 610 "SB03OY.f"
	tempr = cst[0] * x12[0] - cst[1] * x12[1] + snt * x11[0];
#line 611 "SB03OY.f"
	tempi = cst[0] * x12[1] + cst[1] * x12[0] + snt * x11[1];
#line 612 "SB03OY.f"
	s[(s_dim1 << 1) + 1] = dt[0] * tempr + dt[1] * tempi;
#line 613 "SB03OY.f"
	s[(s_dim1 << 1) + 2] = cst[0] * x22[0] - cst[1] * x22[1] + snt * x21[
		0];

/*        Form  X = ( inv( v11 )*p11 )*conjg( Phat' ). */

#line 617 "SB03OY.f"
	tempr = dp[0] * eta;
#line 618 "SB03OY.f"
	tempi = -dp[1] * eta;
#line 619 "SB03OY.f"
	x11[0] = csp[0] * tempr - csp[1] * tempi + snp * delta[0];
#line 620 "SB03OY.f"
	x11[1] = csp[0] * tempi + csp[1] * tempr - snp * delta[1];
#line 621 "SB03OY.f"
	x21[0] = snp * alpha;
#line 622 "SB03OY.f"
	x12[0] = -snp * tempr + csp[0] * delta[0] - csp[1] * delta[1];
#line 623 "SB03OY.f"
	x12[1] = -snp * tempi - csp[0] * delta[1] - csp[1] * delta[0];
#line 624 "SB03OY.f"
	x22[0] = csp[0] * alpha;
#line 625 "SB03OY.f"
	x22[1] = -csp[1] * alpha;

/*        Finally form  A = conjg( That' )*X. */

#line 629 "SB03OY.f"
	tempr = cst[0] * x11[0] - cst[1] * x11[1] - snt * x21[0];
#line 630 "SB03OY.f"
	tempi = cst[0] * x22[1] + cst[1] * x22[0];
#line 631 "SB03OY.f"
	a[a_dim1 + 1] = dt[0] * tempr + dt[1] * tempi;
#line 632 "SB03OY.f"
	tempr = cst[0] * x12[0] - cst[1] * x12[1] - snt * x22[0];
#line 633 "SB03OY.f"
	tempi = cst[0] * x12[1] + cst[1] * x12[0] - snt * x22[0];
#line 634 "SB03OY.f"
	a[(a_dim1 << 1) + 1] = dt[0] * tempr + dt[1] * tempi;
#line 635 "SB03OY.f"
	a[a_dim1 + 2] = 0.;
#line 636 "SB03OY.f"
	a[(a_dim1 << 1) + 2] = cst[0] * x22[0] + cst[1] * x22[1] + snt * x12[
		0];
#line 637 "SB03OY.f"
    } else {

/*        Case op(M) = M. */

/*        Find  X = That*( v11*s11hat*inv( v11 ) ). ( Defer the scaling.) */

#line 643 "SB03OY.f"
	x11[0] = cst[0] * e1 + cst[1] * e2;
#line 644 "SB03OY.f"
	x11[1] = cst[0] * e2 - cst[1] * e1;
#line 645 "SB03OY.f"
	x21[0] = -snt * e1;
#line 646 "SB03OY.f"
	x21[1] = -snt * e2;
#line 647 "SB03OY.f"
	x12[0] = sgn * (cst[0] * gamma[0] - cst[1] * gamma[1]) + snt * e1;
#line 648 "SB03OY.f"
	x12[1] = sgn * (-cst[0] * gamma[1] - cst[1] * gamma[0]) - snt * e2;
#line 649 "SB03OY.f"
	x22[0] = cst[0] * e1 + cst[1] * e2 - sgn * snt * gamma[0];
#line 650 "SB03OY.f"
	x22[1] = -cst[0] * e2 + cst[1] * e1 + sgn * snt * gamma[1];

/*        Now find  B = X*conjg( That' ). ( Include the scaling here.) */

#line 654 "SB03OY.f"
	s[s_dim1 + 1] = cst[0] * x11[0] - cst[1] * x11[1] + snt * x12[0];
#line 655 "SB03OY.f"
	tempr = cst[0] * x21[0] - cst[1] * x21[1] + snt * x22[0];
#line 656 "SB03OY.f"
	tempi = cst[0] * x21[1] + cst[1] * x21[0] + snt * x22[1];
#line 657 "SB03OY.f"
	s[s_dim1 + 2] = dt[0] * tempr - dt[1] * tempi;
#line 658 "SB03OY.f"
	tempr = cst[0] * x12[0] + cst[1] * x12[1] - snt * x11[0];
#line 659 "SB03OY.f"
	tempi = cst[0] * x12[1] - cst[1] * x12[0] - snt * x11[1];
#line 660 "SB03OY.f"
	s[(s_dim1 << 1) + 1] = dt[0] * tempr + dt[1] * tempi;
#line 661 "SB03OY.f"
	s[(s_dim1 << 1) + 2] = cst[0] * x22[0] + cst[1] * x22[1] - snt * x21[
		0];

/*        Form  X = Phat*( p11*inv( v11 ) ). */

#line 665 "SB03OY.f"
	tempr = dp[0] * eta;
#line 666 "SB03OY.f"
	tempi = -dp[1] * eta;
#line 667 "SB03OY.f"
	x11[0] = csp[0] * alpha;
#line 668 "SB03OY.f"
	x11[1] = csp[1] * alpha;
#line 669 "SB03OY.f"
	x21[0] = snp * alpha;
#line 670 "SB03OY.f"
	x12[0] = csp[0] * delta[0] + csp[1] * delta[1] - snp * tempr;
#line 671 "SB03OY.f"
	x12[1] = -csp[0] * delta[1] + csp[1] * delta[0] - snp * tempi;
#line 672 "SB03OY.f"
	x22[0] = csp[0] * tempr + csp[1] * tempi + snp * delta[0];
#line 673 "SB03OY.f"
	x22[1] = csp[0] * tempi - csp[1] * tempr - snp * delta[1];

/*        Finally form  A = X*conjg( That' ). */

#line 677 "SB03OY.f"
	a[a_dim1 + 1] = cst[0] * x11[0] - cst[1] * x11[1] + snt * x12[0];
#line 678 "SB03OY.f"
	a[a_dim1 + 2] = 0.;
#line 679 "SB03OY.f"
	a[(a_dim1 << 1) + 1] = cst[0] * x12[0] + cst[1] * x12[1] - snt * x11[
		0];
#line 680 "SB03OY.f"
	tempr = cst[0] * x22[0] + cst[1] * x22[1] - snt * x21[0];
#line 681 "SB03OY.f"
	tempi = cst[0] * x22[1] - cst[1] * x22[0];
#line 682 "SB03OY.f"
	a[(a_dim1 << 1) + 2] = dt[0] * tempr + dt[1] * tempi;
#line 683 "SB03OY.f"
    }

#line 685 "SB03OY.f"
    if (*scale != 1.) {
#line 686 "SB03OY.f"
	a[a_dim1 + 1] = *scale * a[a_dim1 + 1];
#line 687 "SB03OY.f"
	a[(a_dim1 << 1) + 1] = *scale * a[(a_dim1 << 1) + 1];
#line 688 "SB03OY.f"
	a[(a_dim1 << 1) + 2] = *scale * a[(a_dim1 << 1) + 2];
#line 689 "SB03OY.f"
    }

#line 691 "SB03OY.f"
    return 0;
/* *** Last line of SB03OY *** */
} /* sb03oy_ */

