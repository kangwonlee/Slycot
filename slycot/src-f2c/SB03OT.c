#line 1 "SB03OT.f"
/* SB03OT.f -- translated by f2c (version 20100827).
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

#line 1 "SB03OT.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static doublereal c_b19 = -1.;
static doublereal c_b47 = 1.;
static integer c__3 = 3;
static integer c__0 = 0;

/* Subroutine */ int sb03ot_(logical *discr, logical *ltrans, integer *n, 
	doublereal *s, integer *lds, doublereal *r__, integer *ldr, 
	doublereal *scale, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, s_dim1, s_offset, i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static doublereal a[4]	/* was [2][2] */, b[4]	/* was [2][2] */;
    static integer j, k;
    static doublereal u[4]	/* was [2][2] */, d1, d2;
    static integer j1, j2, j3, k1, k2, k3;
    static doublereal t1, t2, t3, t4, v1, v2, v3, v4, dr, eps, sum, tau1, 
	    tau2;
    static integer isgn;
    static logical cont;
    static doublereal temp, smin;
    static logical tbyt;
    extern /* Subroutine */ int mb04nd_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    static doublereal alpha;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), mb04od_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen);
    static integer infom;
    extern /* Subroutine */ int sb03or_(logical *, logical *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *), dcopy_(integer 
	    *, doublereal *, integer *, doublereal *, integer *), dswap_(
	    integer *, doublereal *, integer *, doublereal *, integer *), 
	    sb03oy_(logical *, logical *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), dtrmm_(char *, char *, char *, char *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static integer ksize;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dtrmv_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    static integer kount;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    static doublereal scaloc;
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal absskk, bignum, smlnum;


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

/*     To solve for X = op(U)'*op(U) either the stable non-negative */
/*     definite continuous-time Lyapunov equation */
/*                                   2 */
/*        op(S)'*X + X*op(S) = -scale *op(R)'*op(R)                   (1) */

/*     or the convergent non-negative definite discrete-time Lyapunov */
/*     equation */
/*                                   2 */
/*        op(S)'*X*op(S) - X = -scale *op(R)'*op(R)                   (2) */

/*     where op(K) = K or K' (i.e., the transpose of the matrix K), S is */
/*     an N-by-N block upper triangular matrix with one-by-one or */
/*     two-by-two blocks on the diagonal, R is an N-by-N upper triangular */
/*     matrix, and scale is an output scale factor, set less than or */
/*     equal to 1 to avoid overflow in X. */

/*     In the case of equation (1) the matrix S must be stable (that */
/*     is, all the eigenvalues of S must have negative real parts), */
/*     and for equation (2) the matrix S must be convergent (that is, */
/*     all the eigenvalues of S must lie inside the unit circle). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DISCR   LOGICAL */
/*             Specifies the type of Lyapunov equation to be solved as */
/*             follows: */
/*             = .TRUE. :  Equation (2), discrete-time case; */
/*             = .FALSE.:  Equation (1), continuous-time case. */

/*     LTRANS  LOGICAL */
/*             Specifies the form of op(K) to be used, as follows: */
/*             = .FALSE.:  op(K) = K    (No transpose); */
/*             = .TRUE. :  op(K) = K**T (Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices S and R.  N >= 0. */

/*     S       (input) DOUBLE PRECISION array of dimension (LDS,N) */
/*             The leading N-by-N upper Hessenberg part of this array */
/*             must contain the block upper triangular matrix. */
/*             The elements below the upper Hessenberg part of the array */
/*             S are not referenced. The 2-by-2 blocks must only */
/*             correspond to complex conjugate pairs of eigenvalues (not */
/*             to real eigenvalues). */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,N). */

/*     R       (input/output) DOUBLE PRECISION array of dimension (LDR,N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*             array must contain the upper triangular matrix R. */
/*             On exit, the leading N-by-N upper triangular part of this */
/*             array contains the upper triangular matrix U. */
/*             The strict lower triangle of R is not referenced. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (4*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
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
/*             = 2:  if the matrix S is not stable (that is, one or more */
/*                   of the eigenvalues of S has a non-negative real */
/*                   part), if DISCR = .FALSE., or not convergent (that */
/*                   is, one or more of the eigenvalues of S lies outside */
/*                   the unit circle), if DISCR = .TRUE.; */
/*             = 3:  if the matrix S has two or more consecutive non-zero */
/*                   elements on the first sub-diagonal, so that there is */
/*                   a block larger than 2-by-2 on the diagonal; */
/*             = 4:  if the matrix S has a 2-by-2 diagonal block with */
/*                   real eigenvalues instead of a complex conjugate */
/*                   pair. */

/*     METHOD */

/*     The method used by the routine is based on a variant of the */
/*     Bartels and Stewart backward substitution method [1], that finds */
/*     the Cholesky factor op(U) directly without first finding X and */
/*     without the need to form the normal matrix op(R)'*op(R) [2]. */

/*     The continuous-time Lyapunov equation in the canonical form */
/*                                                        2 */
/*       op(S)'*op(U)'*op(U) + op(U)'*op(U)*op(S) = -scale *op(R)'*op(R), */

/*     or the discrete-time Lyapunov equation in the canonical form */
/*                                                        2 */
/*       op(S)'*op(U)'*op(U)*op(S) - op(U)'*op(U) = -scale *op(R)'*op(R), */

/*     where U and R are upper triangular, is solved for U. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W. */
/*         Solution of the matrix equation  A'X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     [2] Hammarling, S.J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-325, 1982. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     FURTHER COMMENTS */

/*     The Lyapunov equation may be very ill-conditioned. In particular */
/*     if S is only just stable (or convergent) then the Lyapunov */
/*     equation will be ill-conditioned. "Large" elements in U relative */
/*     to those of S and R, or a "small" value for scale, is a symptom */
/*     of ill-conditioning. A condition estimate can be computed using */
/*     SLICOT Library routine SB03MD. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Supersedes Release 2.0 routine SB03CZ by Sven Hammarling, */
/*     NAG Ltd, United Kingdom, Oct. 1986. */
/*     Partly based on SB03CZ and PLYAP1 by A. Varga, University of */
/*     Bochum, May 1992. */

/*     REVISIONS */

/*     Dec. 1997, April 1998, May 1999, Feb. 2004. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 215 "SB03OT.f"
    /* Parameter adjustments */
#line 215 "SB03OT.f"
    s_dim1 = *lds;
#line 215 "SB03OT.f"
    s_offset = 1 + s_dim1;
#line 215 "SB03OT.f"
    s -= s_offset;
#line 215 "SB03OT.f"
    r_dim1 = *ldr;
#line 215 "SB03OT.f"
    r_offset = 1 + r_dim1;
#line 215 "SB03OT.f"
    r__ -= r_offset;
#line 215 "SB03OT.f"
    --dwork;
#line 215 "SB03OT.f"

#line 215 "SB03OT.f"
    /* Function Body */
#line 215 "SB03OT.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 219 "SB03OT.f"
    if (*n < 0) {
#line 220 "SB03OT.f"
	*info = -3;
#line 221 "SB03OT.f"
    } else if (*lds < max(1,*n)) {
#line 222 "SB03OT.f"
	*info = -5;
#line 223 "SB03OT.f"
    } else if (*ldr < max(1,*n)) {
#line 224 "SB03OT.f"
	*info = -7;
#line 225 "SB03OT.f"
    }

#line 227 "SB03OT.f"
    if (*info != 0) {

/*        Error return. */

#line 231 "SB03OT.f"
	i__1 = -(*info);
#line 231 "SB03OT.f"
	xerbla_("SB03OT", &i__1, (ftnlen)6);
#line 232 "SB03OT.f"
	return 0;
#line 233 "SB03OT.f"
    }

#line 235 "SB03OT.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 239 "SB03OT.f"
    if (*n == 0) {
#line 239 "SB03OT.f"
	return 0;
#line 239 "SB03OT.f"
    }

/*     Set constants to control overflow. */

#line 244 "SB03OT.f"
    eps = dlamch_("P", (ftnlen)1);
#line 245 "SB03OT.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 246 "SB03OT.f"
    bignum = 1. / smlnum;
#line 247 "SB03OT.f"
    dlabad_(&smlnum, &bignum);
#line 248 "SB03OT.f"
    smlnum = smlnum * (doublereal) (*n * *n) / eps;
#line 249 "SB03OT.f"
    bignum = 1. / smlnum;

/* Computing MAX */
#line 251 "SB03OT.f"
    d__1 = smlnum, d__2 = eps * dlanhs_("Max", n, &s[s_offset], lds, &dwork[1]
	    , (ftnlen)3);
#line 251 "SB03OT.f"
    smin = max(d__1,d__2);
#line 252 "SB03OT.f"
    infom = 0;

/*     Start the solution. Most of the comments refer to notation and */
/*     equations in sections 5 and 10 of the second reference above. */

/*     Determine whether or not the current block is two-by-two. */
/*     K gives the position of the start of the current block and */
/*     TBYT is true if the block is two-by-two. */

#line 261 "SB03OT.f"
    cont = ! (*discr);
#line 262 "SB03OT.f"
    isgn = 1;
#line 263 "SB03OT.f"
    if (! (*ltrans)) {

/*        Case op(M) = M. */

#line 267 "SB03OT.f"
	kount = 1;

#line 269 "SB03OT.f"
L10:
/*        WHILE( KOUNT.LE.N )LOOP */
#line 271 "SB03OT.f"
	if (kount <= *n) {
#line 272 "SB03OT.f"
	    k = kount;
#line 273 "SB03OT.f"
	    if (kount >= *n) {
#line 274 "SB03OT.f"
		tbyt = FALSE_;
#line 275 "SB03OT.f"
		++kount;
#line 276 "SB03OT.f"
	    } else if (s[k + 1 + k * s_dim1] == 0.) {
#line 277 "SB03OT.f"
		tbyt = FALSE_;
#line 278 "SB03OT.f"
		++kount;
#line 279 "SB03OT.f"
	    } else {
#line 280 "SB03OT.f"
		tbyt = TRUE_;
#line 281 "SB03OT.f"
		if (k + 1 < *n) {
#line 282 "SB03OT.f"
		    if (s[k + 2 + (k + 1) * s_dim1] != 0.) {
#line 283 "SB03OT.f"
			*info = 3;
#line 284 "SB03OT.f"
			return 0;
#line 285 "SB03OT.f"
		    }
#line 286 "SB03OT.f"
		}
#line 287 "SB03OT.f"
		kount += 2;
#line 288 "SB03OT.f"
	    }
#line 289 "SB03OT.f"
	    if (tbyt) {

/*              Solve the two-by-two Lyapunov equation (6.1) or (10.19), */
/*              using the routine SB03OY. */

#line 294 "SB03OT.f"
		b[0] = s[k + k * s_dim1];
#line 295 "SB03OT.f"
		b[1] = s[k + 1 + k * s_dim1];
#line 296 "SB03OT.f"
		b[2] = s[k + (k + 1) * s_dim1];
#line 297 "SB03OT.f"
		b[3] = s[k + 1 + (k + 1) * s_dim1];
#line 298 "SB03OT.f"
		u[0] = r__[k + k * r_dim1];
#line 299 "SB03OT.f"
		u[2] = r__[k + (k + 1) * r_dim1];
#line 300 "SB03OT.f"
		u[3] = r__[k + 1 + (k + 1) * r_dim1];

#line 302 "SB03OT.f"
		sb03oy_(discr, ltrans, &isgn, b, &c__2, u, &c__2, a, &c__2, &
			scaloc, info);
#line 304 "SB03OT.f"
		if (*info > 1) {
#line 304 "SB03OT.f"
		    return 0;
#line 304 "SB03OT.f"
		}
#line 306 "SB03OT.f"
		infom = max(*info,infom);
#line 307 "SB03OT.f"
		if (scaloc != 1.) {

#line 309 "SB03OT.f"
		    i__1 = *n;
#line 309 "SB03OT.f"
		    for (j = 1; j <= i__1; ++j) {
#line 310 "SB03OT.f"
			dscal_(&j, &scaloc, &r__[j * r_dim1 + 1], &c__1);
#line 311 "SB03OT.f"
/* L20: */
#line 311 "SB03OT.f"
		    }

#line 313 "SB03OT.f"
		    *scale *= scaloc;
#line 314 "SB03OT.f"
		}
#line 315 "SB03OT.f"
		r__[k + k * r_dim1] = u[0];
#line 316 "SB03OT.f"
		r__[k + (k + 1) * r_dim1] = u[2];
#line 317 "SB03OT.f"
		r__[k + 1 + (k + 1) * r_dim1] = u[3];

/*              If we are not at the end of S then set up and solve */
/*              equation (6.2) or (10.20). */

/*              Note that  SB03OY  returns  ( u11*s11*inv( u11 ) ) in  B */
/*              and returns scaled alpha in  A.  ksize is the order of */
/*              the remainder of  S.  k1, k2 and k3  point to the start */
/*              of vectors in  DWORK. */

#line 327 "SB03OT.f"
		if (kount <= *n) {
#line 328 "SB03OT.f"
		    ksize = *n - k - 1;
#line 329 "SB03OT.f"
		    k1 = ksize + 1;
#line 330 "SB03OT.f"
		    k2 = ksize + k1;
#line 331 "SB03OT.f"
		    k3 = ksize + k2;

/*                 Form the right-hand side of (6.2) or (10.20), the */
/*                 first column in DWORK( 1 ) ,..., DWORK( n - k - 1 ) */
/*                 the second in DWORK( n - k ) ,..., */
/*                 DWORK( 2*( n - k - 1 ) ). */

#line 338 "SB03OT.f"
		    dcopy_(&ksize, &r__[k + (k + 2) * r_dim1], ldr, &dwork[1],
			     &c__1);
#line 339 "SB03OT.f"
		    dcopy_(&ksize, &r__[k + 1 + (k + 2) * r_dim1], ldr, &
			    dwork[k1], &c__1);
#line 340 "SB03OT.f"
		    dtrmm_("Right", "Upper", "No transpose", "Non-unit", &
			    ksize, &c__2, &c_b19, a, &c__2, &dwork[1], &ksize,
			     (ftnlen)5, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 343 "SB03OT.f"
		    if (cont) {
#line 344 "SB03OT.f"
			d__1 = -r__[k + k * r_dim1];
#line 344 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + (k + 2) * s_dim1], lds, &
				dwork[1], &c__1);
#line 346 "SB03OT.f"
			d__1 = -r__[k + (k + 1) * r_dim1];
#line 346 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + 1 + (k + 2) * s_dim1], 
				lds, &dwork[1], &c__1);
#line 348 "SB03OT.f"
			d__1 = -r__[k + 1 + (k + 1) * r_dim1];
#line 348 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + 1 + (k + 2) * s_dim1], 
				lds, &dwork[k1], &c__1);
#line 350 "SB03OT.f"
		    } else {
#line 351 "SB03OT.f"
			d__1 = -r__[k + k * r_dim1] * b[0];
#line 351 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + (k + 2) * s_dim1], lds, &
				dwork[1], &c__1);
#line 353 "SB03OT.f"
			d__1 = -(r__[k + (k + 1) * r_dim1] * b[0] + r__[k + 1 
				+ (k + 1) * r_dim1] * b[1]);
#line 353 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + 1 + (k + 2) * s_dim1], 
				lds, &dwork[1], &c__1);
#line 355 "SB03OT.f"
			d__1 = -r__[k + k * r_dim1] * b[2];
#line 355 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + (k + 2) * s_dim1], lds, &
				dwork[k1], &c__1);
#line 357 "SB03OT.f"
			d__1 = -(r__[k + (k + 1) * r_dim1] * b[2] + r__[k + 1 
				+ (k + 1) * r_dim1] * b[3]);
#line 357 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + 1 + (k + 2) * s_dim1], 
				lds, &dwork[k1], &c__1);
#line 360 "SB03OT.f"
		    }

/*                 SB03OR  solves the Sylvester equations. The solution */
/*                 is overwritten on DWORK. */

#line 365 "SB03OT.f"
		    sb03or_(discr, ltrans, &ksize, &c__2, &s[k + 2 + (k + 2) *
			     s_dim1], lds, b, &c__2, &dwork[1], &ksize, &
			    scaloc, info);
#line 367 "SB03OT.f"
		    infom = max(*info,infom);
#line 368 "SB03OT.f"
		    if (scaloc != 1.) {

#line 370 "SB03OT.f"
			i__1 = *n;
#line 370 "SB03OT.f"
			for (j = 1; j <= i__1; ++j) {
#line 371 "SB03OT.f"
			    dscal_(&j, &scaloc, &r__[j * r_dim1 + 1], &c__1);
#line 372 "SB03OT.f"
/* L30: */
#line 372 "SB03OT.f"
			}

#line 374 "SB03OT.f"
			*scale *= scaloc;
#line 375 "SB03OT.f"
		    }

/*                 Copy the solution into the next  2*( n - k - 1 ) */
/*                 elements of  DWORK. */

#line 380 "SB03OT.f"
		    i__1 = ksize << 1;
#line 380 "SB03OT.f"
		    dcopy_(&i__1, &dwork[1], &c__1, &dwork[k2], &c__1);

/*                 Now form the matrix  Rhat  of equation (6.4) or */
/*                 (10.22). Note that (10.22) is incorrect, so here we */
/*                 implement a corrected version of (10.22). */

#line 386 "SB03OT.f"
		    if (cont) {

/*                    Swap the two rows of R with DWORK. */

#line 390 "SB03OT.f"
			dswap_(&ksize, &dwork[1], &c__1, &r__[k + (k + 2) * 
				r_dim1], ldr);
#line 391 "SB03OT.f"
			dswap_(&ksize, &dwork[k1], &c__1, &r__[k + 1 + (k + 2)
				 * r_dim1], ldr);

/*                    1st column: */

#line 395 "SB03OT.f"
			d__1 = -a[0];
#line 395 "SB03OT.f"
			daxpy_(&ksize, &d__1, &dwork[k2], &c__1, &dwork[1], &
				c__1);
#line 397 "SB03OT.f"
			d__1 = -a[2];
#line 397 "SB03OT.f"
			daxpy_(&ksize, &d__1, &dwork[k3], &c__1, &dwork[1], &
				c__1);

/*                    2nd column: */

#line 402 "SB03OT.f"
			d__1 = -a[3];
#line 402 "SB03OT.f"
			daxpy_(&ksize, &d__1, &dwork[k3], &c__1, &dwork[k1], &
				c__1);
#line 404 "SB03OT.f"
		    } else {

/*                    Form  v = S1'*u + s*u11', overwriting  v  on DWORK. */

/*                    Compute  S1'*u,  first multiplying by the */
/*                    triangular part of  S1. */

#line 411 "SB03OT.f"
			dtrmm_("Left", "Upper", "Transpose", "Non-unit", &
				ksize, &c__2, &c_b47, &s[k + 2 + (k + 2) * 
				s_dim1], lds, &dwork[1], &ksize, (ftnlen)4, (
				ftnlen)5, (ftnlen)9, (ftnlen)8);

/*                    Then multiply by the subdiagonal of  S1  and add in */
/*                    to the above result. */

#line 418 "SB03OT.f"
			j1 = k1;
#line 419 "SB03OT.f"
			j2 = k + 2;

#line 421 "SB03OT.f"
			i__1 = ksize - 1;
#line 421 "SB03OT.f"
			for (j = 1; j <= i__1; ++j) {
#line 422 "SB03OT.f"
			    if (s[j2 + 1 + j2 * s_dim1] != 0.) {
#line 423 "SB03OT.f"
				dwork[j] = s[j2 + 1 + j2 * s_dim1] * dwork[k2 
					+ j] + dwork[j];
#line 424 "SB03OT.f"
				dwork[j1] = s[j2 + 1 + j2 * s_dim1] * dwork[
					k3 + j] + dwork[j1];
#line 426 "SB03OT.f"
			    }
#line 427 "SB03OT.f"
			    ++j1;
#line 428 "SB03OT.f"
			    ++j2;
#line 429 "SB03OT.f"
/* L40: */
#line 429 "SB03OT.f"
			}

/*                    Add in s*u11'. */

#line 433 "SB03OT.f"
			daxpy_(&ksize, &r__[k + k * r_dim1], &s[k + (k + 2) * 
				s_dim1], lds, &dwork[1], &c__1);
#line 435 "SB03OT.f"
			daxpy_(&ksize, &r__[k + (k + 1) * r_dim1], &s[k + 1 + 
				(k + 2) * s_dim1], lds, &dwork[1], &c__1);
#line 437 "SB03OT.f"
			daxpy_(&ksize, &r__[k + 1 + (k + 1) * r_dim1], &s[k + 
				1 + (k + 2) * s_dim1], lds, &dwork[k1], &c__1)
				;

/*                    Next recover r from R, swapping r with u. */

#line 442 "SB03OT.f"
			dswap_(&ksize, &dwork[k2], &c__1, &r__[k + (k + 2) * 
				r_dim1], ldr);
#line 443 "SB03OT.f"
			dswap_(&ksize, &dwork[k3], &c__1, &r__[k + 1 + (k + 2)
				 * r_dim1], ldr);

/*                    Now we perform the QR factorization. */

/*                    ( a ) = Q*( t ), */
/*                    ( b ) */

/*                    and form */

/*                    ( p' ) = Q'*( r' ). */
/*                    ( y' )      ( v' ) */

/*                    y  is then the correct vector to use in (10.22). */
/*                    Note that  a  is upper triangular and that  t  and */
/*                    p  are not required. */

#line 459 "SB03OT.f"
			dlarfg_(&c__3, a, b, &c__1, &tau1);
#line 460 "SB03OT.f"
			v1 = b[0];
#line 461 "SB03OT.f"
			t1 = tau1 * v1;
#line 462 "SB03OT.f"
			v2 = b[1];
#line 463 "SB03OT.f"
			t2 = tau1 * v2;
#line 464 "SB03OT.f"
			sum = a[2] + v1 * b[2] + v2 * b[3];
#line 465 "SB03OT.f"
			b[2] -= sum * t1;
#line 466 "SB03OT.f"
			b[3] -= sum * t2;
#line 467 "SB03OT.f"
			dlarfg_(&c__3, &a[3], &b[2], &c__1, &tau2);
#line 468 "SB03OT.f"
			v3 = b[2];
#line 469 "SB03OT.f"
			t3 = tau2 * v3;
#line 470 "SB03OT.f"
			v4 = b[3];
#line 471 "SB03OT.f"
			t4 = tau2 * v4;
#line 472 "SB03OT.f"
			j1 = k1;
#line 473 "SB03OT.f"
			j2 = k2;
#line 474 "SB03OT.f"
			j3 = k3;

#line 476 "SB03OT.f"
			i__1 = ksize;
#line 476 "SB03OT.f"
			for (j = 1; j <= i__1; ++j) {
#line 477 "SB03OT.f"
			    sum = dwork[j2] + v1 * dwork[j] + v2 * dwork[j1];
#line 478 "SB03OT.f"
			    d1 = dwork[j] - sum * t1;
#line 479 "SB03OT.f"
			    d2 = dwork[j1] - sum * t2;
#line 480 "SB03OT.f"
			    sum = dwork[j3] + v3 * d1 + v4 * d2;
#line 481 "SB03OT.f"
			    dwork[j] = d1 - sum * t3;
#line 482 "SB03OT.f"
			    dwork[j1] = d2 - sum * t4;
#line 483 "SB03OT.f"
			    ++j1;
#line 484 "SB03OT.f"
			    ++j2;
#line 485 "SB03OT.f"
			    ++j3;
#line 486 "SB03OT.f"
/* L50: */
#line 486 "SB03OT.f"
			}

#line 488 "SB03OT.f"
		    }

/*                 Now update  R1  to give  Rhat. */

#line 492 "SB03OT.f"
		    dcopy_(&ksize, &dwork[1], &c__1, &dwork[k2], &c__1);
#line 493 "SB03OT.f"
		    dcopy_(&ksize, &dwork[k1], &c__1, &dwork[k3], &c__1);
#line 494 "SB03OT.f"
		    dcopy_(&ksize, &dwork[k3], &c__1, &dwork[2], &c__2);
#line 495 "SB03OT.f"
		    dcopy_(&ksize, &dwork[k2], &c__1, &dwork[1], &c__2);
#line 496 "SB03OT.f"
		    mb04od_("Full", &ksize, &c__0, &c__2, &r__[k + 2 + (k + 2)
			     * r_dim1], ldr, &dwork[1], &c__2, &dwork[1], &
			    c__1, &dwork[1], &c__1, &dwork[k2], &dwork[k3], (
			    ftnlen)4);
#line 499 "SB03OT.f"
		}
#line 500 "SB03OT.f"
	    } else {

/*              1-by-1 block. */

/*              Make sure S is stable or convergent and find u11 in */
/*              equation (5.13) or (10.15). */

#line 507 "SB03OT.f"
		if (*discr) {
#line 508 "SB03OT.f"
		    absskk = (d__1 = s[k + k * s_dim1], abs(d__1));
#line 509 "SB03OT.f"
		    if (absskk - 1. >= 0.) {
#line 510 "SB03OT.f"
			*info = 2;
#line 511 "SB03OT.f"
			return 0;
#line 512 "SB03OT.f"
		    }
#line 513 "SB03OT.f"
		    temp = sqrt((1. - absskk) * (absskk + 1.));
#line 514 "SB03OT.f"
		} else {
#line 515 "SB03OT.f"
		    if (s[k + k * s_dim1] >= 0.) {
#line 516 "SB03OT.f"
			*info = 2;
#line 517 "SB03OT.f"
			return 0;
#line 518 "SB03OT.f"
		    }
#line 519 "SB03OT.f"
		    temp = sqrt((d__1 = s[k + k * s_dim1] * 2., abs(d__1)));
#line 520 "SB03OT.f"
		}

#line 522 "SB03OT.f"
		scaloc = 1.;
#line 523 "SB03OT.f"
		if (temp < smin) {
#line 524 "SB03OT.f"
		    temp = smin;
#line 525 "SB03OT.f"
		    infom = 1;
#line 526 "SB03OT.f"
		}
#line 527 "SB03OT.f"
		dr = (d__1 = r__[k + k * r_dim1], abs(d__1));
#line 528 "SB03OT.f"
		if (temp < 1. && dr > 1.) {
#line 529 "SB03OT.f"
		    if (dr > bignum * temp) {
#line 529 "SB03OT.f"
			scaloc = 1. / dr;
#line 529 "SB03OT.f"
		    }
#line 531 "SB03OT.f"
		}
#line 532 "SB03OT.f"
		alpha = d_sign(&temp, &r__[k + k * r_dim1]);
#line 533 "SB03OT.f"
		r__[k + k * r_dim1] /= alpha;
#line 534 "SB03OT.f"
		if (scaloc != 1.) {

#line 536 "SB03OT.f"
		    i__1 = *n;
#line 536 "SB03OT.f"
		    for (j = 1; j <= i__1; ++j) {
#line 537 "SB03OT.f"
			dscal_(&j, &scaloc, &r__[j * r_dim1 + 1], &c__1);
#line 538 "SB03OT.f"
/* L60: */
#line 538 "SB03OT.f"
		    }

#line 540 "SB03OT.f"
		    *scale *= scaloc;
#line 541 "SB03OT.f"
		}

/*              If we are not at the end of  S  then set up and solve */
/*              equation (5.14) or (10.16).  ksize is the order of the */
/*              remainder of  S.  k1 and k2 point to the start of vectors */
/*              in  DWORK. */

#line 548 "SB03OT.f"
		if (kount <= *n) {
#line 549 "SB03OT.f"
		    ksize = *n - k;
#line 550 "SB03OT.f"
		    k1 = ksize + 1;
#line 551 "SB03OT.f"
		    k2 = ksize + k1;

/*                 Form the right-hand side in DWORK( 1 ),..., */
/*                 DWORK( n - k ). */

#line 556 "SB03OT.f"
		    dcopy_(&ksize, &r__[k + (k + 1) * r_dim1], ldr, &dwork[1],
			     &c__1);
#line 557 "SB03OT.f"
		    d__1 = -alpha;
#line 557 "SB03OT.f"
		    dscal_(&ksize, &d__1, &dwork[1], &c__1);
#line 558 "SB03OT.f"
		    if (cont) {
#line 559 "SB03OT.f"
			d__1 = -r__[k + k * r_dim1];
#line 559 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + (k + 1) * s_dim1], lds, &
				dwork[1], &c__1);
#line 561 "SB03OT.f"
		    } else {
#line 562 "SB03OT.f"
			d__1 = -s[k + k * s_dim1] * r__[k + k * r_dim1];
#line 562 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + (k + 1) * s_dim1], lds, &
				dwork[1], &c__1);
#line 564 "SB03OT.f"
		    }

/*                 SB03OR solves the Sylvester equations. The solution is */
/*                 overwritten on  DWORK. */

#line 569 "SB03OT.f"
		    sb03or_(discr, ltrans, &ksize, &c__1, &s[k + 1 + (k + 1) *
			     s_dim1], lds, &s[k + k * s_dim1], &c__1, &dwork[
			    1], &ksize, &scaloc, info);
#line 571 "SB03OT.f"
		    infom = max(*info,infom);
#line 572 "SB03OT.f"
		    if (scaloc != 1.) {

#line 574 "SB03OT.f"
			i__1 = *n;
#line 574 "SB03OT.f"
			for (j = 1; j <= i__1; ++j) {
#line 575 "SB03OT.f"
			    dscal_(&j, &scaloc, &r__[j * r_dim1 + 1], &c__1);
#line 576 "SB03OT.f"
/* L70: */
#line 576 "SB03OT.f"
			}

#line 578 "SB03OT.f"
			*scale *= scaloc;
#line 579 "SB03OT.f"
		    }

/*                 Copy the solution into the next  ( n - k ) elements */
/*                 of  DWORK,  copy the solution back into  R  and copy */
/*                 the row of  R  back into  DWORK. */

#line 585 "SB03OT.f"
		    dcopy_(&ksize, &dwork[1], &c__1, &dwork[k1], &c__1);
#line 586 "SB03OT.f"
		    dswap_(&ksize, &dwork[1], &c__1, &r__[k + (k + 1) * 
			    r_dim1], ldr);

/*                 Now form the matrix  Rhat  of equation (5.15) or */
/*                 (10.17), first computing  y  in  DWORK,  and then */
/*                 updating  R1. */

#line 592 "SB03OT.f"
		    if (cont) {
#line 593 "SB03OT.f"
			d__1 = -alpha;
#line 593 "SB03OT.f"
			daxpy_(&ksize, &d__1, &dwork[k1], &c__1, &dwork[1], &
				c__1);
#line 594 "SB03OT.f"
		    } else {

/*                    First form  lambda( 1 )*r  and then add in */
/*                    alpha*u11*s. */

#line 599 "SB03OT.f"
			d__1 = -s[k + k * s_dim1];
#line 599 "SB03OT.f"
			dscal_(&ksize, &d__1, &dwork[1], &c__1);
#line 600 "SB03OT.f"
			d__1 = alpha * r__[k + k * r_dim1];
#line 600 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k + (k + 1) * s_dim1], lds, &
				dwork[1], &c__1);

/*                    Now form  alpha*S1'*u,  first multiplying by the */
/*                    sub-diagonal of  S1  and then the triangular part */
/*                    of  S1,  and add the result in DWORK. */

#line 607 "SB03OT.f"
			j1 = k + 1;

#line 609 "SB03OT.f"
			i__1 = ksize - 1;
#line 609 "SB03OT.f"
			for (j = 1; j <= i__1; ++j) {
#line 610 "SB03OT.f"
			    if (s[j1 + 1 + j1 * s_dim1] != 0.) {
#line 610 "SB03OT.f"
				dwork[j] = alpha * s[j1 + 1 + j1 * s_dim1] * 
					dwork[k1 + j] + dwork[j];
#line 610 "SB03OT.f"
			    }
#line 612 "SB03OT.f"
			    ++j1;
#line 613 "SB03OT.f"
/* L80: */
#line 613 "SB03OT.f"
			}

#line 615 "SB03OT.f"
			dtrmv_("Upper", "Transpose", "Non-unit", &ksize, &s[k 
				+ 1 + (k + 1) * s_dim1], lds, &dwork[k1], &
				c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 617 "SB03OT.f"
			daxpy_(&ksize, &alpha, &dwork[k1], &c__1, &dwork[1], &
				c__1);
#line 618 "SB03OT.f"
		    }
#line 619 "SB03OT.f"
		    mb04od_("Full", &ksize, &c__0, &c__1, &r__[k + 1 + (k + 1)
			     * r_dim1], ldr, &dwork[1], &c__1, &dwork[1], &
			    c__1, &dwork[1], &c__1, &dwork[k2], &dwork[k1], (
			    ftnlen)4);
#line 622 "SB03OT.f"
		}
#line 623 "SB03OT.f"
	    }
#line 624 "SB03OT.f"
	    goto L10;
#line 625 "SB03OT.f"
	}
/*        END WHILE 10 */

#line 628 "SB03OT.f"
    } else {

/*        Case op(M) = M'. */

#line 632 "SB03OT.f"
	kount = *n;

#line 634 "SB03OT.f"
L90:
/*        WHILE( KOUNT.GE.1 )LOOP */
#line 636 "SB03OT.f"
	if (kount >= 1) {
#line 637 "SB03OT.f"
	    k = kount;
#line 638 "SB03OT.f"
	    if (kount == 1) {
#line 639 "SB03OT.f"
		tbyt = FALSE_;
#line 640 "SB03OT.f"
		--kount;
#line 641 "SB03OT.f"
	    } else if (s[k + (k - 1) * s_dim1] == 0.) {
#line 642 "SB03OT.f"
		tbyt = FALSE_;
#line 643 "SB03OT.f"
		--kount;
#line 644 "SB03OT.f"
	    } else {
#line 645 "SB03OT.f"
		tbyt = TRUE_;
#line 646 "SB03OT.f"
		--k;
#line 647 "SB03OT.f"
		if (k > 1) {
#line 648 "SB03OT.f"
		    if (s[k + (k - 1) * s_dim1] != 0.) {
#line 649 "SB03OT.f"
			*info = 3;
#line 650 "SB03OT.f"
			return 0;
#line 651 "SB03OT.f"
		    }
#line 652 "SB03OT.f"
		}
#line 653 "SB03OT.f"
		kount += -2;
#line 654 "SB03OT.f"
	    }
#line 655 "SB03OT.f"
	    if (tbyt) {

/*              Solve the two-by-two Lyapunov equation corresponding to */
/*              (6.1) or (10.19), using the routine SB03OY. */

#line 660 "SB03OT.f"
		b[0] = s[k + k * s_dim1];
#line 661 "SB03OT.f"
		b[1] = s[k + 1 + k * s_dim1];
#line 662 "SB03OT.f"
		b[2] = s[k + (k + 1) * s_dim1];
#line 663 "SB03OT.f"
		b[3] = s[k + 1 + (k + 1) * s_dim1];
#line 664 "SB03OT.f"
		u[0] = r__[k + k * r_dim1];
#line 665 "SB03OT.f"
		u[2] = r__[k + (k + 1) * r_dim1];
#line 666 "SB03OT.f"
		u[3] = r__[k + 1 + (k + 1) * r_dim1];

#line 668 "SB03OT.f"
		sb03oy_(discr, ltrans, &isgn, b, &c__2, u, &c__2, a, &c__2, &
			scaloc, info);
#line 670 "SB03OT.f"
		if (*info > 1) {
#line 670 "SB03OT.f"
		    return 0;
#line 670 "SB03OT.f"
		}
#line 672 "SB03OT.f"
		infom = max(*info,infom);
#line 673 "SB03OT.f"
		if (scaloc != 1.) {

#line 675 "SB03OT.f"
		    i__1 = *n;
#line 675 "SB03OT.f"
		    for (j = 1; j <= i__1; ++j) {
#line 676 "SB03OT.f"
			dscal_(&j, &scaloc, &r__[j * r_dim1 + 1], &c__1);
#line 677 "SB03OT.f"
/* L100: */
#line 677 "SB03OT.f"
		    }

#line 679 "SB03OT.f"
		    *scale *= scaloc;
#line 680 "SB03OT.f"
		}
#line 681 "SB03OT.f"
		r__[k + k * r_dim1] = u[0];
#line 682 "SB03OT.f"
		r__[k + (k + 1) * r_dim1] = u[2];
#line 683 "SB03OT.f"
		r__[k + 1 + (k + 1) * r_dim1] = u[3];

/*              If we are not at the front of S then set up and solve */
/*              equation corresponding to (6.2) or (10.20). */

/*              Note that  SB03OY  returns  ( inv( u11 )*s11*u11 ) in  B */
/*              and returns scaled alpha, alpha = inv( u11 )*r11, in  A. */
/*              ksize is the order of the remainder leading part of  S. */
/*              k1, k2 and k3 point to the start of vectors in  DWORK. */

#line 693 "SB03OT.f"
		if (kount >= 1) {
#line 694 "SB03OT.f"
		    ksize = k - 1;
#line 695 "SB03OT.f"
		    k1 = ksize + 1;
#line 696 "SB03OT.f"
		    k2 = ksize + k1;
#line 697 "SB03OT.f"
		    k3 = ksize + k2;

/*                 Form the right-hand side of equations corresponding to */
/*                 (6.2) or (10.20), the first column in DWORK( 1 ) ,..., */
/*                 DWORK( k - 1 ) the second in DWORK( k ) ,..., */
/*                 DWORK( 2*( k - 1 ) ). */

#line 704 "SB03OT.f"
		    dcopy_(&ksize, &r__[k * r_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 705 "SB03OT.f"
		    dcopy_(&ksize, &r__[(k + 1) * r_dim1 + 1], &c__1, &dwork[
			    k1], &c__1);
#line 706 "SB03OT.f"
		    dtrmm_("Right", "Upper", "Transpose", "Non-unit", &ksize, 
			    &c__2, &c_b19, a, &c__2, &dwork[1], &ksize, (
			    ftnlen)5, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 708 "SB03OT.f"
		    if (cont) {
#line 709 "SB03OT.f"
			d__1 = -r__[k + k * r_dim1];
#line 709 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k * s_dim1 + 1], &c__1, &
				dwork[1], &c__1);
#line 710 "SB03OT.f"
			d__1 = -r__[k + (k + 1) * r_dim1];
#line 710 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k * s_dim1 + 1], &c__1, &
				dwork[k1], &c__1);
#line 712 "SB03OT.f"
			d__1 = -r__[k + 1 + (k + 1) * r_dim1];
#line 712 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[(k + 1) * s_dim1 + 1], &c__1,
				 &dwork[k1], &c__1);
#line 714 "SB03OT.f"
		    } else {
#line 715 "SB03OT.f"
			d__1 = -(r__[k + k * r_dim1] * b[0] + r__[k + (k + 1) 
				* r_dim1] * b[2]);
#line 715 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k * s_dim1 + 1], &c__1, &
				dwork[1], &c__1);
#line 717 "SB03OT.f"
			d__1 = -r__[k + 1 + (k + 1) * r_dim1] * b[2];
#line 717 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[(k + 1) * s_dim1 + 1], &c__1,
				 &dwork[1], &c__1);
#line 719 "SB03OT.f"
			d__1 = -(r__[k + k * r_dim1] * b[1] + r__[k + (k + 1) 
				* r_dim1] * b[3]);
#line 719 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k * s_dim1 + 1], &c__1, &
				dwork[k1], &c__1);
#line 721 "SB03OT.f"
			d__1 = -r__[k + 1 + (k + 1) * r_dim1] * b[3];
#line 721 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[(k + 1) * s_dim1 + 1], &c__1,
				 &dwork[k1], &c__1);
#line 723 "SB03OT.f"
		    }

/*                 SB03OR  solves the Sylvester equations. The solution */
/*                 is overwritten on DWORK. */

#line 728 "SB03OT.f"
		    sb03or_(discr, ltrans, &ksize, &c__2, &s[s_offset], lds, 
			    b, &c__2, &dwork[1], &ksize, &scaloc, info);
#line 730 "SB03OT.f"
		    infom = max(*info,infom);
#line 731 "SB03OT.f"
		    if (scaloc != 1.) {

#line 733 "SB03OT.f"
			i__1 = *n;
#line 733 "SB03OT.f"
			for (j = 1; j <= i__1; ++j) {
#line 734 "SB03OT.f"
			    dscal_(&j, &scaloc, &r__[j * r_dim1 + 1], &c__1);
#line 735 "SB03OT.f"
/* L110: */
#line 735 "SB03OT.f"
			}

#line 737 "SB03OT.f"
			*scale *= scaloc;
#line 738 "SB03OT.f"
		    }

/*                 Copy the solution into the next  2*( k - 1 ) elements */
/*                 of  DWORK. */

#line 743 "SB03OT.f"
		    i__1 = ksize << 1;
#line 743 "SB03OT.f"
		    dcopy_(&i__1, &dwork[1], &c__1, &dwork[k2], &c__1);

/*                 Now form the matrix  Rhat  of equation corresponding */
/*                 to (6.4) or (10.22) (corrected version). */

#line 748 "SB03OT.f"
		    if (cont) {

/*                    Swap the two columns of R with DWORK. */

#line 752 "SB03OT.f"
			dswap_(&ksize, &dwork[1], &c__1, &r__[k * r_dim1 + 1],
				 &c__1);
#line 753 "SB03OT.f"
			dswap_(&ksize, &dwork[k1], &c__1, &r__[(k + 1) * 
				r_dim1 + 1], &c__1);

/*                    1st column: */

#line 757 "SB03OT.f"
			d__1 = -a[0];
#line 757 "SB03OT.f"
			daxpy_(&ksize, &d__1, &dwork[k2], &c__1, &dwork[1], &
				c__1);

/*                    2nd column: */

#line 762 "SB03OT.f"
			d__1 = -a[2];
#line 762 "SB03OT.f"
			daxpy_(&ksize, &d__1, &dwork[k2], &c__1, &dwork[k1], &
				c__1);
#line 764 "SB03OT.f"
			d__1 = -a[3];
#line 764 "SB03OT.f"
			daxpy_(&ksize, &d__1, &dwork[k3], &c__1, &dwork[k1], &
				c__1);
#line 766 "SB03OT.f"
		    } else {

/*                    Form  v = S1*u + s*u11, overwriting  v  on DWORK. */

/*                    Compute  S1*u,  first multiplying by the triangular */
/*                    part of  S1. */

#line 773 "SB03OT.f"
			dtrmm_("Left", "Upper", "No transpose", "Non-unit", &
				ksize, &c__2, &c_b47, &s[s_offset], lds, &
				dwork[1], &ksize, (ftnlen)4, (ftnlen)5, (
				ftnlen)12, (ftnlen)8);

/*                    Then multiply by the subdiagonal of  S1  and add in */
/*                    to the above result. */

#line 780 "SB03OT.f"
			j1 = k1;

#line 782 "SB03OT.f"
			i__1 = ksize;
#line 782 "SB03OT.f"
			for (j = 2; j <= i__1; ++j) {
#line 783 "SB03OT.f"
			    ++j1;
#line 784 "SB03OT.f"
			    if (s[j + (j - 1) * s_dim1] != 0.) {
#line 785 "SB03OT.f"
				dwork[j] = s[j + (j - 1) * s_dim1] * dwork[k2 
					+ j - 2] + dwork[j];
#line 786 "SB03OT.f"
				dwork[j1] = s[j + (j - 1) * s_dim1] * dwork[
					k3 + j - 2] + dwork[j1];
#line 788 "SB03OT.f"
			    }
#line 789 "SB03OT.f"
/* L120: */
#line 789 "SB03OT.f"
			}

/*                    Add in s*u11. */

#line 793 "SB03OT.f"
			daxpy_(&ksize, &r__[k + k * r_dim1], &s[k * s_dim1 + 
				1], &c__1, &dwork[1], &c__1);
#line 794 "SB03OT.f"
			daxpy_(&ksize, &r__[k + (k + 1) * r_dim1], &s[k * 
				s_dim1 + 1], &c__1, &dwork[k1], &c__1);
#line 796 "SB03OT.f"
			daxpy_(&ksize, &r__[k + 1 + (k + 1) * r_dim1], &s[(k 
				+ 1) * s_dim1 + 1], &c__1, &dwork[k1], &c__1);

/*                    Next recover r from R, swapping r with u. */

#line 801 "SB03OT.f"
			dswap_(&ksize, &dwork[k2], &c__1, &r__[k * r_dim1 + 1]
				, &c__1);
#line 802 "SB03OT.f"
			dswap_(&ksize, &dwork[k3], &c__1, &r__[(k + 1) * 
				r_dim1 + 1], &c__1);

/*                    Now we perform the QL factorization. */

/*                    ( a' ) = Q*( t ), */
/*                    ( b' ) */

/*                    and form */

/*                    ( p' ) = Q'*( r' ). */
/*                    ( y' )      ( v' ) */

/*                    y  is then the correct vector to use in the */
/*                    relation corresponding to (10.22). */
/*                    Note that  a  is upper triangular and that  t  and */
/*                    p  are not required. */

#line 819 "SB03OT.f"
			dlarfg_(&c__3, &a[3], &b[1], &c__2, &tau1);
#line 820 "SB03OT.f"
			v1 = b[1];
#line 821 "SB03OT.f"
			t1 = tau1 * v1;
#line 822 "SB03OT.f"
			v2 = b[3];
#line 823 "SB03OT.f"
			t2 = tau1 * v2;
#line 824 "SB03OT.f"
			sum = a[2] + v1 * b[0] + v2 * b[2];
#line 825 "SB03OT.f"
			b[0] -= sum * t1;
#line 826 "SB03OT.f"
			b[2] -= sum * t2;
#line 827 "SB03OT.f"
			dlarfg_(&c__3, a, b, &c__2, &tau2);
#line 828 "SB03OT.f"
			v3 = b[0];
#line 829 "SB03OT.f"
			t3 = tau2 * v3;
#line 830 "SB03OT.f"
			v4 = b[2];
#line 831 "SB03OT.f"
			t4 = tau2 * v4;
#line 832 "SB03OT.f"
			j1 = k1;
#line 833 "SB03OT.f"
			j2 = k2;
#line 834 "SB03OT.f"
			j3 = k3;

#line 836 "SB03OT.f"
			i__1 = ksize;
#line 836 "SB03OT.f"
			for (j = 1; j <= i__1; ++j) {
#line 837 "SB03OT.f"
			    sum = dwork[j3] + v1 * dwork[j] + v2 * dwork[j1];
#line 838 "SB03OT.f"
			    d1 = dwork[j] - sum * t1;
#line 839 "SB03OT.f"
			    d2 = dwork[j1] - sum * t2;
#line 840 "SB03OT.f"
			    sum = dwork[j2] + v3 * d1 + v4 * d2;
#line 841 "SB03OT.f"
			    dwork[j] = d1 - sum * t3;
#line 842 "SB03OT.f"
			    dwork[j1] = d2 - sum * t4;
#line 843 "SB03OT.f"
			    ++j1;
#line 844 "SB03OT.f"
			    ++j2;
#line 845 "SB03OT.f"
			    ++j3;
#line 846 "SB03OT.f"
/* L130: */
#line 846 "SB03OT.f"
			}

#line 848 "SB03OT.f"
		    }

/*                 Now update  R1  to give  Rhat. */

#line 852 "SB03OT.f"
		    mb04nd_("Full", &ksize, &c__0, &c__2, &r__[r_offset], ldr,
			     &dwork[1], &ksize, &dwork[1], &c__1, &dwork[1], &
			    c__1, &dwork[k2], &dwork[k3], (ftnlen)4);
#line 855 "SB03OT.f"
		}
#line 856 "SB03OT.f"
	    } else {

/*              1-by-1 block. */

/*              Make sure S is stable or convergent and find u11 in */
/*              equation corresponding to (5.13) or (10.15). */

#line 863 "SB03OT.f"
		if (*discr) {
#line 864 "SB03OT.f"
		    absskk = (d__1 = s[k + k * s_dim1], abs(d__1));
#line 865 "SB03OT.f"
		    if (absskk - 1. >= 0.) {
#line 866 "SB03OT.f"
			*info = 2;
#line 867 "SB03OT.f"
			return 0;
#line 868 "SB03OT.f"
		    }
#line 869 "SB03OT.f"
		    temp = sqrt((1. - absskk) * (absskk + 1.));
#line 870 "SB03OT.f"
		} else {
#line 871 "SB03OT.f"
		    if (s[k + k * s_dim1] >= 0.) {
#line 872 "SB03OT.f"
			*info = 2;
#line 873 "SB03OT.f"
			return 0;
#line 874 "SB03OT.f"
		    }
#line 875 "SB03OT.f"
		    temp = sqrt((d__1 = s[k + k * s_dim1] * 2., abs(d__1)));
#line 876 "SB03OT.f"
		}

#line 878 "SB03OT.f"
		scaloc = 1.;
#line 879 "SB03OT.f"
		if (temp < smin) {
#line 880 "SB03OT.f"
		    temp = smin;
#line 881 "SB03OT.f"
		    infom = 1;
#line 882 "SB03OT.f"
		}
#line 883 "SB03OT.f"
		dr = (d__1 = r__[k + k * r_dim1], abs(d__1));
#line 884 "SB03OT.f"
		if (temp < 1. && dr > 1.) {
#line 885 "SB03OT.f"
		    if (dr > bignum * temp) {
#line 885 "SB03OT.f"
			scaloc = 1. / dr;
#line 885 "SB03OT.f"
		    }
#line 887 "SB03OT.f"
		}
#line 888 "SB03OT.f"
		alpha = d_sign(&temp, &r__[k + k * r_dim1]);
#line 889 "SB03OT.f"
		r__[k + k * r_dim1] /= alpha;
#line 890 "SB03OT.f"
		if (scaloc != 1.) {

#line 892 "SB03OT.f"
		    i__1 = *n;
#line 892 "SB03OT.f"
		    for (j = 1; j <= i__1; ++j) {
#line 893 "SB03OT.f"
			dscal_(&j, &scaloc, &r__[j * r_dim1 + 1], &c__1);
#line 894 "SB03OT.f"
/* L140: */
#line 894 "SB03OT.f"
		    }

#line 896 "SB03OT.f"
		    *scale *= scaloc;
#line 897 "SB03OT.f"
		}

/*              If we are not at the front of  S  then set up and solve */
/*              equation corresponding to (5.14) or (10.16).  ksize is */
/*              the order of the remainder leading part of  S.  k1 and k2 */
/*              point to the start of vectors in  DWORK. */

#line 904 "SB03OT.f"
		if (kount >= 1) {
#line 905 "SB03OT.f"
		    ksize = k - 1;
#line 906 "SB03OT.f"
		    k1 = ksize + 1;
#line 907 "SB03OT.f"
		    k2 = ksize + k1;

/*                 Form the right-hand side in DWORK( 1 ),..., */
/*                 DWORK( k - 1 ). */

#line 912 "SB03OT.f"
		    dcopy_(&ksize, &r__[k * r_dim1 + 1], &c__1, &dwork[1], &
			    c__1);
#line 913 "SB03OT.f"
		    d__1 = -alpha;
#line 913 "SB03OT.f"
		    dscal_(&ksize, &d__1, &dwork[1], &c__1);
#line 914 "SB03OT.f"
		    if (cont) {
#line 915 "SB03OT.f"
			d__1 = -r__[k + k * r_dim1];
#line 915 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k * s_dim1 + 1], &c__1, &
				dwork[1], &c__1);
#line 916 "SB03OT.f"
		    } else {
#line 917 "SB03OT.f"
			d__1 = -s[k + k * s_dim1] * r__[k + k * r_dim1];
#line 917 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k * s_dim1 + 1], &c__1, &
				dwork[1], &c__1);
#line 919 "SB03OT.f"
		    }

/*                 SB03OR solves the Sylvester equations. The solution is */
/*                 overwritten on  DWORK. */

#line 924 "SB03OT.f"
		    sb03or_(discr, ltrans, &ksize, &c__1, &s[s_offset], lds, &
			    s[k + k * s_dim1], &c__1, &dwork[1], &ksize, &
			    scaloc, info);
#line 926 "SB03OT.f"
		    infom = max(*info,infom);
#line 927 "SB03OT.f"
		    if (scaloc != 1.) {

#line 929 "SB03OT.f"
			i__1 = *n;
#line 929 "SB03OT.f"
			for (j = 1; j <= i__1; ++j) {
#line 930 "SB03OT.f"
			    dscal_(&j, &scaloc, &r__[j * r_dim1 + 1], &c__1);
#line 931 "SB03OT.f"
/* L150: */
#line 931 "SB03OT.f"
			}

#line 933 "SB03OT.f"
			*scale *= scaloc;
#line 934 "SB03OT.f"
		    }

/*                 Copy the solution into the next  ( k - 1 ) elements */
/*                 of  DWORK,  copy the solution back into  R  and copy */
/*                 the column of  R  back into  DWORK. */

#line 940 "SB03OT.f"
		    dcopy_(&ksize, &dwork[1], &c__1, &dwork[k1], &c__1);
#line 941 "SB03OT.f"
		    dswap_(&ksize, &dwork[1], &c__1, &r__[k * r_dim1 + 1], &
			    c__1);

/*                 Now form the matrix  Rhat  of equation corresponding */
/*                 to (5.15) or (10.17), first computing  y  in  DWORK, */
/*                 and then updating  R1. */

#line 947 "SB03OT.f"
		    if (cont) {
#line 948 "SB03OT.f"
			d__1 = -alpha;
#line 948 "SB03OT.f"
			daxpy_(&ksize, &d__1, &dwork[k1], &c__1, &dwork[1], &
				c__1);
#line 949 "SB03OT.f"
		    } else {

/*                    First form  lambda( 1 )*r  and then add in */
/*                    alpha*u11*s. */

#line 954 "SB03OT.f"
			d__1 = -s[k + k * s_dim1];
#line 954 "SB03OT.f"
			dscal_(&ksize, &d__1, &dwork[1], &c__1);
#line 955 "SB03OT.f"
			d__1 = alpha * r__[k + k * r_dim1];
#line 955 "SB03OT.f"
			daxpy_(&ksize, &d__1, &s[k * s_dim1 + 1], &c__1, &
				dwork[1], &c__1);

/*                    Now form  alpha*S1*u,  first multiplying by the */
/*                    sub-diagonal of  S1  and then the triangular part */
/*                    of  S1,  and add the result in DWORK. */

#line 962 "SB03OT.f"
			i__1 = ksize;
#line 962 "SB03OT.f"
			for (j = 2; j <= i__1; ++j) {
#line 963 "SB03OT.f"
			    if (s[j + (j - 1) * s_dim1] != 0.) {
#line 963 "SB03OT.f"
				dwork[j] = alpha * s[j + (j - 1) * s_dim1] * 
					dwork[k1 + j - 2] + dwork[j];
#line 963 "SB03OT.f"
			    }
#line 965 "SB03OT.f"
/* L160: */
#line 965 "SB03OT.f"
			}

#line 967 "SB03OT.f"
			dtrmv_("Upper", "No transpose", "Non-unit", &ksize, &
				s[s_offset], lds, &dwork[k1], &c__1, (ftnlen)
				5, (ftnlen)12, (ftnlen)8);
#line 969 "SB03OT.f"
			daxpy_(&ksize, &alpha, &dwork[k1], &c__1, &dwork[1], &
				c__1);
#line 970 "SB03OT.f"
		    }
#line 971 "SB03OT.f"
		    mb04nd_("Full", &ksize, &c__0, &c__1, &r__[r_offset], ldr,
			     &dwork[1], &ksize, &dwork[1], &c__1, &dwork[1], &
			    c__1, &dwork[k2], &dwork[k1], (ftnlen)4);
#line 974 "SB03OT.f"
		}
#line 975 "SB03OT.f"
	    }
#line 976 "SB03OT.f"
	    goto L90;
#line 977 "SB03OT.f"
	}
/*        END WHILE 90 */

#line 980 "SB03OT.f"
    }
#line 981 "SB03OT.f"
    *info = infom;
#line 982 "SB03OT.f"
    return 0;
/* *** Last line of SB03OT *** */
} /* sb03ot_ */

