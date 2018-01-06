#line 1 "SB03OR.f"
/* SB03OR.f -- translated by f2c (version 20100827).
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

#line 1 "SB03OR.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_false = FALSE_;
static integer c__2 = 2;

/* Subroutine */ int sb03or_(logical *discr, logical *ltrans, integer *n, 
	integer *m, doublereal *s, integer *lds, doublereal *a, integer *lda, 
	doublereal *c__, integer *ldc, doublereal *scale, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, s_dim1, s_offset, i__1, i__2;

    /* Local variables */
    static integer j, l;
    static doublereal x[4]	/* was [2][2] */;
    static integer l1, l2;
    static doublereal g11, g12, g21, g22;
    static integer dl;
    static doublereal at[4]	/* was [2][2] */, vec[4]	/* was [2][2] 
	    */;
    static integer l2p1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer isgn;
    static logical tbyt;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer infom;
    extern /* Subroutine */ int sb04px_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer lnext;
    static doublereal xnorm;
    extern /* Subroutine */ int dlasy2_(logical *, logical *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal scaloc;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To compute the solution of the Sylvester equations */

/*        op(S)'*X + X*op(A) = scale*C, if DISCR = .FALSE.  or */

/*        op(S)'*X*op(A) - X = scale*C, if DISCR = .TRUE. */

/*     where op(K) = K or K' (i.e., the transpose of the matrix K), S is */
/*     an N-by-N block upper triangular matrix with one-by-one and */
/*     two-by-two blocks on the diagonal, A is an M-by-M matrix (M = 1 or */
/*     M = 2), X and C are each N-by-M matrices, and scale is an output */
/*     scale factor, set less than or equal to 1 to avoid overflow in X. */
/*     The solution X is overwritten on C. */

/*     SB03OR  is a service routine for the Lyapunov solver  SB03OT. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DISCR   LOGICAL */
/*             Specifies the equation to be solved: */
/*             = .FALSE.:  op(S)'*X + X*op(A) = scale*C; */
/*             = .TRUE. :  op(S)'*X*op(A) - X = scale*C. */

/*     LTRANS  LOGICAL */
/*             Specifies the form of op(K) to be used, as follows: */
/*             = .FALSE.:  op(K) = K    (No transpose); */
/*             = .TRUE. :  op(K) = K**T (Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix  S  and also the number of rows of */
/*             matrices  X and C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The order of the matrix  A  and also the number of columns */
/*             of matrices  X and C.  M = 1 or M = 2. */

/*     S       (input) DOUBLE PRECISION array, dimension (LDS,N) */
/*             The leading  N-by-N  upper Hessenberg part of the array  S */
/*             must contain the block upper triangular matrix. The */
/*             elements below the upper Hessenberg part of the array  S */
/*             are not referenced.  The array  S  must not contain */
/*             diagonal blocks larger than two-by-two and the two-by-two */
/*             blocks must only correspond to complex conjugate pairs of */
/*             eigenvalues, not to real eigenvalues. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,N). */

/*     A       (input) DOUBLE PRECISION array, dimension (LDS,M) */
/*             The leading  M-by-M  part of this array must contain a */
/*             given matrix, where M = 1 or M = 2. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= M. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,M) */
/*             On entry, C must contain an N-by-M matrix, where M = 1 or */
/*             M = 2. */
/*             On exit, C contains the N-by-M matrix X, the solution of */
/*             the Sylvester equation. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if DISCR = .FALSE., and S and -A have common */
/*                   eigenvalues, or if DISCR = .TRUE., and S and A have */
/*                   eigenvalues whose product is equal to unity; */
/*                   a solution has been computed using slightly */
/*                   perturbed values. */

/*     METHOD */

/*     The LAPACK scheme for solving Sylvester equations is adapted. */

/*     REFERENCES */

/*     [1] Hammarling, S.J. */
/*         Numerical solution of the stable, non-negative definite */
/*         Lyapunov equation. */
/*         IMA J. Num. Anal., 2, pp. 303-325, 1982. */

/*     NUMERICAL ASPECTS */
/*                               2 */
/*     The algorithm requires 0(N M) operations and is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Supersedes Release 2.0 routines SB03CW and SB03CX by */
/*     Sven Hammarling, NAG Ltd, United Kingdom, Oct. 1986. */
/*     Partly based on routine PLYAP4 by A. Varga, University of Bochum, */
/*     May 1992. */

/*     REVISIONS */

/*     December 1997, April 1998, May 1999, April 2000. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 167 "SB03OR.f"
    /* Parameter adjustments */
#line 167 "SB03OR.f"
    s_dim1 = *lds;
#line 167 "SB03OR.f"
    s_offset = 1 + s_dim1;
#line 167 "SB03OR.f"
    s -= s_offset;
#line 167 "SB03OR.f"
    a_dim1 = *lda;
#line 167 "SB03OR.f"
    a_offset = 1 + a_dim1;
#line 167 "SB03OR.f"
    a -= a_offset;
#line 167 "SB03OR.f"
    c_dim1 = *ldc;
#line 167 "SB03OR.f"
    c_offset = 1 + c_dim1;
#line 167 "SB03OR.f"
    c__ -= c_offset;
#line 167 "SB03OR.f"

#line 167 "SB03OR.f"
    /* Function Body */
#line 167 "SB03OR.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 171 "SB03OR.f"
    if (*n < 0) {
#line 172 "SB03OR.f"
	*info = -3;
#line 173 "SB03OR.f"
    } else if (! (*m == 1 || *m == 2)) {
#line 174 "SB03OR.f"
	*info = -4;
#line 175 "SB03OR.f"
    } else if (*lds < max(1,*n)) {
#line 176 "SB03OR.f"
	*info = -6;
#line 177 "SB03OR.f"
    } else if (*lda < *m) {
#line 178 "SB03OR.f"
	*info = -8;
#line 179 "SB03OR.f"
    } else if (*ldc < max(1,*n)) {
#line 180 "SB03OR.f"
	*info = -10;
#line 181 "SB03OR.f"
    }

#line 183 "SB03OR.f"
    if (*info != 0) {

/*        Error return. */

#line 187 "SB03OR.f"
	i__1 = -(*info);
#line 187 "SB03OR.f"
	xerbla_("SB03OR", &i__1, (ftnlen)6);
#line 188 "SB03OR.f"
	return 0;
#line 189 "SB03OR.f"
    }

#line 191 "SB03OR.f"
    *scale = 1.;

/*     Quick return if possible. */

#line 195 "SB03OR.f"
    if (*n == 0) {
#line 195 "SB03OR.f"
	return 0;
#line 195 "SB03OR.f"
    }

#line 198 "SB03OR.f"
    isgn = 1;
#line 199 "SB03OR.f"
    tbyt = *m == 2;
#line 200 "SB03OR.f"
    infom = 0;

/*     Construct A'. */

#line 204 "SB03OR.f"
    at[0] = a[a_dim1 + 1];
#line 205 "SB03OR.f"
    if (tbyt) {
#line 206 "SB03OR.f"
	at[2] = a[a_dim1 + 2];
#line 207 "SB03OR.f"
	at[1] = a[(a_dim1 << 1) + 1];
#line 208 "SB03OR.f"
	at[3] = a[(a_dim1 << 1) + 2];
#line 209 "SB03OR.f"
    }

#line 211 "SB03OR.f"
    if (*ltrans) {

/*        Start row loop (index = L). */
/*        L1 (L2) : row index of the first (last) row of X(L). */

#line 216 "SB03OR.f"
	lnext = *n;

#line 218 "SB03OR.f"
	for (l = *n; l >= 1; --l) {
#line 219 "SB03OR.f"
	    if (l > lnext) {
#line 219 "SB03OR.f"
		goto L20;
#line 219 "SB03OR.f"
	    }
#line 221 "SB03OR.f"
	    l1 = l;
#line 222 "SB03OR.f"
	    l2 = l;
#line 223 "SB03OR.f"
	    if (l > 1) {
#line 224 "SB03OR.f"
		if (s[l + (l - 1) * s_dim1] != 0.) {
#line 224 "SB03OR.f"
		    --l1;
#line 224 "SB03OR.f"
		}
#line 226 "SB03OR.f"
		lnext = l1 - 1;
#line 227 "SB03OR.f"
	    }
#line 228 "SB03OR.f"
	    dl = l2 - l1 + 1;
/* Computing MIN */
#line 229 "SB03OR.f"
	    i__1 = l2 + 1;
#line 229 "SB03OR.f"
	    l2p1 = min(i__1,*n);

#line 231 "SB03OR.f"
	    if (*discr) {

/*              Solve  S*X*A' - X = scale*C. */

/*              The L-th block of X is determined from */

/*              S(L,L)*X(L)*A' - X(L) = C(L) - R(L), */

/*              where */

/*                      N */
/*              R(L) = SUM [S(L,J)*X(J)] * A' . */
/*                    J=L+1 */

#line 245 "SB03OR.f"
		i__1 = *n - l2;
#line 245 "SB03OR.f"
		g11 = -ddot_(&i__1, &s[l1 + l2p1 * s_dim1], lds, &c__[l2p1 + 
			c_dim1], &c__1);
#line 246 "SB03OR.f"
		if (tbyt) {
#line 247 "SB03OR.f"
		    i__1 = *n - l2;
#line 247 "SB03OR.f"
		    g12 = -ddot_(&i__1, &s[l1 + l2p1 * s_dim1], lds, &c__[
			    l2p1 + (c_dim1 << 1)], &c__1);
#line 249 "SB03OR.f"
		    vec[0] = c__[l1 + c_dim1] + g11 * at[0] + g12 * at[1];
#line 250 "SB03OR.f"
		    vec[2] = c__[l1 + (c_dim1 << 1)] + g11 * at[2] + g12 * at[
			    3];
#line 251 "SB03OR.f"
		} else {
#line 252 "SB03OR.f"
		    vec[0] = c__[l1 + c_dim1] + g11 * at[0];
#line 253 "SB03OR.f"
		}
#line 254 "SB03OR.f"
		if (dl != 1) {
#line 255 "SB03OR.f"
		    i__1 = *n - l2;
#line 255 "SB03OR.f"
		    g21 = -ddot_(&i__1, &s[l2 + l2p1 * s_dim1], lds, &c__[
			    l2p1 + c_dim1], &c__1);
#line 257 "SB03OR.f"
		    if (tbyt) {
#line 258 "SB03OR.f"
			i__1 = *n - l2;
#line 258 "SB03OR.f"
			g22 = -ddot_(&i__1, &s[l2 + l2p1 * s_dim1], lds, &c__[
				l2p1 + (c_dim1 << 1)], &c__1);
#line 260 "SB03OR.f"
			vec[1] = c__[l2 + c_dim1] + g21 * at[0] + g22 * at[1];
#line 262 "SB03OR.f"
			vec[3] = c__[l2 + (c_dim1 << 1)] + g21 * at[2] + g22 *
				 at[3];
#line 264 "SB03OR.f"
		    } else {
#line 265 "SB03OR.f"
			vec[1] = c__[l2 + c_dim1] + g21 * at[0];
#line 266 "SB03OR.f"
		    }
#line 267 "SB03OR.f"
		}
#line 268 "SB03OR.f"
		i__1 = -isgn;
#line 268 "SB03OR.f"
		sb04px_(&c_false, &c_false, &i__1, &dl, m, &s[l1 + l1 * 
			s_dim1], lds, at, &c__2, vec, &c__2, &scaloc, x, &
			c__2, &xnorm, info);
#line 271 "SB03OR.f"
	    } else {

/*              Solve  S*X + X*A' = scale*C. */

/*              The L-th block of X is determined from */

/*              S(L,L)*X(L) + X(L)*A' = C(L) - R(L), */

/*              where */
/*                       N */
/*              R(L) =  SUM S(L,J)*X(J) . */
/*                     J=L+1 */

#line 284 "SB03OR.f"
		i__1 = *n - l2;
#line 284 "SB03OR.f"
		vec[0] = c__[l1 + c_dim1] - ddot_(&i__1, &s[l1 + l2p1 * 
			s_dim1], lds, &c__[l2p1 + c_dim1], &c__1);
#line 287 "SB03OR.f"
		if (tbyt) {
#line 287 "SB03OR.f"
		    i__1 = *n - l2;
#line 287 "SB03OR.f"
		    vec[2] = c__[l1 + (c_dim1 << 1)] - ddot_(&i__1, &s[l1 + 
			    l2p1 * s_dim1], lds, &c__[l2p1 + (c_dim1 << 1)], &
			    c__1);
#line 287 "SB03OR.f"
		}

#line 292 "SB03OR.f"
		if (dl != 1) {
#line 293 "SB03OR.f"
		    i__1 = *n - l2;
#line 293 "SB03OR.f"
		    vec[1] = c__[l2 + c_dim1] - ddot_(&i__1, &s[l2 + l2p1 * 
			    s_dim1], lds, &c__[l2p1 + c_dim1], &c__1);
#line 296 "SB03OR.f"
		    if (tbyt) {
#line 296 "SB03OR.f"
			i__1 = *n - l2;
#line 296 "SB03OR.f"
			vec[3] = c__[l2 + (c_dim1 << 1)] - ddot_(&i__1, &s[l2 
				+ l2p1 * s_dim1], lds, &c__[l2p1 + (c_dim1 << 
				1)], &c__1);
#line 296 "SB03OR.f"
		    }
#line 300 "SB03OR.f"
		}
#line 301 "SB03OR.f"
		dlasy2_(&c_false, &c_false, &isgn, &dl, m, &s[l1 + l1 * 
			s_dim1], lds, at, &c__2, vec, &c__2, &scaloc, x, &
			c__2, &xnorm, info);
#line 304 "SB03OR.f"
	    }
#line 305 "SB03OR.f"
	    infom = max(*info,infom);
#line 306 "SB03OR.f"
	    if (scaloc != 1.) {

#line 308 "SB03OR.f"
		i__1 = *m;
#line 308 "SB03OR.f"
		for (j = 1; j <= i__1; ++j) {
#line 309 "SB03OR.f"
		    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 310 "SB03OR.f"
/* L10: */
#line 310 "SB03OR.f"
		}

#line 312 "SB03OR.f"
		*scale *= scaloc;
#line 313 "SB03OR.f"
	    }
#line 314 "SB03OR.f"
	    c__[l1 + c_dim1] = x[0];
#line 315 "SB03OR.f"
	    if (tbyt) {
#line 315 "SB03OR.f"
		c__[l1 + (c_dim1 << 1)] = x[2];
#line 315 "SB03OR.f"
	    }
#line 316 "SB03OR.f"
	    if (dl != 1) {
#line 317 "SB03OR.f"
		c__[l2 + c_dim1] = x[1];
#line 318 "SB03OR.f"
		if (tbyt) {
#line 318 "SB03OR.f"
		    c__[l2 + (c_dim1 << 1)] = x[3];
#line 318 "SB03OR.f"
		}
#line 319 "SB03OR.f"
	    }
#line 320 "SB03OR.f"
L20:
#line 320 "SB03OR.f"
	    ;
#line 320 "SB03OR.f"
	}

#line 322 "SB03OR.f"
    } else {

/*        Start row loop (index = L). */
/*        L1 (L2) : row index of the first (last) row of X(L). */

#line 327 "SB03OR.f"
	lnext = 1;

#line 329 "SB03OR.f"
	i__1 = *n;
#line 329 "SB03OR.f"
	for (l = 1; l <= i__1; ++l) {
#line 330 "SB03OR.f"
	    if (l < lnext) {
#line 330 "SB03OR.f"
		goto L40;
#line 330 "SB03OR.f"
	    }
#line 332 "SB03OR.f"
	    l1 = l;
#line 333 "SB03OR.f"
	    l2 = l;
#line 334 "SB03OR.f"
	    if (l < *n) {
#line 335 "SB03OR.f"
		if (s[l + 1 + l * s_dim1] != 0.) {
#line 335 "SB03OR.f"
		    ++l2;
#line 335 "SB03OR.f"
		}
#line 337 "SB03OR.f"
		lnext = l2 + 1;
#line 338 "SB03OR.f"
	    }
#line 339 "SB03OR.f"
	    dl = l2 - l1 + 1;

#line 341 "SB03OR.f"
	    if (*discr) {

/*              Solve  A'*X'*S - X' = scale*C'. */

/*              The L-th block of X is determined from */

/*              A'*X(L)'*S(L,L) - X(L)' = C(L)' - R(L), */

/*              where */

/*                          L-1 */
/*              R(L) = A' * SUM [X(J)'*S(J,L)] . */
/*                          J=1 */

#line 355 "SB03OR.f"
		i__2 = l1 - 1;
#line 355 "SB03OR.f"
		g11 = -ddot_(&i__2, &c__[c_offset], &c__1, &s[l1 * s_dim1 + 1]
			, &c__1);
#line 356 "SB03OR.f"
		if (tbyt) {
#line 357 "SB03OR.f"
		    i__2 = l1 - 1;
#line 357 "SB03OR.f"
		    g21 = -ddot_(&i__2, &c__[(c_dim1 << 1) + 1], &c__1, &s[l1 
			    * s_dim1 + 1], &c__1);
#line 358 "SB03OR.f"
		    vec[0] = c__[l1 + c_dim1] + at[0] * g11 + at[2] * g21;
#line 359 "SB03OR.f"
		    vec[1] = c__[l1 + (c_dim1 << 1)] + at[1] * g11 + at[3] * 
			    g21;
#line 360 "SB03OR.f"
		} else {
#line 361 "SB03OR.f"
		    vec[0] = c__[l1 + c_dim1] + at[0] * g11;
#line 362 "SB03OR.f"
		}
#line 363 "SB03OR.f"
		if (dl != 1) {
#line 364 "SB03OR.f"
		    i__2 = l1 - 1;
#line 364 "SB03OR.f"
		    g12 = -ddot_(&i__2, &c__[c_offset], &c__1, &s[l2 * s_dim1 
			    + 1], &c__1);
#line 365 "SB03OR.f"
		    if (tbyt) {
#line 366 "SB03OR.f"
			i__2 = l1 - 1;
#line 366 "SB03OR.f"
			g22 = -ddot_(&i__2, &c__[(c_dim1 << 1) + 1], &c__1, &
				s[l2 * s_dim1 + 1], &c__1);
#line 367 "SB03OR.f"
			vec[2] = c__[l2 + c_dim1] + at[0] * g12 + at[2] * g22;
#line 369 "SB03OR.f"
			vec[3] = c__[l2 + (c_dim1 << 1)] + at[1] * g12 + at[3]
				 * g22;
#line 371 "SB03OR.f"
		    } else {
#line 372 "SB03OR.f"
			vec[2] = c__[l2 + c_dim1] + at[0] * g12;
#line 373 "SB03OR.f"
		    }
#line 374 "SB03OR.f"
		}
#line 375 "SB03OR.f"
		i__2 = -isgn;
#line 375 "SB03OR.f"
		sb04px_(&c_false, &c_false, &i__2, m, &dl, at, &c__2, &s[l1 + 
			l1 * s_dim1], lds, vec, &c__2, &scaloc, x, &c__2, &
			xnorm, info);
#line 378 "SB03OR.f"
	    } else {

/*              Solve  A'*X' + X'*S = scale*C'. */

/*              The L-th block of X is determined from */

/*              A'*X(L)' + X(L)'*S(L,L) = C(L)' - R(L), */

/*              where */
/*                     L-1 */
/*              R(L) = SUM [X(J)'*S(J,L)]. */
/*                     J=1 */

#line 391 "SB03OR.f"
		i__2 = l1 - 1;
#line 391 "SB03OR.f"
		vec[0] = c__[l1 + c_dim1] - ddot_(&i__2, &c__[c_offset], &
			c__1, &s[l1 * s_dim1 + 1], &c__1);
#line 393 "SB03OR.f"
		if (tbyt) {
#line 393 "SB03OR.f"
		    i__2 = l1 - 1;
#line 393 "SB03OR.f"
		    vec[1] = c__[l1 + (c_dim1 << 1)] - ddot_(&i__2, &c__[(
			    c_dim1 << 1) + 1], &c__1, &s[l1 * s_dim1 + 1], &
			    c__1);
#line 393 "SB03OR.f"
		}

#line 397 "SB03OR.f"
		if (dl != 1) {
#line 398 "SB03OR.f"
		    i__2 = l1 - 1;
#line 398 "SB03OR.f"
		    vec[2] = c__[l2 + c_dim1] - ddot_(&i__2, &c__[c_offset], &
			    c__1, &s[l2 * s_dim1 + 1], &c__1);
#line 400 "SB03OR.f"
		    if (tbyt) {
#line 400 "SB03OR.f"
			i__2 = l1 - 1;
#line 400 "SB03OR.f"
			vec[3] = c__[l2 + (c_dim1 << 1)] - ddot_(&i__2, &c__[(
				c_dim1 << 1) + 1], &c__1, &s[l2 * s_dim1 + 1],
				 &c__1);
#line 400 "SB03OR.f"
		    }
#line 403 "SB03OR.f"
		}
#line 404 "SB03OR.f"
		dlasy2_(&c_false, &c_false, &isgn, m, &dl, at, &c__2, &s[l1 + 
			l1 * s_dim1], lds, vec, &c__2, &scaloc, x, &c__2, &
			xnorm, info);
#line 407 "SB03OR.f"
	    }
#line 408 "SB03OR.f"
	    infom = max(*info,infom);
#line 409 "SB03OR.f"
	    if (scaloc != 1.) {

#line 411 "SB03OR.f"
		i__2 = *m;
#line 411 "SB03OR.f"
		for (j = 1; j <= i__2; ++j) {
#line 412 "SB03OR.f"
		    dscal_(n, &scaloc, &c__[j * c_dim1 + 1], &c__1);
#line 413 "SB03OR.f"
/* L30: */
#line 413 "SB03OR.f"
		}

#line 415 "SB03OR.f"
		*scale *= scaloc;
#line 416 "SB03OR.f"
	    }
#line 417 "SB03OR.f"
	    c__[l1 + c_dim1] = x[0];
#line 418 "SB03OR.f"
	    if (tbyt) {
#line 418 "SB03OR.f"
		c__[l1 + (c_dim1 << 1)] = x[1];
#line 418 "SB03OR.f"
	    }
#line 419 "SB03OR.f"
	    if (dl != 1) {
#line 420 "SB03OR.f"
		c__[l2 + c_dim1] = x[2];
#line 421 "SB03OR.f"
		if (tbyt) {
#line 421 "SB03OR.f"
		    c__[l2 + (c_dim1 << 1)] = x[3];
#line 421 "SB03OR.f"
		}
#line 422 "SB03OR.f"
	    }
#line 423 "SB03OR.f"
L40:
#line 423 "SB03OR.f"
	    ;
#line 423 "SB03OR.f"
	}
#line 424 "SB03OR.f"
    }

#line 426 "SB03OR.f"
    *info = infom;
#line 427 "SB03OR.f"
    return 0;
/* *** Last line of SB03OR *** */
} /* sb03or_ */

