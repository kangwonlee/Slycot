#line 1 "SB03OU.f"
/* SB03OU.f -- translated by f2c (version 20100827).
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

#line 1 "SB03OU.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = 0.;

/* Subroutine */ int sb03ou_(logical *discr, logical *ltrans, integer *n, 
	integer *m, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *tau, doublereal *u, integer *ldu, doublereal *scale, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, u_dim1, u_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, k, l, mn;
    extern /* Subroutine */ int sb03ot_(logical *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *), dcopy_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), dgeqrf_(integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dgerqf_(integer *, integer *, doublereal *, integer *,
	     doublereal *, doublereal *, integer *, integer *), dlacpy_(char *
	    , integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer wrkopt;


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
/*        op(A)'*X + X*op(A) = -scale *op(B)'*op(B)                   (1) */

/*     or the convergent non-negative definite discrete-time Lyapunov */
/*     equation */
/*                                   2 */
/*        op(A)'*X*op(A) - X = -scale *op(B)'*op(B)                   (2) */

/*     where op(K) = K or K' (i.e., the transpose of the matrix K), A is */
/*     an N-by-N matrix in real Schur form, op(B) is an M-by-N matrix, */
/*     U is an upper triangular matrix containing the Cholesky factor of */
/*     the solution matrix X, X = op(U)'*op(U), and scale is an output */
/*     scale factor, set less than or equal to 1 to avoid overflow in X. */
/*     If matrix B has full rank then the solution matrix X will be */
/*     positive-definite and hence the Cholesky factor U will be */
/*     nonsingular, but if B is rank deficient then X may only be */
/*     positive semi-definite and U will be singular. */

/*     In the case of equation (1) the matrix A must be stable (that */
/*     is, all the eigenvalues of A must have negative real parts), */
/*     and for equation (2) the matrix A must be convergent (that is, */
/*     all the eigenvalues of A must lie inside the unit circle). */

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
/*             The order of the matrix A and the number of columns in */
/*             matrix op(B).  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of rows in matrix op(B).  M >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N upper Hessenberg part of this array */
/*             must contain a real Schur form matrix S. The elements */
/*             below the upper Hessenberg part of the array A are not */
/*             referenced. The 2-by-2 blocks must only correspond to */
/*             complex conjugate pairs of eigenvalues (not to real */
/*             eigenvalues). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             if LTRANS = .FALSE., and dimension (LDB,M), if */
/*             LTRANS = .TRUE.. */
/*             On entry, if LTRANS = .FALSE., the leading M-by-N part of */
/*             this array must contain the coefficient matrix B of the */
/*             equation. */
/*             On entry, if LTRANS = .TRUE., the leading N-by-M part of */
/*             this array must contain the coefficient matrix B of the */
/*             equation. */
/*             On exit, if LTRANS = .FALSE., the leading */
/*             MIN(M,N)-by-MIN(M,N) upper triangular part of this array */
/*             contains the upper triangular matrix R (as defined in */
/*             METHOD), and the M-by-MIN(M,N) strictly lower triangular */
/*             part together with the elements of the array TAU are */
/*             overwritten by details of the matrix P (also defined in */
/*             METHOD). When M < N, columns (M+1),...,N of the array B */
/*             are overwritten by the matrix Z (see METHOD). */
/*             On exit, if LTRANS = .TRUE., the leading */
/*             MIN(M,N)-by-MIN(M,N) upper triangular part of */
/*             B(1:N,M-N+1), if M >= N, or of B(N-M+1:N,1:M), if M < N, */
/*             contains the upper triangular matrix R (as defined in */
/*             METHOD), and the remaining elements (below the diagonal */
/*             of R) together with the elements of the array TAU are */
/*             overwritten by details of the matrix P (also defined in */
/*             METHOD). When M < N, rows 1,...,(N-M) of the array B */
/*             are overwritten by the matrix Z (see METHOD). */

/*     LDB     INTEGER */
/*             The leading dimension of array B. */
/*             LDB >= MAX(1,M), if LTRANS = .FALSE., */
/*             LDB >= MAX(1,N), if LTRANS = .TRUE.. */

/*     TAU     (output) DOUBLE PRECISION array of dimension (MIN(N,M)) */
/*             This array contains the scalar factors of the elementary */
/*             reflectors defining the matrix P. */

/*     U       (output) DOUBLE PRECISION array of dimension (LDU,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the Cholesky factor of the solution matrix X of */
/*             the problem, X = op(U)'*op(U). */
/*             The array U may be identified with B in the calling */
/*             statement, if B is properly dimensioned, and the */
/*             intermediate results returned in B are not needed. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,N). */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, or INFO = 1, DWORK(1) returns the */
/*             optimal value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. LDWORK >= MAX(1,4*N). */
/*             For optimum performance LDWORK should sometimes be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the Lyapunov equation is (nearly) singular */
/*                   (warning indicator); */
/*                   if DISCR = .FALSE., this means that while the matrix */
/*                   A has computed eigenvalues with negative real parts, */
/*                   it is only just stable in the sense that small */
/*                   perturbations in A can make one or more of the */
/*                   eigenvalues have a non-negative real part; */
/*                   if DISCR = .TRUE., this means that while the matrix */
/*                   A has computed eigenvalues inside the unit circle, */
/*                   it is nevertheless only just convergent, in the */
/*                   sense that small perturbations in A can make one or */
/*                   more of the eigenvalues lie outside the unit circle; */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrix A is unchanged); */
/*             = 2:  if matrix A is not stable (that is, one or more of */
/*                   the eigenvalues of A has a non-negative real part), */
/*                   if DISCR = .FALSE., or not convergent (that is, one */
/*                   or more of the eigenvalues of A lies outside the */
/*                   unit circle), if DISCR = .TRUE.; */
/*             = 3:  if matrix A has two or more consecutive non-zero */
/*                   elements on the first sub-diagonal, so that there is */
/*                   a block larger than 2-by-2 on the diagonal; */
/*             = 4:  if matrix A has a 2-by-2 diagonal block with real */
/*                   eigenvalues instead of a complex conjugate pair. */

/*     METHOD */

/*     The method used by the routine is based on the Bartels and */
/*     Stewart method [1], except that it finds the upper triangular */
/*     matrix U directly without first finding X and without the need */
/*     to form the normal matrix op(B)'*op(B) [2]. */

/*     If LTRANS = .FALSE., the matrix B is factored as */

/*        B = P ( R ),  M >= N,   B = P ( R  Z ),  M < N, */
/*              ( 0 ) */

/*     (QR factorization), where P is an M-by-M orthogonal matrix and */
/*     R is a square upper triangular matrix. */

/*     If LTRANS = .TRUE., the matrix B is factored as */

/*        B = ( 0  R ) P,  M >= N,  B = ( Z ) P,  M < N, */
/*                                      ( R ) */

/*     (RQ factorization), where P is an M-by-M orthogonal matrix and */
/*     R is a square upper triangular matrix. */

/*     These factorizations are used to solve the continuous-time */
/*     Lyapunov equation in the canonical form */
/*                                                        2 */
/*       op(A)'*op(U)'*op(U) + op(U)'*op(U)*op(A) = -scale *op(F)'*op(F), */

/*     or the discrete-time Lyapunov equation in the canonical form */
/*                                                        2 */
/*       op(A)'*op(U)'*op(U)*op(A) - op(U)'*op(U) = -scale *op(F)'*op(F), */

/*     where U and F are N-by-N upper triangular matrices, and */

/*        F = R,                                  if M >= N, or */

/*        F = ( R ),    if LTRANS = .FALSE.,  or */
/*            ( 0 ) */

/*        F = ( 0  Z ), if LTRANS = .TRUE.,       if M < N. */
/*            ( 0  R ) */

/*     The canonical equation is solved for U. */

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

/*     The Lyapunov equation may be very ill-conditioned. In particular, */
/*     if A is only just stable (or convergent) then the Lyapunov */
/*     equation will be ill-conditioned. "Large" elements in U relative */
/*     to those of A and B, or a "small" value for scale, are symptoms */
/*     of ill-conditioning. A condition estimate can be computed using */
/*     SLICOT Library routine SB03MD. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Supersedes Release 2.0 routine SB03CZ by Sven Hammarling, */
/*     NAG Ltd, United Kingdom. */
/*     Partly based on routine PLYAPS by A. Varga, University of Bochum, */
/*     May 1992. */

/*     REVISIONS */

/*     Dec. 1997, April 1998, May 1999. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 280 "SB03OU.f"
    /* Parameter adjustments */
#line 280 "SB03OU.f"
    a_dim1 = *lda;
#line 280 "SB03OU.f"
    a_offset = 1 + a_dim1;
#line 280 "SB03OU.f"
    a -= a_offset;
#line 280 "SB03OU.f"
    b_dim1 = *ldb;
#line 280 "SB03OU.f"
    b_offset = 1 + b_dim1;
#line 280 "SB03OU.f"
    b -= b_offset;
#line 280 "SB03OU.f"
    --tau;
#line 280 "SB03OU.f"
    u_dim1 = *ldu;
#line 280 "SB03OU.f"
    u_offset = 1 + u_dim1;
#line 280 "SB03OU.f"
    u -= u_offset;
#line 280 "SB03OU.f"
    --dwork;
#line 280 "SB03OU.f"

#line 280 "SB03OU.f"
    /* Function Body */
#line 280 "SB03OU.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 284 "SB03OU.f"
    if (*n < 0) {
#line 285 "SB03OU.f"
	*info = -3;
#line 286 "SB03OU.f"
    } else if (*m < 0) {
#line 287 "SB03OU.f"
	*info = -4;
#line 288 "SB03OU.f"
    } else if (*lda < max(1,*n)) {
#line 289 "SB03OU.f"
	*info = -6;
#line 290 "SB03OU.f"
    } else if (*ldb < max(1,*m) && ! (*ltrans) || *ldb < max(1,*n) && *ltrans)
	     {
#line 292 "SB03OU.f"
	*info = -8;
#line 293 "SB03OU.f"
    } else if (*ldu < max(1,*n)) {
#line 294 "SB03OU.f"
	*info = -11;
#line 295 "SB03OU.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 295 "SB03OU.f"
	i__1 = 1, i__2 = *n << 2;
#line 295 "SB03OU.f"
	if (*ldwork < max(i__1,i__2)) {
#line 296 "SB03OU.f"
	    *info = -14;
#line 297 "SB03OU.f"
	}
#line 297 "SB03OU.f"
    }

#line 299 "SB03OU.f"
    if (*info != 0) {

/*        Error return. */

#line 303 "SB03OU.f"
	i__1 = -(*info);
#line 303 "SB03OU.f"
	xerbla_("SB03OU", &i__1, (ftnlen)6);
#line 304 "SB03OU.f"
	return 0;
#line 305 "SB03OU.f"
    }

/*     Quick return if possible. */

#line 309 "SB03OU.f"
    mn = min(*n,*m);
#line 310 "SB03OU.f"
    if (mn == 0) {
#line 311 "SB03OU.f"
	*scale = 1.;
#line 312 "SB03OU.f"
	dwork[1] = 1.;
#line 313 "SB03OU.f"
	return 0;
#line 314 "SB03OU.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 322 "SB03OU.f"
    if (*ltrans) {

/*        Case op(K) = K'. */

/*        Perform the RQ factorization of B. */
/*        Workspace: need   N; */
/*                   prefer N*NB. */

#line 330 "SB03OU.f"
	dgerqf_(n, m, &b[b_offset], ldb, &tau[1], &dwork[1], ldwork, info);

/*        The triangular matrix F is constructed in the array U so that */
/*        U can share the same memory as B. */

#line 335 "SB03OU.f"
	if (*m >= *n) {
#line 336 "SB03OU.f"
	    dlacpy_("Upper", &mn, n, &b[(*m - *n + 1) * b_dim1 + 1], ldb, &u[
		    u_offset], ldu, (ftnlen)5);
#line 337 "SB03OU.f"
	} else {

#line 339 "SB03OU.f"
	    for (i__ = *m; i__ >= 1; --i__) {
#line 340 "SB03OU.f"
		i__1 = *n - *m + i__;
#line 340 "SB03OU.f"
		dcopy_(&i__1, &b[i__ * b_dim1 + 1], &c__1, &u[(*n - *m + i__) 
			* u_dim1 + 1], &c__1);
#line 341 "SB03OU.f"
/* L10: */
#line 341 "SB03OU.f"
	    }

#line 343 "SB03OU.f"
	    i__1 = *n - *m;
#line 343 "SB03OU.f"
	    dlaset_("Full", n, &i__1, &c_b8, &c_b8, &u[u_offset], ldu, (
		    ftnlen)4);
#line 344 "SB03OU.f"
	}
#line 345 "SB03OU.f"
    } else {

/*        Case op(K) = K. */

/*        Perform the QR factorization of B. */
/*        Workspace: need   N; */
/*                   prefer N*NB. */

#line 353 "SB03OU.f"
	dgeqrf_(m, n, &b[b_offset], ldb, &tau[1], &dwork[1], ldwork, info);
#line 354 "SB03OU.f"
	dlacpy_("Upper", &mn, n, &b[b_offset], ldb, &u[u_offset], ldu, (
		ftnlen)5);
#line 355 "SB03OU.f"
	if (*m < *n) {
#line 355 "SB03OU.f"
	    i__1 = *n - *m;
#line 355 "SB03OU.f"
	    i__2 = *n - *m;
#line 355 "SB03OU.f"
	    dlaset_("Upper", &i__1, &i__2, &c_b8, &c_b8, &u[*m + 1 + (*m + 1) 
		    * u_dim1], ldu, (ftnlen)5);
#line 355 "SB03OU.f"
	}
#line 358 "SB03OU.f"
    }
#line 359 "SB03OU.f"
    wrkopt = (integer) dwork[1];

/*     Solve the canonical Lyapunov equation */
/*                                                      2 */
/*     op(A)'*op(U)'*op(U) + op(U)'*op(U)*op(A) = -scale *op(F)'*op(F), */

/*     or */
/*                                                      2 */
/*     op(A)'*op(U)'*op(U)*op(A) - op(U)'*op(U) = -scale *op(F)'*op(F) */

/*     for U. */

#line 371 "SB03OU.f"
    sb03ot_(discr, ltrans, n, &a[a_offset], lda, &u[u_offset], ldu, scale, &
	    dwork[1], info);
#line 373 "SB03OU.f"
    if (*info != 0 && *info != 1) {
#line 373 "SB03OU.f"
	return 0;
#line 373 "SB03OU.f"
    }

/*     Make the diagonal elements of U non-negative. */

#line 378 "SB03OU.f"
    if (*ltrans) {

#line 380 "SB03OU.f"
	i__1 = *n;
#line 380 "SB03OU.f"
	for (j = 1; j <= i__1; ++j) {
#line 381 "SB03OU.f"
	    if (u[j + j * u_dim1] < 0.) {

#line 383 "SB03OU.f"
		i__2 = j;
#line 383 "SB03OU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 384 "SB03OU.f"
		    u[i__ + j * u_dim1] = -u[i__ + j * u_dim1];
#line 385 "SB03OU.f"
/* L20: */
#line 385 "SB03OU.f"
		}

#line 387 "SB03OU.f"
	    }
#line 388 "SB03OU.f"
/* L30: */
#line 388 "SB03OU.f"
	}

#line 390 "SB03OU.f"
    } else {
#line 391 "SB03OU.f"
	k = 1;

#line 393 "SB03OU.f"
	i__1 = *n;
#line 393 "SB03OU.f"
	for (j = 1; j <= i__1; ++j) {
#line 394 "SB03OU.f"
	    dwork[k] = u[j + j * u_dim1];
#line 395 "SB03OU.f"
	    l = 1;

#line 397 "SB03OU.f"
	    i__2 = j;
#line 397 "SB03OU.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 398 "SB03OU.f"
		if (dwork[l] < 0.) {
#line 398 "SB03OU.f"
		    u[i__ + j * u_dim1] = -u[i__ + j * u_dim1];
#line 398 "SB03OU.f"
		}
#line 399 "SB03OU.f"
		++l;
#line 400 "SB03OU.f"
/* L40: */
#line 400 "SB03OU.f"
	    }

#line 402 "SB03OU.f"
	    ++k;
#line 403 "SB03OU.f"
/* L50: */
#line 403 "SB03OU.f"
	}

#line 405 "SB03OU.f"
    }

/* Computing MAX */
#line 407 "SB03OU.f"
    i__1 = wrkopt, i__2 = *n << 2;
#line 407 "SB03OU.f"
    dwork[1] = (doublereal) max(i__1,i__2);
#line 408 "SB03OU.f"
    return 0;
/* *** Last line of SB03OU *** */
} /* sb03ou_ */

