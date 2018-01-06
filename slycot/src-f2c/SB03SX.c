#line 1 "SB03SX.f"
/* SB03SX.f -- translated by f2c (version 20100827).
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

#line 1 "SB03SX.f"
/* Table of constant values */

static doublereal c_b24 = 0.;
static doublereal c_b25 = 1.;
static doublereal c_b26 = .5;

/* Subroutine */ int sb03sx_(char *trana, char *uplo, char *lyapun, integer *
	n, doublereal *xanorm, doublereal *t, integer *ldt, doublereal *u, 
	integer *ldu, doublereal *r__, integer *ldr, doublereal *ferr, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen trana_len, ftnlen uplo_len, ftnlen lyapun_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, t_dim1, t_offset, u_dim1, u_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, ij, nn;
    static doublereal est;
    static integer kase;
    static doublereal temp;
    static integer itmp, info2;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen), dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb01ru_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), sb03mx_(char *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen);
    static logical lower;
    static char uplow[1];
    extern /* Subroutine */ int dlacon_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *), xerbla_(char *, integer *, 
	    ftnlen);
    static logical update;
    static char tranat[1];
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);
    static logical notrna;


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

/*     To estimate a forward error bound for the solution X of a real */
/*     discrete-time Lyapunov matrix equation, */

/*            op(A)'*X*op(A) - X = C, */

/*     where op(A) = A or A' (A**T) and C is symmetric (C = C**T). The */
/*     matrix A, the right hand side C, and the solution X are N-by-N. */
/*     An absolute residual matrix, which takes into account the rounding */
/*     errors in forming it, is given in the array R. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     UPLO    CHARACTER*1 */
/*             Specifies which part of the symmetric matrix R is to be */
/*             used, as follows: */
/*             = 'U':  Upper triangular part; */
/*             = 'L':  Lower triangular part. */

/*     LYAPUN  CHARACTER*1 */
/*             Specifies whether or not the original Lyapunov equations */
/*             should be solved, as follows: */
/*             = 'O':  Solve the original Lyapunov equations, updating */
/*                     the right-hand sides and solutions with the */
/*                     matrix U, e.g., X <-- U'*X*U; */
/*             = 'R':  Solve reduced Lyapunov equations only, without */
/*                     updating the right-hand sides and solutions. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and R.  N >= 0. */

/*     XANORM  (input) DOUBLE PRECISION */
/*             The absolute (maximal) norm of the symmetric solution */
/*             matrix X of the Lyapunov equation.  XANORM >= 0. */

/*     T       (input) DOUBLE PRECISION array, dimension (LDT,N) */
/*             The leading N-by-N upper Hessenberg part of this array */
/*             must contain the upper quasi-triangular matrix T in Schur */
/*             canonical form from a Schur factorization of A. */

/*     LDT     INTEGER */
/*             The leading dimension of array T.  LDT >= MAX(1,N). */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,N) */
/*             The leading N-by-N part of this array must contain the */
/*             orthogonal matrix U from a real Schur factorization of A. */
/*             If LYAPUN = 'R', the array U is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of array U. */
/*             LDU >= 1,        if LYAPUN = 'R'; */
/*             LDU >= MAX(1,N), if LYAPUN = 'O'. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR,N) */
/*             On entry, if UPLO = 'U', the leading N-by-N upper */
/*             triangular part of this array must contain the upper */
/*             triangular part of the absolute residual matrix R, with */
/*             bounds on rounding errors added. */
/*             On entry, if UPLO = 'L', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangular part of the absolute residual matrix R, with */
/*             bounds on rounding errors added. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the symmetric absolute residual matrix R (with bounds on */
/*             rounding errors added), fully stored. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     FERR    (output) DOUBLE PRECISION */
/*             An estimated forward error bound for the solution X. */
/*             If XTRUE is the true solution, FERR bounds the magnitude */
/*             of the largest entry in (X - XTRUE) divided by the */
/*             magnitude of the largest entry in X. */
/*             If N = 0 or XANORM = 0, FERR is set to 0, without any */
/*             calculations. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 0,            if N = 0; */
/*             LDWORK >= MAX(3,2*N*N), if N > 0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = N+1:  if T has almost reciprocal eigenvalues; perturbed */
/*                   values were used to solve Lyapunov equations (but */
/*                   the matrix T is unchanged). */

/*     METHOD */

/*     The forward error bound is estimated using a practical error bound */
/*     similar to the one proposed in [1], based on the 1-norm estimator */
/*     in [2]. */

/*     REFERENCES */

/*     [1] Higham, N.J. */
/*         Perturbation theory and backward error for AX-XB=C. */
/*         BIT, vol. 33, pp. 124-136, 1993. */

/*     [2] Higham, N.J. */
/*         FORTRAN codes for estimating the one-norm of a real or */
/*         complex matrix, with applications to condition estimation. */
/*         ACM Trans. Math. Softw., 14, pp. 381-396, 1988. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     FURTHER COMMENTS */

/*     The option LYAPUN = 'R' may occasionally produce slightly worse */
/*     or better estimates, and it is much faster than the option 'O'. */
/*     The routine can be also used as a final step in estimating a */
/*     forward error bound for the solution of a discrete-time algebraic */
/*     matrix Riccati equation. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998. Partly based on DDLSVX (and then SB03SD) by P. Petkov, */
/*     Tech. University of Sofia, March 1998 (and December 1998). */

/*     REVISIONS */

/*     February 6, 1999, V. Sima, Katholieke Univ. Leuven, Belgium. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2004. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form. */

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
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

#line 212 "SB03SX.f"
    /* Parameter adjustments */
#line 212 "SB03SX.f"
    t_dim1 = *ldt;
#line 212 "SB03SX.f"
    t_offset = 1 + t_dim1;
#line 212 "SB03SX.f"
    t -= t_offset;
#line 212 "SB03SX.f"
    u_dim1 = *ldu;
#line 212 "SB03SX.f"
    u_offset = 1 + u_dim1;
#line 212 "SB03SX.f"
    u -= u_offset;
#line 212 "SB03SX.f"
    r_dim1 = *ldr;
#line 212 "SB03SX.f"
    r_offset = 1 + r_dim1;
#line 212 "SB03SX.f"
    r__ -= r_offset;
#line 212 "SB03SX.f"
    --iwork;
#line 212 "SB03SX.f"
    --dwork;
#line 212 "SB03SX.f"

#line 212 "SB03SX.f"
    /* Function Body */
#line 212 "SB03SX.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 213 "SB03SX.f"
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

#line 215 "SB03SX.f"
    nn = *n * *n;
#line 216 "SB03SX.f"
    *info = 0;
#line 217 "SB03SX.f"
    if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || lsame_(trana,
	     "C", (ftnlen)1, (ftnlen)1))) {
#line 219 "SB03SX.f"
	*info = -1;
#line 220 "SB03SX.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 222 "SB03SX.f"
	*info = -2;
#line 223 "SB03SX.f"
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
#line 224 "SB03SX.f"
	*info = -3;
#line 225 "SB03SX.f"
    } else if (*n < 0) {
#line 226 "SB03SX.f"
	*info = -4;
#line 227 "SB03SX.f"
    } else if (*xanorm < 0.) {
#line 228 "SB03SX.f"
	*info = -5;
#line 229 "SB03SX.f"
    } else if (*ldt < max(1,*n)) {
#line 230 "SB03SX.f"
	*info = -7;
#line 231 "SB03SX.f"
    } else if (*ldu < 1 || update && *ldu < *n) {
#line 232 "SB03SX.f"
	*info = -9;
#line 233 "SB03SX.f"
    } else if (*ldr < max(1,*n)) {
#line 234 "SB03SX.f"
	*info = -11;
#line 235 "SB03SX.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 235 "SB03SX.f"
	i__1 = 3, i__2 = nn << 1;
#line 235 "SB03SX.f"
	if (*ldwork < 0 || *ldwork < max(i__1,i__2) && *n > 0) {
#line 237 "SB03SX.f"
	    *info = -15;
#line 238 "SB03SX.f"
	}
#line 238 "SB03SX.f"
    }

#line 240 "SB03SX.f"
    if (*info != 0) {
#line 241 "SB03SX.f"
	i__1 = -(*info);
#line 241 "SB03SX.f"
	xerbla_("SB03SX", &i__1, (ftnlen)6);
#line 242 "SB03SX.f"
	return 0;
#line 243 "SB03SX.f"
    }

/*     Quick return if possible. */

#line 247 "SB03SX.f"
    *ferr = 0.;
#line 248 "SB03SX.f"
    if (*n == 0 || *xanorm == 0.) {
#line 248 "SB03SX.f"
	return 0;
#line 248 "SB03SX.f"
    }

#line 251 "SB03SX.f"
    itmp = nn + 1;

#line 253 "SB03SX.f"
    if (notrna) {
#line 254 "SB03SX.f"
	*(unsigned char *)tranat = 'T';
#line 255 "SB03SX.f"
    } else {
#line 256 "SB03SX.f"
	*(unsigned char *)tranat = 'N';
#line 257 "SB03SX.f"
    }

/*     Fill in the remaining triangle of the symmetric residual matrix. */

#line 261 "SB03SX.f"
    ma02ed_(uplo, n, &r__[r_offset], ldr, (ftnlen)1);

#line 263 "SB03SX.f"
    kase = 0;

/*     REPEAT */
#line 266 "SB03SX.f"
L10:
#line 267 "SB03SX.f"
    dlacon_(&nn, &dwork[itmp], &dwork[1], &iwork[1], &est, &kase);
#line 268 "SB03SX.f"
    if (kase != 0) {

/*        Select the triangular part of symmetric matrix to be used. */

#line 272 "SB03SX.f"
	if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp], (ftnlen)
		6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &dwork[1], n, &
		dwork[itmp], (ftnlen)6, (ftnlen)5)) {
#line 276 "SB03SX.f"
	    *(unsigned char *)uplow = 'U';
#line 277 "SB03SX.f"
	    lower = FALSE_;
#line 278 "SB03SX.f"
	} else {
#line 279 "SB03SX.f"
	    *(unsigned char *)uplow = 'L';
#line 280 "SB03SX.f"
	    lower = TRUE_;
#line 281 "SB03SX.f"
	}

#line 283 "SB03SX.f"
	if (kase == 2) {
#line 284 "SB03SX.f"
	    ij = 0;
#line 285 "SB03SX.f"
	    if (lower) {

/*              Scale the lower triangular part of symmetric matrix */
/*              by the residual matrix. */

#line 290 "SB03SX.f"
		i__1 = *n;
#line 290 "SB03SX.f"
		for (j = 1; j <= i__1; ++j) {
#line 291 "SB03SX.f"
		    i__2 = *n;
#line 291 "SB03SX.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 292 "SB03SX.f"
			++ij;
#line 293 "SB03SX.f"
			dwork[ij] *= r__[i__ + j * r_dim1];
#line 294 "SB03SX.f"
/* L20: */
#line 294 "SB03SX.f"
		    }
#line 295 "SB03SX.f"
		    ij += j;
#line 296 "SB03SX.f"
/* L30: */
#line 296 "SB03SX.f"
		}
#line 297 "SB03SX.f"
	    } else {

/*              Scale the upper triangular part of symmetric matrix */
/*              by the residual matrix. */

#line 302 "SB03SX.f"
		i__1 = *n;
#line 302 "SB03SX.f"
		for (j = 1; j <= i__1; ++j) {
#line 303 "SB03SX.f"
		    i__2 = j;
#line 303 "SB03SX.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 304 "SB03SX.f"
			++ij;
#line 305 "SB03SX.f"
			dwork[ij] *= r__[i__ + j * r_dim1];
#line 306 "SB03SX.f"
/* L40: */
#line 306 "SB03SX.f"
		    }
#line 307 "SB03SX.f"
		    ij = ij + *n - j;
#line 308 "SB03SX.f"
/* L50: */
#line 308 "SB03SX.f"
		}
#line 309 "SB03SX.f"
	    }
#line 310 "SB03SX.f"
	}

#line 312 "SB03SX.f"
	if (update) {

/*           Transform the right-hand side: RHS := U'*RHS*U. */

#line 316 "SB03SX.f"
	    mb01ru_(uplow, "Transpose", n, n, &c_b24, &c_b25, &dwork[1], n, &
		    u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &info2,
		     (ftnlen)1, (ftnlen)9);
#line 318 "SB03SX.f"
	    i__1 = *n + 1;
#line 318 "SB03SX.f"
	    dscal_(n, &c_b26, &dwork[1], &i__1);
#line 319 "SB03SX.f"
	}
#line 320 "SB03SX.f"
	ma02ed_(uplow, n, &dwork[1], n, (ftnlen)1);

#line 322 "SB03SX.f"
	if (kase == 2) {

/*           Solve op(T)'*Y*op(T) - Y = scale*RHS. */

#line 326 "SB03SX.f"
	    sb03mx_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &dwork[
		    itmp], &info2, (ftnlen)1);
#line 328 "SB03SX.f"
	} else {

/*           Solve op(T)*W*op(T)' - W = scale*RHS. */

#line 332 "SB03SX.f"
	    sb03mx_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
		    dwork[itmp], &info2, (ftnlen)1);
#line 334 "SB03SX.f"
	}

#line 336 "SB03SX.f"
	if (info2 > 0) {
#line 336 "SB03SX.f"
	    *info = *n + 1;
#line 336 "SB03SX.f"
	}

#line 339 "SB03SX.f"
	if (update) {

/*           Transform back to obtain the solution: Z := U*Z*U', with */
/*           Z = Y or Z = W. */

#line 344 "SB03SX.f"
	    mb01ru_(uplow, "No transpose", n, n, &c_b24, &c_b25, &dwork[1], n,
		     &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &
		    info2, (ftnlen)1, (ftnlen)12);
#line 346 "SB03SX.f"
	    i__1 = *n + 1;
#line 346 "SB03SX.f"
	    dscal_(n, &c_b26, &dwork[1], &i__1);
#line 347 "SB03SX.f"
	}

#line 349 "SB03SX.f"
	if (kase == 1) {
#line 350 "SB03SX.f"
	    ij = 0;
#line 351 "SB03SX.f"
	    if (lower) {

/*              Scale the lower triangular part of symmetric matrix */
/*              by the residual matrix. */

#line 356 "SB03SX.f"
		i__1 = *n;
#line 356 "SB03SX.f"
		for (j = 1; j <= i__1; ++j) {
#line 357 "SB03SX.f"
		    i__2 = *n;
#line 357 "SB03SX.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 358 "SB03SX.f"
			++ij;
#line 359 "SB03SX.f"
			dwork[ij] *= r__[i__ + j * r_dim1];
#line 360 "SB03SX.f"
/* L60: */
#line 360 "SB03SX.f"
		    }
#line 361 "SB03SX.f"
		    ij += j;
#line 362 "SB03SX.f"
/* L70: */
#line 362 "SB03SX.f"
		}
#line 363 "SB03SX.f"
	    } else {

/*              Scale the upper triangular part of symmetric matrix */
/*              by the residual matrix. */

#line 368 "SB03SX.f"
		i__1 = *n;
#line 368 "SB03SX.f"
		for (j = 1; j <= i__1; ++j) {
#line 369 "SB03SX.f"
		    i__2 = j;
#line 369 "SB03SX.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 370 "SB03SX.f"
			++ij;
#line 371 "SB03SX.f"
			dwork[ij] *= r__[i__ + j * r_dim1];
#line 372 "SB03SX.f"
/* L80: */
#line 372 "SB03SX.f"
		    }
#line 373 "SB03SX.f"
		    ij = ij + *n - j;
#line 374 "SB03SX.f"
/* L90: */
#line 374 "SB03SX.f"
		}
#line 375 "SB03SX.f"
	    }
#line 376 "SB03SX.f"
	}

/*        Fill in the remaining triangle of the symmetric matrix. */

#line 380 "SB03SX.f"
	ma02ed_(uplow, n, &dwork[1], n, (ftnlen)1);
#line 381 "SB03SX.f"
	goto L10;
#line 382 "SB03SX.f"
    }

/*     UNTIL KASE = 0 */

/*     Compute the estimate of the relative error. */

#line 388 "SB03SX.f"
    temp = *xanorm * scale;
#line 389 "SB03SX.f"
    if (temp > est) {
#line 390 "SB03SX.f"
	*ferr = est / temp;
#line 391 "SB03SX.f"
    } else {
#line 392 "SB03SX.f"
	*ferr = 1.;
#line 393 "SB03SX.f"
    }

#line 395 "SB03SX.f"
    return 0;

/* *** Last line of SB03SX *** */
} /* sb03sx_ */

