#line 1 "SB03QX.f"
/* SB03QX.f -- translated by f2c (version 20100827).
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

#line 1 "SB03QX.f"
/* Table of constant values */

static doublereal c_b24 = 0.;
static doublereal c_b25 = 1.;
static doublereal c_b26 = .5;

/* Subroutine */ int sb03qx_(char *trana, char *uplo, char *lyapun, integer *
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
	    integer *, ftnlen, ftnlen), sb03my_(char *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen);
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
/*     continuous-time Lyapunov matrix equation, */

/*            op(A)'*X + X*op(A) = C, */

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
/*             The length of the array DWORK.  LDWORK >= 2*N*N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = N+1:  if the matrices T and -T' have common or very */
/*                   close eigenvalues; perturbed values were used to */
/*                   solve Lyapunov equations (but the matrix T is */
/*                   unchanged). */

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
/*     forward error bound for the solution of a continuous-time */
/*     algebraic matrix Riccati equation. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998. Partly based on DGLSVX (and then SB03QD) by P. Petkov, */
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

#line 211 "SB03QX.f"
    /* Parameter adjustments */
#line 211 "SB03QX.f"
    t_dim1 = *ldt;
#line 211 "SB03QX.f"
    t_offset = 1 + t_dim1;
#line 211 "SB03QX.f"
    t -= t_offset;
#line 211 "SB03QX.f"
    u_dim1 = *ldu;
#line 211 "SB03QX.f"
    u_offset = 1 + u_dim1;
#line 211 "SB03QX.f"
    u -= u_offset;
#line 211 "SB03QX.f"
    r_dim1 = *ldr;
#line 211 "SB03QX.f"
    r_offset = 1 + r_dim1;
#line 211 "SB03QX.f"
    r__ -= r_offset;
#line 211 "SB03QX.f"
    --iwork;
#line 211 "SB03QX.f"
    --dwork;
#line 211 "SB03QX.f"

#line 211 "SB03QX.f"
    /* Function Body */
#line 211 "SB03QX.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 212 "SB03QX.f"
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

#line 214 "SB03QX.f"
    nn = *n * *n;
#line 215 "SB03QX.f"
    *info = 0;
#line 216 "SB03QX.f"
    if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || lsame_(trana,
	     "C", (ftnlen)1, (ftnlen)1))) {
#line 218 "SB03QX.f"
	*info = -1;
#line 219 "SB03QX.f"
    } else if (! (lsame_(uplo, "L", (ftnlen)1, (ftnlen)1) || lsame_(uplo, 
	    "U", (ftnlen)1, (ftnlen)1))) {
#line 221 "SB03QX.f"
	*info = -2;
#line 222 "SB03QX.f"
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
#line 223 "SB03QX.f"
	*info = -3;
#line 224 "SB03QX.f"
    } else if (*n < 0) {
#line 225 "SB03QX.f"
	*info = -4;
#line 226 "SB03QX.f"
    } else if (*xanorm < 0.) {
#line 227 "SB03QX.f"
	*info = -5;
#line 228 "SB03QX.f"
    } else if (*ldt < max(1,*n)) {
#line 229 "SB03QX.f"
	*info = -7;
#line 230 "SB03QX.f"
    } else if (*ldu < 1 || update && *ldu < *n) {
#line 231 "SB03QX.f"
	*info = -9;
#line 232 "SB03QX.f"
    } else if (*ldr < max(1,*n)) {
#line 233 "SB03QX.f"
	*info = -11;
#line 234 "SB03QX.f"
    } else if (*ldwork < nn << 1) {
#line 235 "SB03QX.f"
	*info = -15;
#line 236 "SB03QX.f"
    }

#line 238 "SB03QX.f"
    if (*info != 0) {
#line 239 "SB03QX.f"
	i__1 = -(*info);
#line 239 "SB03QX.f"
	xerbla_("SB03QX", &i__1, (ftnlen)6);
#line 240 "SB03QX.f"
	return 0;
#line 241 "SB03QX.f"
    }

/*     Quick return if possible. */

#line 245 "SB03QX.f"
    *ferr = 0.;
#line 246 "SB03QX.f"
    if (*n == 0 || *xanorm == 0.) {
#line 246 "SB03QX.f"
	return 0;
#line 246 "SB03QX.f"
    }

#line 249 "SB03QX.f"
    itmp = nn + 1;

#line 251 "SB03QX.f"
    if (notrna) {
#line 252 "SB03QX.f"
	*(unsigned char *)tranat = 'T';
#line 253 "SB03QX.f"
    } else {
#line 254 "SB03QX.f"
	*(unsigned char *)tranat = 'N';
#line 255 "SB03QX.f"
    }

/*     Fill in the remaining triangle of the symmetric residual matrix. */

#line 259 "SB03QX.f"
    ma02ed_(uplo, n, &r__[r_offset], ldr, (ftnlen)1);

#line 261 "SB03QX.f"
    kase = 0;

/*     REPEAT */
#line 264 "SB03QX.f"
L10:
#line 265 "SB03QX.f"
    dlacon_(&nn, &dwork[itmp], &dwork[1], &iwork[1], &est, &kase);
#line 266 "SB03QX.f"
    if (kase != 0) {

/*        Select the triangular part of symmetric matrix to be used. */

#line 270 "SB03QX.f"
	if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp], (ftnlen)
		6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &dwork[1], n, &
		dwork[itmp], (ftnlen)6, (ftnlen)5)) {
#line 274 "SB03QX.f"
	    *(unsigned char *)uplow = 'U';
#line 275 "SB03QX.f"
	    lower = FALSE_;
#line 276 "SB03QX.f"
	} else {
#line 277 "SB03QX.f"
	    *(unsigned char *)uplow = 'L';
#line 278 "SB03QX.f"
	    lower = TRUE_;
#line 279 "SB03QX.f"
	}

#line 281 "SB03QX.f"
	if (kase == 2) {
#line 282 "SB03QX.f"
	    ij = 0;
#line 283 "SB03QX.f"
	    if (lower) {

/*              Scale the lower triangular part of symmetric matrix */
/*              by the residual matrix. */

#line 288 "SB03QX.f"
		i__1 = *n;
#line 288 "SB03QX.f"
		for (j = 1; j <= i__1; ++j) {
#line 289 "SB03QX.f"
		    i__2 = *n;
#line 289 "SB03QX.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 290 "SB03QX.f"
			++ij;
#line 291 "SB03QX.f"
			dwork[ij] *= r__[i__ + j * r_dim1];
#line 292 "SB03QX.f"
/* L20: */
#line 292 "SB03QX.f"
		    }
#line 293 "SB03QX.f"
		    ij += j;
#line 294 "SB03QX.f"
/* L30: */
#line 294 "SB03QX.f"
		}
#line 295 "SB03QX.f"
	    } else {

/*              Scale the upper triangular part of symmetric matrix */
/*              by the residual matrix. */

#line 300 "SB03QX.f"
		i__1 = *n;
#line 300 "SB03QX.f"
		for (j = 1; j <= i__1; ++j) {
#line 301 "SB03QX.f"
		    i__2 = j;
#line 301 "SB03QX.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 302 "SB03QX.f"
			++ij;
#line 303 "SB03QX.f"
			dwork[ij] *= r__[i__ + j * r_dim1];
#line 304 "SB03QX.f"
/* L40: */
#line 304 "SB03QX.f"
		    }
#line 305 "SB03QX.f"
		    ij = ij + *n - j;
#line 306 "SB03QX.f"
/* L50: */
#line 306 "SB03QX.f"
		}
#line 307 "SB03QX.f"
	    }
#line 308 "SB03QX.f"
	}

#line 310 "SB03QX.f"
	if (update) {

/*           Transform the right-hand side: RHS := U'*RHS*U. */

#line 314 "SB03QX.f"
	    mb01ru_(uplow, "Transpose", n, n, &c_b24, &c_b25, &dwork[1], n, &
		    u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &info2,
		     (ftnlen)1, (ftnlen)9);
#line 316 "SB03QX.f"
	    i__1 = *n + 1;
#line 316 "SB03QX.f"
	    dscal_(n, &c_b26, &dwork[1], &i__1);
#line 317 "SB03QX.f"
	}
#line 318 "SB03QX.f"
	ma02ed_(uplow, n, &dwork[1], n, (ftnlen)1);

#line 320 "SB03QX.f"
	if (kase == 2) {

/*           Solve op(T)'*Y + Y*op(T) = scale*RHS. */

#line 324 "SB03QX.f"
	    sb03my_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &info2,
		     (ftnlen)1);
#line 325 "SB03QX.f"
	} else {

/*           Solve op(T)*W + W*op(T)' = scale*RHS. */

#line 329 "SB03QX.f"
	    sb03my_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
		    info2, (ftnlen)1);
#line 330 "SB03QX.f"
	}

#line 332 "SB03QX.f"
	if (info2 > 0) {
#line 332 "SB03QX.f"
	    *info = *n + 1;
#line 332 "SB03QX.f"
	}

#line 335 "SB03QX.f"
	if (update) {

/*           Transform back to obtain the solution: Z := U*Z*U', with */
/*           Z = Y or Z = W. */

#line 340 "SB03QX.f"
	    mb01ru_(uplow, "No transpose", n, n, &c_b24, &c_b25, &dwork[1], n,
		     &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &
		    info2, (ftnlen)1, (ftnlen)12);
#line 342 "SB03QX.f"
	    i__1 = *n + 1;
#line 342 "SB03QX.f"
	    dscal_(n, &c_b26, &dwork[1], &i__1);
#line 343 "SB03QX.f"
	}

#line 345 "SB03QX.f"
	if (kase == 1) {
#line 346 "SB03QX.f"
	    ij = 0;
#line 347 "SB03QX.f"
	    if (lower) {

/*              Scale the lower triangular part of symmetric matrix */
/*              by the residual matrix. */

#line 352 "SB03QX.f"
		i__1 = *n;
#line 352 "SB03QX.f"
		for (j = 1; j <= i__1; ++j) {
#line 353 "SB03QX.f"
		    i__2 = *n;
#line 353 "SB03QX.f"
		    for (i__ = j; i__ <= i__2; ++i__) {
#line 354 "SB03QX.f"
			++ij;
#line 355 "SB03QX.f"
			dwork[ij] *= r__[i__ + j * r_dim1];
#line 356 "SB03QX.f"
/* L60: */
#line 356 "SB03QX.f"
		    }
#line 357 "SB03QX.f"
		    ij += j;
#line 358 "SB03QX.f"
/* L70: */
#line 358 "SB03QX.f"
		}
#line 359 "SB03QX.f"
	    } else {

/*              Scale the upper triangular part of symmetric matrix */
/*              by the residual matrix. */

#line 364 "SB03QX.f"
		i__1 = *n;
#line 364 "SB03QX.f"
		for (j = 1; j <= i__1; ++j) {
#line 365 "SB03QX.f"
		    i__2 = j;
#line 365 "SB03QX.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 366 "SB03QX.f"
			++ij;
#line 367 "SB03QX.f"
			dwork[ij] *= r__[i__ + j * r_dim1];
#line 368 "SB03QX.f"
/* L80: */
#line 368 "SB03QX.f"
		    }
#line 369 "SB03QX.f"
		    ij = ij + *n - j;
#line 370 "SB03QX.f"
/* L90: */
#line 370 "SB03QX.f"
		}
#line 371 "SB03QX.f"
	    }
#line 372 "SB03QX.f"
	}

/*        Fill in the remaining triangle of the symmetric matrix. */

#line 376 "SB03QX.f"
	ma02ed_(uplow, n, &dwork[1], n, (ftnlen)1);
#line 377 "SB03QX.f"
	goto L10;
#line 378 "SB03QX.f"
    }

/*     UNTIL KASE = 0 */

/*     Compute the estimate of the relative error. */

#line 384 "SB03QX.f"
    temp = *xanorm * scale;
#line 385 "SB03QX.f"
    if (temp > est) {
#line 386 "SB03QX.f"
	*ferr = est / temp;
#line 387 "SB03QX.f"
    } else {
#line 388 "SB03QX.f"
	*ferr = 1.;
#line 389 "SB03QX.f"
    }

#line 391 "SB03QX.f"
    return 0;

/* *** Last line of SB03QX *** */
} /* sb03qx_ */

