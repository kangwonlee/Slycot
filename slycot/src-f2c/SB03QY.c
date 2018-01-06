#line 1 "SB03QY.f"
/* SB03QY.f -- translated by f2c (version 20100827).
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

#line 1 "SB03QY.f"
/* Table of constant values */

static doublereal c_b21 = 0.;
static doublereal c_b22 = 1.;
static doublereal c_b23 = .5;

/* Subroutine */ int sb03qy_(char *job, char *trana, char *lyapun, integer *n,
	 doublereal *t, integer *ldt, doublereal *u, integer *ldu, doublereal 
	*x, integer *ldx, doublereal *sep, doublereal *thnorm, integer *iwork,
	 doublereal *dwork, integer *ldwork, integer *info, ftnlen job_len, 
	ftnlen trana_len, ftnlen lyapun_len)
{
    /* System generated locals */
    integer t_dim1, t_offset, u_dim1, u_offset, x_dim1, x_offset, i__1;

    /* Local variables */
    static integer nn;
    static doublereal est;
    static integer kase, itmp;
    static char uplo[1];
    static integer info2;
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
    static logical wants, wantt;
    extern /* Subroutine */ int dsyr2k_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacon_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
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

/*     To estimate the separation between the matrices op(A) and -op(A)', */

/*     sep(op(A),-op(A)') = min norm(op(A)'*X + X*op(A))/norm(X) */
/*                        = 1 / norm(inv(Omega)) */

/*     and/or the 1-norm of Theta, where op(A) = A or A' (A**T), and */
/*     Omega and Theta are linear operators associated to the real */
/*     continuous-time Lyapunov matrix equation */

/*            op(A)'*X + X*op(A) = C, */

/*     defined by */

/*     Omega(W) = op(A)'*W + W*op(A), */
/*     Theta(W) = inv(Omega(op(W)'*X + X*op(W))). */

/*     The 1-norm condition estimators are used. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'S':  Compute the separation only; */
/*             = 'T':  Compute the norm of Theta only; */
/*             = 'B':  Compute both the separation and the norm of Theta. */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

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
/*             The order of the matrices A and X.  N >= 0. */

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

/*     X       (input) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading N-by-N part of this array must contain the */
/*             solution matrix X of the Lyapunov equation (reduced */
/*             Lyapunov equation if LYAPUN = 'R'). */
/*             If JOB = 'S', the array X is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of array X. */
/*             LDX >= 1,        if JOB = 'S'; */
/*             LDX >= MAX(1,N), if JOB = 'T' or 'B'. */

/*     SEP     (output) DOUBLE PRECISION */
/*             If JOB = 'S' or JOB = 'B', and INFO >= 0, SEP contains the */
/*             estimated separation of the matrices op(A) and -op(A)'. */
/*             If JOB = 'T' or N = 0, SEP is not referenced. */

/*     THNORM  (output) DOUBLE PRECISION */
/*             If JOB = 'T' or JOB = 'B', and INFO >= 0, THNORM contains */
/*             the estimated 1-norm of operator Theta. */
/*             If JOB = 'S' or N = 0, THNORM is not referenced. */

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

/*     SEP is defined as the separation of op(A) and -op(A)': */

/*            sep( op(A), -op(A)' ) = sigma_min( K ) */

/*     where sigma_min(K) is the smallest singular value of the */
/*     N*N-by-N*N matrix */

/*        K = kprod( I(N), op(A)' ) + kprod( op(A)', I(N) ). */

/*     I(N) is an N-by-N identity matrix, and kprod denotes the Kronecker */
/*     product. The routine estimates sigma_min(K) by the reciprocal of */
/*     an estimate of the 1-norm of inverse(K), computed as suggested in */
/*     [1]. This involves the solution of several continuous-time */
/*     Lyapunov equations, either direct or transposed. The true */
/*     reciprocal 1-norm of inverse(K) cannot differ from sigma_min(K) by */
/*     more than a factor of N. */
/*     The 1-norm of Theta is estimated similarly. */

/*     REFERENCES */

/*     [1] Higham, N.J. */
/*         FORTRAN codes for estimating the one-norm of a real or */
/*         complex matrix, with applications to condition estimation. */
/*         ACM Trans. Math. Softw., 14, pp. 381-396, 1988. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     FURTHER COMMENTS */

/*     When SEP is zero, the routine returns immediately, with THNORM */
/*     (if requested) not set. In this case, the equation is singular. */
/*     The option LYAPUN = 'R' may occasionally produce slightly worse */
/*     or better estimates, and it is much faster than the option 'O'. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998. Partly based on DGLSVX (and then SB03QD) by P. Petkov, */
/*     Tech. University of Sofia, March 1998 (and December 1998). */

/*     REVISIONS */

/*     February 13, 1999, V. Sima, Katholieke Univ. Leuven, Belgium. */
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

#line 222 "SB03QY.f"
    /* Parameter adjustments */
#line 222 "SB03QY.f"
    t_dim1 = *ldt;
#line 222 "SB03QY.f"
    t_offset = 1 + t_dim1;
#line 222 "SB03QY.f"
    t -= t_offset;
#line 222 "SB03QY.f"
    u_dim1 = *ldu;
#line 222 "SB03QY.f"
    u_offset = 1 + u_dim1;
#line 222 "SB03QY.f"
    u -= u_offset;
#line 222 "SB03QY.f"
    x_dim1 = *ldx;
#line 222 "SB03QY.f"
    x_offset = 1 + x_dim1;
#line 222 "SB03QY.f"
    x -= x_offset;
#line 222 "SB03QY.f"
    --iwork;
#line 222 "SB03QY.f"
    --dwork;
#line 222 "SB03QY.f"

#line 222 "SB03QY.f"
    /* Function Body */
#line 222 "SB03QY.f"
    wants = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 223 "SB03QY.f"
    wantt = lsame_(job, "T", (ftnlen)1, (ftnlen)1);
#line 224 "SB03QY.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 225 "SB03QY.f"
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

#line 227 "SB03QY.f"
    nn = *n * *n;
#line 228 "SB03QY.f"
    *info = 0;
#line 229 "SB03QY.f"
    if (! (wants || wantt || lsame_(job, "B", (ftnlen)1, (ftnlen)1))) {
#line 230 "SB03QY.f"
	*info = -1;
#line 231 "SB03QY.f"
    } else if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1))) {
#line 233 "SB03QY.f"
	*info = -2;
#line 234 "SB03QY.f"
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
#line 235 "SB03QY.f"
	*info = -3;
#line 236 "SB03QY.f"
    } else if (*n < 0) {
#line 237 "SB03QY.f"
	*info = -4;
#line 238 "SB03QY.f"
    } else if (*ldt < max(1,*n)) {
#line 239 "SB03QY.f"
	*info = -6;
#line 240 "SB03QY.f"
    } else if (*ldu < 1 || update && *ldu < *n) {
#line 241 "SB03QY.f"
	*info = -8;
#line 242 "SB03QY.f"
    } else if (*ldx < 1 || ! wants && *ldx < *n) {
#line 243 "SB03QY.f"
	*info = -10;
#line 244 "SB03QY.f"
    } else if (*ldwork < nn << 1) {
#line 245 "SB03QY.f"
	*info = -15;
#line 246 "SB03QY.f"
    }

#line 248 "SB03QY.f"
    if (*info != 0) {
#line 249 "SB03QY.f"
	i__1 = -(*info);
#line 249 "SB03QY.f"
	xerbla_("SB03QY", &i__1, (ftnlen)6);
#line 250 "SB03QY.f"
	return 0;
#line 251 "SB03QY.f"
    }

/*     Quick return if possible. */

#line 255 "SB03QY.f"
    if (*n == 0) {
#line 255 "SB03QY.f"
	return 0;
#line 255 "SB03QY.f"
    }

#line 258 "SB03QY.f"
    itmp = nn + 1;

#line 260 "SB03QY.f"
    if (notrna) {
#line 261 "SB03QY.f"
	*(unsigned char *)tranat = 'T';
#line 262 "SB03QY.f"
    } else {
#line 263 "SB03QY.f"
	*(unsigned char *)tranat = 'N';
#line 264 "SB03QY.f"
    }

#line 266 "SB03QY.f"
    if (! wantt) {

/*        Estimate sep(op(A),-op(A)'). */
/*        Workspace:  2*N*N. */

#line 271 "SB03QY.f"
	kase = 0;

/*        REPEAT */
#line 274 "SB03QY.f"
L10:
#line 275 "SB03QY.f"
	dlacon_(&nn, &dwork[itmp], &dwork[1], &iwork[1], &est, &kase);
#line 276 "SB03QY.f"
	if (kase != 0) {

/*           Select the triangular part of symmetric matrix to be used. */

#line 280 "SB03QY.f"
	    if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp], (
		    ftnlen)6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &
		    dwork[1], n, &dwork[itmp], (ftnlen)6, (ftnlen)5)) {
#line 284 "SB03QY.f"
		*(unsigned char *)uplo = 'U';
#line 285 "SB03QY.f"
	    } else {
#line 286 "SB03QY.f"
		*(unsigned char *)uplo = 'L';
#line 287 "SB03QY.f"
	    }

#line 289 "SB03QY.f"
	    if (update) {

/*              Transform the right-hand side: RHS := U'*RHS*U. */

#line 293 "SB03QY.f"
		mb01ru_(uplo, "Transpose", n, n, &c_b21, &c_b22, &dwork[1], n,
			 &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &
			info2, (ftnlen)1, (ftnlen)9);
#line 296 "SB03QY.f"
		i__1 = *n + 1;
#line 296 "SB03QY.f"
		dscal_(n, &c_b23, &dwork[1], &i__1);
#line 297 "SB03QY.f"
	    }
#line 298 "SB03QY.f"
	    ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);

#line 300 "SB03QY.f"
	    if (kase == 1) {

/*              Solve op(T)'*Y + Y*op(T) = scale*RHS. */

#line 304 "SB03QY.f"
		sb03my_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			info2, (ftnlen)1);
#line 305 "SB03QY.f"
	    } else {

/*              Solve op(T)*W + W*op(T)' = scale*RHS. */

#line 309 "SB03QY.f"
		sb03my_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			info2, (ftnlen)1);
#line 310 "SB03QY.f"
	    }

#line 312 "SB03QY.f"
	    if (info2 > 0) {
#line 312 "SB03QY.f"
		*info = *n + 1;
#line 312 "SB03QY.f"
	    }

#line 315 "SB03QY.f"
	    if (update) {

/*              Transform back to obtain the solution: Z := U*Z*U', with */
/*              Z = Y or Z = W. */

#line 320 "SB03QY.f"
		mb01ru_(uplo, "No transpose", n, n, &c_b21, &c_b22, &dwork[1],
			 n, &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &
			nn, &info2, (ftnlen)1, (ftnlen)12);
#line 323 "SB03QY.f"
		i__1 = *n + 1;
#line 323 "SB03QY.f"
		dscal_(n, &c_b23, &dwork[1], &i__1);

/*              Fill in the remaining triangle of the symmetric matrix. */

#line 327 "SB03QY.f"
		ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);
#line 328 "SB03QY.f"
	    }

#line 330 "SB03QY.f"
	    goto L10;
#line 331 "SB03QY.f"
	}
/*        UNTIL KASE = 0 */

#line 334 "SB03QY.f"
	if (est > scale) {
#line 335 "SB03QY.f"
	    *sep = scale / est;
#line 336 "SB03QY.f"
	} else {
#line 337 "SB03QY.f"
	    bignum = 1. / dlamch_("Safe minimum", (ftnlen)12);
#line 338 "SB03QY.f"
	    if (scale < est * bignum) {
#line 339 "SB03QY.f"
		*sep = scale / est;
#line 340 "SB03QY.f"
	    } else {
#line 341 "SB03QY.f"
		*sep = bignum;
#line 342 "SB03QY.f"
	    }
#line 343 "SB03QY.f"
	}

/*        Return if the equation is singular. */

#line 347 "SB03QY.f"
	if (*sep == 0.) {
#line 347 "SB03QY.f"
	    return 0;
#line 347 "SB03QY.f"
	}
#line 349 "SB03QY.f"
    }

#line 351 "SB03QY.f"
    if (! wants) {

/*        Estimate norm(Theta). */
/*        Workspace:  2*N*N. */

#line 356 "SB03QY.f"
	kase = 0;

/*        REPEAT */
#line 359 "SB03QY.f"
L20:
#line 360 "SB03QY.f"
	dlacon_(&nn, &dwork[itmp], &dwork[1], &iwork[1], &est, &kase);
#line 361 "SB03QY.f"
	if (kase != 0) {

/*           Select the triangular part of symmetric matrix to be used. */

#line 365 "SB03QY.f"
	    if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp], (
		    ftnlen)6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &
		    dwork[1], n, &dwork[itmp], (ftnlen)6, (ftnlen)5)) {
#line 369 "SB03QY.f"
		*(unsigned char *)uplo = 'U';
#line 370 "SB03QY.f"
	    } else {
#line 371 "SB03QY.f"
		*(unsigned char *)uplo = 'L';
#line 372 "SB03QY.f"
	    }

/*           Fill in the remaining triangle of the symmetric matrix. */

#line 376 "SB03QY.f"
	    ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);

/*           Compute RHS = op(W)'*X + X*op(W). */

#line 380 "SB03QY.f"
	    dsyr2k_(uplo, tranat, n, n, &c_b22, &dwork[1], n, &x[x_offset], 
		    ldx, &c_b21, &dwork[itmp], n, (ftnlen)1, (ftnlen)1);
#line 382 "SB03QY.f"
	    dlacpy_(uplo, n, n, &dwork[itmp], n, &dwork[1], n, (ftnlen)1);

#line 384 "SB03QY.f"
	    if (update) {

/*              Transform the right-hand side: RHS := U'*RHS*U. */

#line 388 "SB03QY.f"
		mb01ru_(uplo, "Transpose", n, n, &c_b21, &c_b22, &dwork[1], n,
			 &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &
			info2, (ftnlen)1, (ftnlen)9);
#line 391 "SB03QY.f"
		i__1 = *n + 1;
#line 391 "SB03QY.f"
		dscal_(n, &c_b23, &dwork[1], &i__1);
#line 392 "SB03QY.f"
	    }
#line 393 "SB03QY.f"
	    ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);

#line 395 "SB03QY.f"
	    if (kase == 1) {

/*              Solve op(T)'*Y + Y*op(T) = scale*RHS. */

#line 399 "SB03QY.f"
		sb03my_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			info2, (ftnlen)1);
#line 400 "SB03QY.f"
	    } else {

/*              Solve op(T)*W + W*op(T)' = scale*RHS. */

#line 404 "SB03QY.f"
		sb03my_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			info2, (ftnlen)1);
#line 405 "SB03QY.f"
	    }

#line 407 "SB03QY.f"
	    if (info2 > 0) {
#line 407 "SB03QY.f"
		*info = *n + 1;
#line 407 "SB03QY.f"
	    }

#line 410 "SB03QY.f"
	    if (update) {

/*              Transform back to obtain the solution: Z := U*Z*U', with */
/*              Z = Y or Z = W. */

#line 415 "SB03QY.f"
		mb01ru_(uplo, "No transpose", n, n, &c_b21, &c_b22, &dwork[1],
			 n, &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &
			nn, &info2, (ftnlen)1, (ftnlen)12);
#line 418 "SB03QY.f"
		i__1 = *n + 1;
#line 418 "SB03QY.f"
		dscal_(n, &c_b23, &dwork[1], &i__1);

/*              Fill in the remaining triangle of the symmetric matrix. */

#line 422 "SB03QY.f"
		ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);
#line 423 "SB03QY.f"
	    }

#line 425 "SB03QY.f"
	    goto L20;
#line 426 "SB03QY.f"
	}
/*        UNTIL KASE = 0 */

#line 429 "SB03QY.f"
	if (est < scale) {
#line 430 "SB03QY.f"
	    *thnorm = est / scale;
#line 431 "SB03QY.f"
	} else {
#line 432 "SB03QY.f"
	    bignum = 1. / dlamch_("Safe minimum", (ftnlen)12);
#line 433 "SB03QY.f"
	    if (est < scale * bignum) {
#line 434 "SB03QY.f"
		*thnorm = est / scale;
#line 435 "SB03QY.f"
	    } else {
#line 436 "SB03QY.f"
		*thnorm = bignum;
#line 437 "SB03QY.f"
	    }
#line 438 "SB03QY.f"
	}
#line 439 "SB03QY.f"
    }

#line 441 "SB03QY.f"
    return 0;
/* *** Last line of SB03QY *** */
} /* sb03qy_ */

