#line 1 "SB03SY.f"
/* SB03SY.f -- translated by f2c (version 20100827).
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

#line 1 "SB03SY.f"
/* Table of constant values */

static doublereal c_b21 = 0.;
static doublereal c_b22 = 1.;
static doublereal c_b23 = .5;

/* Subroutine */ int sb03sy_(char *job, char *trana, char *lyapun, integer *n,
	 doublereal *t, integer *ldt, doublereal *u, integer *ldu, doublereal 
	*xa, integer *ldxa, doublereal *sepd, doublereal *thnorm, integer *
	iwork, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	job_len, ftnlen trana_len, ftnlen lyapun_len)
{
    /* System generated locals */
    integer t_dim1, t_offset, u_dim1, u_offset, xa_dim1, xa_offset, i__1, 
	    i__2;

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
	    integer *, ftnlen, ftnlen), sb03mx_(char *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen);
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

/*     To estimate the "separation" between the matrices op(A) and */
/*     op(A)', */

/*     sepd(op(A),op(A)') = min norm(op(A)'*X*op(A) - X)/norm(X) */
/*                        = 1 / norm(inv(Omega)) */

/*     and/or the 1-norm of Theta, where op(A) = A or A' (A**T), and */
/*     Omega and Theta are linear operators associated to the real */
/*     discrete-time Lyapunov matrix equation */

/*            op(A)'*X*op(A) - X = C, */

/*     defined by */

/*     Omega(W) = op(A)'*W*op(A) - W, */
/*     Theta(W) = inv(Omega(op(W)'*X*op(A) + op(A)'*X*op(W))). */

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

/*     XA      (input) DOUBLE PRECISION array, dimension (LDXA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix product X*op(A), if LYAPUN = 'O', or U'*X*U*op(T), */
/*             if LYAPUN = 'R', in the Lyapunov equation. */
/*             If JOB = 'S', the array XA is not referenced. */

/*     LDXA    INTEGER */
/*             The leading dimension of array XA. */
/*             LDXA >= 1,        if JOB = 'S'; */
/*             LDXA >= MAX(1,N), if JOB = 'T' or 'B'. */

/*     SEPD    (output) DOUBLE PRECISION */
/*             If JOB = 'S' or JOB = 'B', and INFO >= 0, SEPD contains */
/*             the estimated quantity sepd(op(A),op(A)'). */
/*             If JOB = 'T' or N = 0, SEPD is not referenced. */

/*     THNORM  (output) DOUBLE PRECISION */
/*             If JOB = 'T' or JOB = 'B', and INFO >= 0, THNORM contains */
/*             the estimated 1-norm of operator Theta. */
/*             If JOB = 'S' or N = 0, THNORM is not referenced. */

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
/*             = N+1:  if T has (almost) reciprocal eigenvalues; */
/*                   perturbed values were used to solve Lyapunov */
/*                   equations (but the matrix T is unchanged). */

/*     METHOD */

/*     SEPD is defined as */

/*            sepd( op(A), op(A)' ) = sigma_min( K ) */

/*     where sigma_min(K) is the smallest singular value of the */
/*     N*N-by-N*N matrix */

/*        K = kprod( op(A)', op(A)' ) - I(N**2). */

/*     I(N**2) is an N*N-by-N*N identity matrix, and kprod denotes the */
/*     Kronecker product. The routine estimates sigma_min(K) by the */
/*     reciprocal of an estimate of the 1-norm of inverse(K), computed as */
/*     suggested in [1]. This involves the solution of several discrete- */
/*     time Lyapunov equations, either direct or transposed. The true */
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

/*     When SEPD is zero, the routine returns immediately, with THNORM */
/*     (if requested) not set. In this case, the equation is singular. */
/*     The option LYAPUN = 'R' may occasionally produce slightly worse */
/*     or better estimates, and it is much faster than the option 'O'. */

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

#line 225 "SB03SY.f"
    /* Parameter adjustments */
#line 225 "SB03SY.f"
    t_dim1 = *ldt;
#line 225 "SB03SY.f"
    t_offset = 1 + t_dim1;
#line 225 "SB03SY.f"
    t -= t_offset;
#line 225 "SB03SY.f"
    u_dim1 = *ldu;
#line 225 "SB03SY.f"
    u_offset = 1 + u_dim1;
#line 225 "SB03SY.f"
    u -= u_offset;
#line 225 "SB03SY.f"
    xa_dim1 = *ldxa;
#line 225 "SB03SY.f"
    xa_offset = 1 + xa_dim1;
#line 225 "SB03SY.f"
    xa -= xa_offset;
#line 225 "SB03SY.f"
    --iwork;
#line 225 "SB03SY.f"
    --dwork;
#line 225 "SB03SY.f"

#line 225 "SB03SY.f"
    /* Function Body */
#line 225 "SB03SY.f"
    wants = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 226 "SB03SY.f"
    wantt = lsame_(job, "T", (ftnlen)1, (ftnlen)1);
#line 227 "SB03SY.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 228 "SB03SY.f"
    update = lsame_(lyapun, "O", (ftnlen)1, (ftnlen)1);

#line 230 "SB03SY.f"
    nn = *n * *n;
#line 231 "SB03SY.f"
    *info = 0;
#line 232 "SB03SY.f"
    if (! (wants || wantt || lsame_(job, "B", (ftnlen)1, (ftnlen)1))) {
#line 233 "SB03SY.f"
	*info = -1;
#line 234 "SB03SY.f"
    } else if (! (notrna || lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1))) {
#line 236 "SB03SY.f"
	*info = -2;
#line 237 "SB03SY.f"
    } else if (! (update || lsame_(lyapun, "R", (ftnlen)1, (ftnlen)1))) {
#line 238 "SB03SY.f"
	*info = -3;
#line 239 "SB03SY.f"
    } else if (*n < 0) {
#line 240 "SB03SY.f"
	*info = -4;
#line 241 "SB03SY.f"
    } else if (*ldt < max(1,*n)) {
#line 242 "SB03SY.f"
	*info = -6;
#line 243 "SB03SY.f"
    } else if (*ldu < 1 || update && *ldu < *n) {
#line 244 "SB03SY.f"
	*info = -8;
#line 245 "SB03SY.f"
    } else if (*ldxa < 1 || ! wants && *ldxa < *n) {
#line 246 "SB03SY.f"
	*info = -10;
#line 247 "SB03SY.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 247 "SB03SY.f"
	i__1 = 3, i__2 = nn << 1;
#line 247 "SB03SY.f"
	if (*ldwork < 0 || *ldwork < max(i__1,i__2) && *n > 0) {
#line 249 "SB03SY.f"
	    *info = -15;
#line 250 "SB03SY.f"
	}
#line 250 "SB03SY.f"
    }

#line 252 "SB03SY.f"
    if (*info != 0) {
#line 253 "SB03SY.f"
	i__1 = -(*info);
#line 253 "SB03SY.f"
	xerbla_("SB03SY", &i__1, (ftnlen)6);
#line 254 "SB03SY.f"
	return 0;
#line 255 "SB03SY.f"
    }

/*     Quick return if possible. */

#line 259 "SB03SY.f"
    if (*n == 0) {
#line 259 "SB03SY.f"
	return 0;
#line 259 "SB03SY.f"
    }

#line 262 "SB03SY.f"
    itmp = nn + 1;

#line 264 "SB03SY.f"
    if (notrna) {
#line 265 "SB03SY.f"
	*(unsigned char *)tranat = 'T';
#line 266 "SB03SY.f"
    } else {
#line 267 "SB03SY.f"
	*(unsigned char *)tranat = 'N';
#line 268 "SB03SY.f"
    }

#line 270 "SB03SY.f"
    if (! wantt) {

/*        Estimate sepd(op(A),op(A)'). */
/*        Workspace:  max(3,2*N*N). */

#line 275 "SB03SY.f"
	kase = 0;

/*        REPEAT */
#line 278 "SB03SY.f"
L10:
#line 279 "SB03SY.f"
	dlacon_(&nn, &dwork[itmp], &dwork[1], &iwork[1], &est, &kase);
#line 280 "SB03SY.f"
	if (kase != 0) {

/*           Select the triangular part of symmetric matrix to be used. */

#line 284 "SB03SY.f"
	    if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp], (
		    ftnlen)6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &
		    dwork[1], n, &dwork[itmp], (ftnlen)6, (ftnlen)5)) {
#line 288 "SB03SY.f"
		*(unsigned char *)uplo = 'U';
#line 289 "SB03SY.f"
	    } else {
#line 290 "SB03SY.f"
		*(unsigned char *)uplo = 'L';
#line 291 "SB03SY.f"
	    }

#line 293 "SB03SY.f"
	    if (update) {

/*              Transform the right-hand side: RHS := U'*RHS*U. */

#line 297 "SB03SY.f"
		mb01ru_(uplo, "Transpose", n, n, &c_b21, &c_b22, &dwork[1], n,
			 &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &
			info2, (ftnlen)1, (ftnlen)9);
#line 300 "SB03SY.f"
		i__1 = *n + 1;
#line 300 "SB03SY.f"
		dscal_(n, &c_b23, &dwork[1], &i__1);
#line 301 "SB03SY.f"
	    }
#line 302 "SB03SY.f"
	    ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);

#line 304 "SB03SY.f"
	    if (kase == 1) {

/*              Solve op(T)'*Y*op(T) - Y = scale*RHS. */

#line 308 "SB03SY.f"
		sb03mx_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[itmp], &info2, (ftnlen)1);
#line 310 "SB03SY.f"
	    } else {

/*              Solve op(T)*W*op(T)' - W = scale*RHS. */

#line 314 "SB03SY.f"
		sb03mx_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[itmp], &info2, (ftnlen)1);
#line 316 "SB03SY.f"
	    }

#line 318 "SB03SY.f"
	    if (info2 > 0) {
#line 318 "SB03SY.f"
		*info = *n + 1;
#line 318 "SB03SY.f"
	    }

#line 321 "SB03SY.f"
	    if (update) {

/*              Transform back to obtain the solution: Z := U*Z*U', with */
/*              Z = Y or Z = W. */

#line 326 "SB03SY.f"
		mb01ru_(uplo, "No transpose", n, n, &c_b21, &c_b22, &dwork[1],
			 n, &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &
			nn, &info2, (ftnlen)1, (ftnlen)12);
#line 329 "SB03SY.f"
		i__1 = *n + 1;
#line 329 "SB03SY.f"
		dscal_(n, &c_b23, &dwork[1], &i__1);

/*              Fill in the remaining triangle of the symmetric matrix. */

#line 333 "SB03SY.f"
		ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);
#line 334 "SB03SY.f"
	    }

#line 336 "SB03SY.f"
	    goto L10;
#line 337 "SB03SY.f"
	}
/*        UNTIL KASE = 0 */

#line 340 "SB03SY.f"
	if (est > scale) {
#line 341 "SB03SY.f"
	    *sepd = scale / est;
#line 342 "SB03SY.f"
	} else {
#line 343 "SB03SY.f"
	    bignum = 1. / dlamch_("Safe minimum", (ftnlen)12);
#line 344 "SB03SY.f"
	    if (scale < est * bignum) {
#line 345 "SB03SY.f"
		*sepd = scale / est;
#line 346 "SB03SY.f"
	    } else {
#line 347 "SB03SY.f"
		*sepd = bignum;
#line 348 "SB03SY.f"
	    }
#line 349 "SB03SY.f"
	}

/*        Return if the equation is singular. */

#line 353 "SB03SY.f"
	if (*sepd == 0.) {
#line 353 "SB03SY.f"
	    return 0;
#line 353 "SB03SY.f"
	}
#line 355 "SB03SY.f"
    }

#line 357 "SB03SY.f"
    if (! wants) {

/*        Estimate norm(Theta). */
/*        Workspace:  max(3,2*N*N). */

#line 362 "SB03SY.f"
	kase = 0;

/*        REPEAT */
#line 365 "SB03SY.f"
L20:
#line 366 "SB03SY.f"
	dlacon_(&nn, &dwork[itmp], &dwork[1], &iwork[1], &est, &kase);
#line 367 "SB03SY.f"
	if (kase != 0) {

/*           Select the triangular part of symmetric matrix to be used. */

#line 371 "SB03SY.f"
	    if (dlansy_("1-norm", "Upper", n, &dwork[1], n, &dwork[itmp], (
		    ftnlen)6, (ftnlen)5) >= dlansy_("1-norm", "Lower", n, &
		    dwork[1], n, &dwork[itmp], (ftnlen)6, (ftnlen)5)) {
#line 375 "SB03SY.f"
		*(unsigned char *)uplo = 'U';
#line 376 "SB03SY.f"
	    } else {
#line 377 "SB03SY.f"
		*(unsigned char *)uplo = 'L';
#line 378 "SB03SY.f"
	    }

/*           Fill in the remaining triangle of the symmetric matrix. */

#line 382 "SB03SY.f"
	    ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);

/*           Compute RHS = op(W)'*X*op(A) + op(A)'*X*op(W). */

#line 386 "SB03SY.f"
	    dsyr2k_(uplo, tranat, n, n, &c_b22, &dwork[1], n, &xa[xa_offset], 
		    ldxa, &c_b21, &dwork[itmp], n, (ftnlen)1, (ftnlen)1);
#line 388 "SB03SY.f"
	    dlacpy_(uplo, n, n, &dwork[itmp], n, &dwork[1], n, (ftnlen)1);

#line 390 "SB03SY.f"
	    if (update) {

/*              Transform the right-hand side: RHS := U'*RHS*U. */

#line 394 "SB03SY.f"
		mb01ru_(uplo, "Transpose", n, n, &c_b21, &c_b22, &dwork[1], n,
			 &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &nn, &
			info2, (ftnlen)1, (ftnlen)9);
#line 397 "SB03SY.f"
		i__1 = *n + 1;
#line 397 "SB03SY.f"
		dscal_(n, &c_b23, &dwork[1], &i__1);
#line 398 "SB03SY.f"
	    }
#line 399 "SB03SY.f"
	    ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);

#line 401 "SB03SY.f"
	    if (kase == 1) {

/*              Solve op(T)'*Y*op(T) - Y = scale*RHS. */

#line 405 "SB03SY.f"
		sb03mx_(trana, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[itmp], &info2, (ftnlen)1);
#line 407 "SB03SY.f"
	    } else {

/*              Solve op(T)*W*op(T)' - W = scale*RHS. */

#line 411 "SB03SY.f"
		sb03mx_(tranat, n, &t[t_offset], ldt, &dwork[1], n, &scale, &
			dwork[itmp], &info2, (ftnlen)1);
#line 413 "SB03SY.f"
	    }

#line 415 "SB03SY.f"
	    if (info2 > 0) {
#line 415 "SB03SY.f"
		*info = *n + 1;
#line 415 "SB03SY.f"
	    }

#line 418 "SB03SY.f"
	    if (update) {

/*              Transform back to obtain the solution: Z := U*Z*U', with */
/*              Z = Y or Z = W. */

#line 423 "SB03SY.f"
		mb01ru_(uplo, "No transpose", n, n, &c_b21, &c_b22, &dwork[1],
			 n, &u[u_offset], ldu, &dwork[1], n, &dwork[itmp], &
			nn, &info2, (ftnlen)1, (ftnlen)12);
#line 426 "SB03SY.f"
		i__1 = *n + 1;
#line 426 "SB03SY.f"
		dscal_(n, &c_b23, &dwork[1], &i__1);

/*              Fill in the remaining triangle of the symmetric matrix. */

#line 430 "SB03SY.f"
		ma02ed_(uplo, n, &dwork[1], n, (ftnlen)1);
#line 431 "SB03SY.f"
	    }

#line 433 "SB03SY.f"
	    goto L20;
#line 434 "SB03SY.f"
	}
/*        UNTIL KASE = 0 */

#line 437 "SB03SY.f"
	if (est < scale) {
#line 438 "SB03SY.f"
	    *thnorm = est / scale;
#line 439 "SB03SY.f"
	} else {
#line 440 "SB03SY.f"
	    bignum = 1. / dlamch_("Safe minimum", (ftnlen)12);
#line 441 "SB03SY.f"
	    if (est < scale * bignum) {
#line 442 "SB03SY.f"
		*thnorm = est / scale;
#line 443 "SB03SY.f"
	    } else {
#line 444 "SB03SY.f"
		*thnorm = bignum;
#line 445 "SB03SY.f"
	    }
#line 446 "SB03SY.f"
	}
#line 447 "SB03SY.f"
    }

#line 449 "SB03SY.f"
    return 0;
/* *** Last line of SB03SY *** */
} /* sb03sy_ */

