#line 1 "SB03PD.f"
/* SB03PD.f -- translated by f2c (version 20100827).
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

#line 1 "SB03PD.f"
/* Table of constant values */

static doublereal c_b15 = 0.;
static doublereal c_b16 = 1.;
static integer c__1 = 1;

/* Subroutine */ int sb03pd_(char *job, char *fact, char *trana, integer *n, 
	doublereal *a, integer *lda, doublereal *u, integer *ldu, doublereal *
	c__, integer *ldc, doublereal *scale, doublereal *sepd, doublereal *
	ferr, doublereal *wr, doublereal *wi, integer *iwork, doublereal *
	dwork, integer *ldwork, integer *info, ftnlen job_len, ftnlen 
	fact_len, ftnlen trana_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, u_dim1, u_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, lwa;
    static doublereal est;
    static integer kase, sdim;
    static logical nota;
    static integer ierr;
    static char uplo[1];
    extern /* Subroutine */ int mb01rd_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen), dgees_(char *, char *, L_fp, integer *
	    , doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, logical *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int sb03mx_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static char notra[1];
    static logical bwork[1], wantx;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scalef;
    extern /* Subroutine */ int dlacon_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    static logical nofact;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    extern logical select_();
    static logical wantbh;
    static integer minwrk;
    static logical wantsp;


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

/*     To solve the real discrete Lyapunov matrix equation */

/*            op(A)'*X*op(A) - X = scale*C */

/*     and/or estimate the quantity, called separation, */

/*         sepd(op(A),op(A)') = min norm(op(A)'*X*op(A) - X)/norm(X) */

/*     where op(A) = A or A' (A**T) and C is symmetric (C = C'). */
/*     (A' denotes the transpose of the matrix A.) A is N-by-N, the right */
/*     hand side C and the solution X are N-by-N, and scale is an output */
/*     scale factor, set less than or equal to 1 to avoid overflow in X. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the computation to be performed, as follows: */
/*             = 'X':  Compute the solution only; */
/*             = 'S':  Compute the separation only; */
/*             = 'B':  Compute both the solution and the separation. */

/*     FACT    CHARACTER*1 */
/*             Specifies whether or not the real Schur factorization */
/*             of the matrix A is supplied on entry, as follows: */
/*             = 'F':  On entry, A and U contain the factors from the */
/*                     real Schur factorization of the matrix A; */
/*             = 'N':  The Schur factorization of A will be computed */
/*                     and the factors will be stored in A and U. */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, X, and C.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. If FACT = 'F', then A contains */
/*             an upper quasi-triangular matrix in Schur canonical form. */
/*             On exit, if INFO = 0 or INFO = N+1, the leading N-by-N */
/*             part of this array contains the upper quasi-triangular */
/*             matrix in Schur canonical form from the Shur factorization */
/*             of A. The contents of array A is not modified if */
/*             FACT = 'F'. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     U       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDU,N) */
/*             If FACT = 'F', then U is an input argument and on entry */
/*             it must contain the orthogonal matrix U from the real */
/*             Schur factorization of A. */
/*             If FACT = 'N', then U is an output argument and on exit, */
/*             if INFO = 0 or INFO = N+1, it contains the orthogonal */
/*             N-by-N matrix from the real Schur factorization of A. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry with JOB = 'X' or 'B', the leading N-by-N part of */
/*             this array must contain the symmetric matrix C. */
/*             On exit with JOB = 'X' or 'B', if INFO = 0 or INFO = N+1, */
/*             the leading N-by-N part of C has been overwritten by the */
/*             symmetric solution matrix X. */
/*             If JOB = 'S', C is not referenced. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= 1,        if JOB = 'S'; */
/*             LDC >= MAX(1,N), otherwise. */

/*     SCALE   (output) DOUBLE PRECISION */
/*             The scale factor, scale, set less than or equal to 1 to */
/*             prevent the solution overflowing. */

/*     SEPD    (output) DOUBLE PRECISION */
/*             If JOB = 'S' or JOB = 'B', and INFO = 0 or INFO = N+1, */
/*             SEPD contains the estimate in the 1-norm of */
/*             sepd(op(A),op(A)'). */
/*             If JOB = 'X' or N = 0, SEPD is not referenced. */

/*     FERR    (output) DOUBLE PRECISION */
/*             If JOB = 'B', and INFO = 0 or INFO = N+1, FERR contains */
/*             an estimated forward error bound for the solution X. */
/*             If XTRUE is the true solution, FERR bounds the relative */
/*             error in the computed solution, measured in the Frobenius */
/*             norm:  norm(X - XTRUE)/norm(XTRUE). */
/*             If JOB = 'X' or JOB = 'S', FERR is not referenced. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             If FACT = 'N', and INFO = 0 or INFO = N+1, WR and WI */
/*             contain the real and imaginary parts, respectively, of the */
/*             eigenvalues of A. */
/*             If FACT = 'F', WR and WI are not referenced. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N*N) */
/*             This array is not referenced if JOB = 'X'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 or INFO = N+1, DWORK(1) returns the */
/*             optimal value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 1 and */
/*             If JOB = 'X' then */
/*                If FACT = 'F', LDWORK >= MAX(N*N,2*N); */
/*                If FACT = 'N', LDWORK >= MAX(N*N,3*N). */
/*             If JOB = 'S' or JOB = 'B' then */
/*                LDWORK >= 2*N*N + 2*N. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, the QR algorithm failed to compute all */
/*                   the eigenvalues (see LAPACK Library routine DGEES); */
/*                   elements i+1:n of WR and WI contain eigenvalues */
/*                   which have converged, and A contains the partially */
/*                   converged Schur form; */
/*             = N+1:  if matrix A has almost reciprocal eigenvalues; */
/*                   perturbed values were used to solve the equation */
/*                   (but the matrix A is unchanged). */

/*     METHOD */

/*     After reducing matrix A to real Schur canonical form (if needed), */
/*     a discrete-time version of the Bartels-Stewart algorithm is used. */
/*     A set of equivalent linear algebraic systems of equations of order */
/*     at most four are formed and solved using Gaussian elimination with */
/*     complete pivoting. */

/*     REFERENCES */

/*     [1] Barraud, A.Y.                   T */
/*         A numerical algorithm to solve A XA - X = Q. */
/*         IEEE Trans. Auto. Contr., AC-22, pp. 883-885, 1977. */

/*     [2] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     FURTHER COMMENTS */

/*     SEPD is defined as */

/*            sepd( op(A), op(A)' ) = sigma_min( T ) */

/*     where sigma_min(T) is the smallest singular value of the */
/*     N*N-by-N*N matrix */

/*        T = kprod( op(A)', op(A)' ) - I(N**2). */

/*     I(N**2) is an N*N-by-N*N identity matrix, and kprod denotes the */
/*     Kronecker product. The program estimates sigma_min(T) by the */
/*     reciprocal of an estimate of the 1-norm of inverse(T). The true */
/*     reciprocal 1-norm of inverse(T) cannot differ from sigma_min(T) by */
/*     more than a factor of N. */

/*     When SEPD is small, small changes in A, C can cause large changes */
/*     in the solution of the equation. An approximate bound on the */
/*     maximum relative error in the computed solution is */

/*                            EPS * norm(A)**2 / SEPD */

/*     where EPS is the machine precision. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Supersedes Release 2.0 routine MB03AD by Control Systems Research */
/*     Group, Kingston Polytechnic, United Kingdom, October 1982. */
/*     Based on DGELPD by P. Petkov, Tech. University of Sofia, September */
/*     1993. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999. */

/*     KEYWORDS */

/*     Lyapunov equation, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

#line 267 "SB03PD.f"
    /* Parameter adjustments */
#line 267 "SB03PD.f"
    a_dim1 = *lda;
#line 267 "SB03PD.f"
    a_offset = 1 + a_dim1;
#line 267 "SB03PD.f"
    a -= a_offset;
#line 267 "SB03PD.f"
    u_dim1 = *ldu;
#line 267 "SB03PD.f"
    u_offset = 1 + u_dim1;
#line 267 "SB03PD.f"
    u -= u_offset;
#line 267 "SB03PD.f"
    c_dim1 = *ldc;
#line 267 "SB03PD.f"
    c_offset = 1 + c_dim1;
#line 267 "SB03PD.f"
    c__ -= c_offset;
#line 267 "SB03PD.f"
    --wr;
#line 267 "SB03PD.f"
    --wi;
#line 267 "SB03PD.f"
    --iwork;
#line 267 "SB03PD.f"
    --dwork;
#line 267 "SB03PD.f"

#line 267 "SB03PD.f"
    /* Function Body */
#line 267 "SB03PD.f"
    wantx = lsame_(job, "X", (ftnlen)1, (ftnlen)1);
#line 268 "SB03PD.f"
    wantsp = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 269 "SB03PD.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 270 "SB03PD.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 271 "SB03PD.f"
    nota = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);

#line 273 "SB03PD.f"
    *info = 0;
#line 274 "SB03PD.f"
    if (! wantbh && ! wantsp && ! wantx) {
#line 275 "SB03PD.f"
	*info = -1;
#line 276 "SB03PD.f"
    } else if (! nofact && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 277 "SB03PD.f"
	*info = -2;
#line 278 "SB03PD.f"
    } else if (! nota && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 280 "SB03PD.f"
	*info = -3;
#line 281 "SB03PD.f"
    } else if (*n < 0) {
#line 282 "SB03PD.f"
	*info = -4;
#line 283 "SB03PD.f"
    } else if (*lda < max(1,*n)) {
#line 284 "SB03PD.f"
	*info = -6;
#line 285 "SB03PD.f"
    } else if (*ldu < max(1,*n)) {
#line 286 "SB03PD.f"
	*info = -8;
#line 287 "SB03PD.f"
    } else if (wantsp && *ldc < 1 || ! wantsp && *ldc < max(1,*n)) {
#line 289 "SB03PD.f"
	*info = -10;
#line 290 "SB03PD.f"
    }

/*     Compute workspace. */

#line 294 "SB03PD.f"
    if (wantx) {
#line 295 "SB03PD.f"
	if (nofact) {
/* Computing MAX */
#line 296 "SB03PD.f"
	    i__1 = *n * *n, i__2 = *n * 3;
#line 296 "SB03PD.f"
	    minwrk = max(i__1,i__2);
#line 297 "SB03PD.f"
	} else {
/* Computing MAX */
#line 298 "SB03PD.f"
	    i__1 = *n * *n, i__2 = *n << 1;
#line 298 "SB03PD.f"
	    minwrk = max(i__1,i__2);
#line 299 "SB03PD.f"
	}
#line 300 "SB03PD.f"
    } else {
#line 301 "SB03PD.f"
	minwrk = (*n << 1) * *n + (*n << 1);
#line 302 "SB03PD.f"
    }
#line 303 "SB03PD.f"
    if (*ldwork < max(1,minwrk)) {
#line 304 "SB03PD.f"
	*info = -18;
#line 305 "SB03PD.f"
    }
#line 306 "SB03PD.f"
    if (*info != 0) {
#line 307 "SB03PD.f"
	i__1 = -(*info);
#line 307 "SB03PD.f"
	xerbla_("SB03PD", &i__1, (ftnlen)6);
#line 308 "SB03PD.f"
	return 0;
#line 309 "SB03PD.f"
    }

/*     Quick return if possible. */

#line 313 "SB03PD.f"
    if (*n == 0) {
#line 314 "SB03PD.f"
	*scale = 1.;
#line 315 "SB03PD.f"
	if (wantbh) {
#line 315 "SB03PD.f"
	    *ferr = 0.;
#line 315 "SB03PD.f"
	}
#line 317 "SB03PD.f"
	dwork[1] = 1.;
#line 318 "SB03PD.f"
	return 0;
#line 319 "SB03PD.f"
    }

#line 321 "SB03PD.f"
    lwa = 0;

#line 323 "SB03PD.f"
    if (nofact) {

/*        Compute the Schur factorization of A. */
/*        Workspace:  need   3*N; */
/*                    prefer larger. */

#line 329 "SB03PD.f"
	dgees_("Vectors", "Not ordered", (L_fp)select_, n, &a[a_offset], lda, 
		&sdim, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], ldwork, 
		bwork, info, (ftnlen)7, (ftnlen)11);
#line 331 "SB03PD.f"
	if (*info > 0) {
#line 331 "SB03PD.f"
	    return 0;
#line 331 "SB03PD.f"
	}
#line 333 "SB03PD.f"
	lwa = (integer) dwork[1];
#line 334 "SB03PD.f"
    }

#line 336 "SB03PD.f"
    if (! wantsp) {

/*        Transform the right-hand side. */
/*        Workspace:  need   N*N. */

#line 341 "SB03PD.f"
	*(unsigned char *)uplo = 'U';
#line 342 "SB03PD.f"
	mb01rd_(uplo, "Transpose", n, n, &c_b15, &c_b16, &c__[c_offset], ldc, 
		&u[u_offset], ldu, &c__[c_offset], ldc, &dwork[1], ldwork, 
		info, (ftnlen)1, (ftnlen)9);

#line 345 "SB03PD.f"
	i__1 = *n;
#line 345 "SB03PD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 346 "SB03PD.f"
	    i__2 = i__ - 1;
#line 346 "SB03PD.f"
	    dcopy_(&i__2, &c__[i__ * c_dim1 + 1], &c__1, &c__[i__ + c_dim1], 
		    ldc);
#line 347 "SB03PD.f"
/* L10: */
#line 347 "SB03PD.f"
	}

/*        Solve the transformed equation. */
/*        Workspace:  2*N. */

#line 352 "SB03PD.f"
	sb03mx_(trana, n, &a[a_offset], lda, &c__[c_offset], ldc, scale, &
		dwork[1], info, (ftnlen)1);
#line 353 "SB03PD.f"
	if (*info > 0) {
#line 353 "SB03PD.f"
	    *info = *n + 1;
#line 353 "SB03PD.f"
	}

/*        Transform back the solution. */

#line 358 "SB03PD.f"
	mb01rd_(uplo, "No transpose", n, n, &c_b15, &c_b16, &c__[c_offset], 
		ldc, &u[u_offset], ldu, &c__[c_offset], ldc, &dwork[1], 
		ldwork, info, (ftnlen)1, (ftnlen)12);

#line 361 "SB03PD.f"
	i__1 = *n;
#line 361 "SB03PD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 362 "SB03PD.f"
	    i__2 = i__ - 1;
#line 362 "SB03PD.f"
	    dcopy_(&i__2, &c__[i__ * c_dim1 + 1], &c__1, &c__[i__ + c_dim1], 
		    ldc);
#line 363 "SB03PD.f"
/* L20: */
#line 363 "SB03PD.f"
	}

#line 365 "SB03PD.f"
    }

#line 367 "SB03PD.f"
    if (! wantx) {

/*        Estimate sepd(op(A),op(A)'). */
/*        Workspace:  2*N*N + 2*N. */

#line 372 "SB03PD.f"
	if (nota) {
#line 373 "SB03PD.f"
	    *(unsigned char *)notra = 'T';
#line 374 "SB03PD.f"
	} else {
#line 375 "SB03PD.f"
	    *(unsigned char *)notra = 'N';
#line 376 "SB03PD.f"
	}

#line 378 "SB03PD.f"
	est = 0.;
#line 379 "SB03PD.f"
	kase = 0;
/*        REPEAT */
#line 381 "SB03PD.f"
L30:
#line 382 "SB03PD.f"
	i__1 = *n * *n;
#line 382 "SB03PD.f"
	dlacon_(&i__1, &dwork[*n * *n + 1], &dwork[1], &iwork[1], &est, &kase)
		;
#line 383 "SB03PD.f"
	if (kase != 0) {
#line 384 "SB03PD.f"
	    if (kase == 1) {
#line 385 "SB03PD.f"
		sb03mx_(trana, n, &a[a_offset], lda, &dwork[1], n, &scalef, &
			dwork[(*n << 1) * *n + 1], &ierr, (ftnlen)1);
#line 387 "SB03PD.f"
	    } else {
#line 388 "SB03PD.f"
		sb03mx_(notra, n, &a[a_offset], lda, &dwork[1], n, &scalef, &
			dwork[(*n << 1) * *n + 1], &ierr, (ftnlen)1);
#line 390 "SB03PD.f"
	    }
#line 391 "SB03PD.f"
	    goto L30;
#line 392 "SB03PD.f"
	}
/*        UNTIL KASE = 0 */

#line 395 "SB03PD.f"
	*sepd = scalef / est;

#line 397 "SB03PD.f"
	if (wantbh) {

/*           Compute the estimate of the relative error. */

/* Computing 2nd power */
#line 401 "SB03PD.f"
	    d__1 = dlanhs_("Frobenius", n, &a[a_offset], lda, &dwork[1], (
		    ftnlen)9);
#line 401 "SB03PD.f"
	    *ferr = dlamch_("Precision", (ftnlen)9) * (d__1 * d__1) / *sepd;
#line 403 "SB03PD.f"
	}
#line 404 "SB03PD.f"
    }

#line 406 "SB03PD.f"
    dwork[1] = (doublereal) max(lwa,minwrk);

#line 408 "SB03PD.f"
    return 0;
/* *** Last line of SB03PD *** */
} /* sb03pd_ */

