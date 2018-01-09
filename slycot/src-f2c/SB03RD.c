#line 1 "SB03RD.f"
/* SB03RD.f -- translated by f2c (version 20100827).
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

#line 1 "SB03RD.f"
/* Table of constant values */

static doublereal c_b15 = 0.;
static doublereal c_b16 = 1.;
static integer c__1 = 1;

/* Subroutine */ int sb03rd_(char *job, char *fact, char *trana, integer *n, 
	doublereal *a, integer *lda, doublereal *u, integer *ldu, doublereal *
	c__, integer *ldc, doublereal *scale, doublereal *sep, doublereal *
	ferr, doublereal *wr, doublereal *wi, integer *iwork, doublereal *
	dwork, integer *ldwork, integer *info, ftnlen job_len, ftnlen 
	fact_len, ftnlen trana_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, c_dim1, c_offset, u_dim1, u_offset, i__1, i__2;

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
    extern /* Subroutine */ int sb03my_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), dcopy_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
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

/*     To solve the real Lyapunov matrix equation */

/*            op(A)'*X + X*op(A) = scale*C */

/*     and/or estimate the separation between the matrices op(A) and */
/*     -op(A)', where op(A) = A or A' (A**T) and C is symmetric (C = C'). */
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

/*     SEP     (output) DOUBLE PRECISION */
/*             If JOB = 'S' or JOB = 'B', and INFO = 0 or INFO = N+1, SEP */
/*             contains the estimated separation of the matrices op(A) */
/*             and -op(A)'. */
/*             If JOB = 'X' or N = 0, SEP is not referenced. */

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
/*                If FACT = 'F', LDWORK >= N*N; */
/*                If FACT = 'N', LDWORK >= MAX(N*N,3*N). */
/*             If JOB = 'S' or JOB = 'B' then */
/*                If FACT = 'F', LDWORK >= 2*N*N; */
/*                If FACT = 'N', LDWORK >= MAX(2*N*N,3*N). */
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
/*             = N+1:  if the matrices A and -A' have common or very */
/*                   close eigenvalues; perturbed values were used to */
/*                   solve the equation (but the matrix A is unchanged). */

/*     METHOD */

/*     After reducing matrix A to real Schur canonical form (if needed), */
/*     the Bartels-Stewart algorithm is used. A set of equivalent linear */
/*     algebraic systems of equations of order at most four are formed */
/*     and solved using Gaussian elimination with complete pivoting. */

/*     REFERENCES */

/*     [1] Bartels, R.H. and Stewart, G.W.  T */
/*         Solution of the matrix equation A X + XB = C. */
/*         Comm. A.C.M., 15, pp. 820-826, 1972. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations. */

/*     FURTHER COMMENTS */

/*     SEP is defined as the separation of op(A) and -op(A)': */

/*            sep( op(A), -op(A)' ) = sigma_min( T ) */

/*     where sigma_min(T) is the smallest singular value of the */
/*     N*N-by-N*N matrix */

/*        T = kprod( I(N), op(A)' ) + kprod( op(A), I(N) ). */

/*     I(N) is an N-by-N identity matrix, and kprod denotes the Kronecker */
/*     product. The program estimates sigma_min(T) by the reciprocal of */
/*     an estimate of the 1-norm of inverse(T). The true reciprocal */
/*     1-norm of inverse(T) cannot differ from sigma_min(T) by more */
/*     than a factor of N. */

/*     When SEP is small, small changes in A, C can cause large changes */
/*     in the solution of the equation. An approximate bound on the */
/*     maximum relative error in the computed solution is */

/*                            EPS * norm(A) / SEP */

/*     where EPS is the machine precision. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, May 1997. */
/*     Supersedes Release 2.0 routine MB03AD by Control Systems Research */
/*     Group, Kingston Polytechnic, United Kingdom, October 1982. */
/*     Based on DGELYP by P. Petkov, Tech. University of Sofia, September */
/*     1993. */

/*     REVISIONS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, May 1999. */

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

#line 259 "SB03RD.f"
    /* Parameter adjustments */
#line 259 "SB03RD.f"
    a_dim1 = *lda;
#line 259 "SB03RD.f"
    a_offset = 1 + a_dim1;
#line 259 "SB03RD.f"
    a -= a_offset;
#line 259 "SB03RD.f"
    u_dim1 = *ldu;
#line 259 "SB03RD.f"
    u_offset = 1 + u_dim1;
#line 259 "SB03RD.f"
    u -= u_offset;
#line 259 "SB03RD.f"
    c_dim1 = *ldc;
#line 259 "SB03RD.f"
    c_offset = 1 + c_dim1;
#line 259 "SB03RD.f"
    c__ -= c_offset;
#line 259 "SB03RD.f"
    --wr;
#line 259 "SB03RD.f"
    --wi;
#line 259 "SB03RD.f"
    --iwork;
#line 259 "SB03RD.f"
    --dwork;
#line 259 "SB03RD.f"

#line 259 "SB03RD.f"
    /* Function Body */
#line 259 "SB03RD.f"
    wantx = lsame_(job, "X", (ftnlen)1, (ftnlen)1);
#line 260 "SB03RD.f"
    wantsp = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 261 "SB03RD.f"
    wantbh = lsame_(job, "B", (ftnlen)1, (ftnlen)1);
#line 262 "SB03RD.f"
    nofact = lsame_(fact, "N", (ftnlen)1, (ftnlen)1);
#line 263 "SB03RD.f"
    nota = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);

#line 265 "SB03RD.f"
    *info = 0;
#line 266 "SB03RD.f"
    if (! wantsp && ! wantbh && ! wantx) {
#line 267 "SB03RD.f"
	*info = -1;
#line 268 "SB03RD.f"
    } else if (! nofact && ! lsame_(fact, "F", (ftnlen)1, (ftnlen)1)) {
#line 269 "SB03RD.f"
	*info = -2;
#line 270 "SB03RD.f"
    } else if (! nota && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! 
	    lsame_(trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 272 "SB03RD.f"
	*info = -3;
#line 273 "SB03RD.f"
    } else if (*n < 0) {
#line 274 "SB03RD.f"
	*info = -4;
#line 275 "SB03RD.f"
    } else if (*lda < max(1,*n)) {
#line 276 "SB03RD.f"
	*info = -6;
#line 277 "SB03RD.f"
    } else if (*ldu < max(1,*n)) {
#line 278 "SB03RD.f"
	*info = -8;
#line 279 "SB03RD.f"
    } else if (wantsp && *ldc < 1 || ! wantsp && *ldc < max(1,*n)) {
#line 281 "SB03RD.f"
	*info = -10;
#line 282 "SB03RD.f"
    }

/*     Compute workspace. */

#line 286 "SB03RD.f"
    if (wantx) {
#line 287 "SB03RD.f"
	if (nofact) {
/* Computing MAX */
#line 288 "SB03RD.f"
	    i__1 = *n * *n, i__2 = *n * 3;
#line 288 "SB03RD.f"
	    minwrk = max(i__1,i__2);
#line 289 "SB03RD.f"
	} else {
#line 290 "SB03RD.f"
	    minwrk = *n * *n;
#line 291 "SB03RD.f"
	}
#line 292 "SB03RD.f"
    } else {
#line 293 "SB03RD.f"
	if (nofact) {
/* Computing MAX */
#line 294 "SB03RD.f"
	    i__1 = (*n << 1) * *n, i__2 = *n * 3;
#line 294 "SB03RD.f"
	    minwrk = max(i__1,i__2);
#line 295 "SB03RD.f"
	} else {
#line 296 "SB03RD.f"
	    minwrk = (*n << 1) * *n;
#line 297 "SB03RD.f"
	}
#line 298 "SB03RD.f"
    }
#line 299 "SB03RD.f"
    if (*ldwork < max(1,minwrk)) {
#line 300 "SB03RD.f"
	*info = -18;
#line 301 "SB03RD.f"
    }

#line 303 "SB03RD.f"
    if (*info != 0) {
#line 304 "SB03RD.f"
	i__1 = -(*info);
#line 304 "SB03RD.f"
	xerbla_("SB03RD", &i__1, (ftnlen)6);
#line 305 "SB03RD.f"
	return 0;
#line 306 "SB03RD.f"
    }

/*     Quick return if possible. */

#line 310 "SB03RD.f"
    if (*n == 0) {
#line 311 "SB03RD.f"
	*scale = 1.;
#line 312 "SB03RD.f"
	if (wantbh) {
#line 312 "SB03RD.f"
	    *ferr = 0.;
#line 312 "SB03RD.f"
	}
#line 314 "SB03RD.f"
	dwork[1] = 1.;
#line 315 "SB03RD.f"
	return 0;
#line 316 "SB03RD.f"
    }

#line 318 "SB03RD.f"
    lwa = 0;

#line 320 "SB03RD.f"
    if (nofact) {

/*        Compute the Schur factorization of A. */
/*        Workspace:  need   3*N; */
/*                    prefer larger. */

#line 326 "SB03RD.f"
	dgees_("Vectors", "Not ordered", (L_fp)select_, n, &a[a_offset], lda, 
		&sdim, &wr[1], &wi[1], &u[u_offset], ldu, &dwork[1], ldwork, 
		bwork, info, (ftnlen)7, (ftnlen)11);
#line 328 "SB03RD.f"
	if (*info > 0) {
#line 328 "SB03RD.f"
	    return 0;
#line 328 "SB03RD.f"
	}
#line 330 "SB03RD.f"
	lwa = (integer) dwork[1];
#line 331 "SB03RD.f"
    }

#line 333 "SB03RD.f"
    if (! wantsp) {

/*        Transform the right-hand side. */
/*        Workspace:  need   N*N. */

#line 338 "SB03RD.f"
	*(unsigned char *)uplo = 'U';
#line 339 "SB03RD.f"
	mb01rd_(uplo, "Transpose", n, n, &c_b15, &c_b16, &c__[c_offset], ldc, 
		&u[u_offset], ldu, &c__[c_offset], ldc, &dwork[1], ldwork, 
		info, (ftnlen)1, (ftnlen)9);

#line 342 "SB03RD.f"
	i__1 = *n;
#line 342 "SB03RD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 343 "SB03RD.f"
	    i__2 = i__ - 1;
#line 343 "SB03RD.f"
	    dcopy_(&i__2, &c__[i__ * c_dim1 + 1], &c__1, &c__[i__ + c_dim1], 
		    ldc);
#line 344 "SB03RD.f"
/* L10: */
#line 344 "SB03RD.f"
	}

/*        Solve the transformed equation. */

#line 348 "SB03RD.f"
	sb03my_(trana, n, &a[a_offset], lda, &c__[c_offset], ldc, scale, info,
		 (ftnlen)1);
#line 349 "SB03RD.f"
	if (*info > 0) {
#line 349 "SB03RD.f"
	    *info = *n + 1;
#line 349 "SB03RD.f"
	}

/*        Transform back the solution. */

#line 354 "SB03RD.f"
	mb01rd_(uplo, "No transpose", n, n, &c_b15, &c_b16, &c__[c_offset], 
		ldc, &u[u_offset], ldu, &c__[c_offset], ldc, &dwork[1], 
		ldwork, info, (ftnlen)1, (ftnlen)12);

#line 357 "SB03RD.f"
	i__1 = *n;
#line 357 "SB03RD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 358 "SB03RD.f"
	    i__2 = i__ - 1;
#line 358 "SB03RD.f"
	    dcopy_(&i__2, &c__[i__ * c_dim1 + 1], &c__1, &c__[i__ + c_dim1], 
		    ldc);
#line 359 "SB03RD.f"
/* L20: */
#line 359 "SB03RD.f"
	}

#line 361 "SB03RD.f"
    }

#line 363 "SB03RD.f"
    if (! wantx) {

/*        Estimate sep(op(A),-op(A)'). */
/*        Workspace:  2*N*N. */

#line 368 "SB03RD.f"
	if (nota) {
#line 369 "SB03RD.f"
	    *(unsigned char *)notra = 'T';
#line 370 "SB03RD.f"
	} else {
#line 371 "SB03RD.f"
	    *(unsigned char *)notra = 'N';
#line 372 "SB03RD.f"
	}

#line 374 "SB03RD.f"
	est = 0.;
#line 375 "SB03RD.f"
	kase = 0;
/*        REPEAT */
#line 377 "SB03RD.f"
L30:
#line 378 "SB03RD.f"
	i__1 = *n * *n;
#line 378 "SB03RD.f"
	dlacon_(&i__1, &dwork[*n * *n + 1], &dwork[1], &iwork[1], &est, &kase)
		;
#line 379 "SB03RD.f"
	if (kase != 0) {
#line 380 "SB03RD.f"
	    if (kase == 1) {
#line 381 "SB03RD.f"
		sb03my_(trana, n, &a[a_offset], lda, &dwork[1], n, &scalef, &
			ierr, (ftnlen)1);
#line 382 "SB03RD.f"
	    } else {
#line 383 "SB03RD.f"
		sb03my_(notra, n, &a[a_offset], lda, &dwork[1], n, &scalef, &
			ierr, (ftnlen)1);
#line 384 "SB03RD.f"
	    }
#line 385 "SB03RD.f"
	    goto L30;
#line 386 "SB03RD.f"
	}
/*        UNTIL KASE = 0 */

#line 389 "SB03RD.f"
	*sep = scalef / est;

#line 391 "SB03RD.f"
	if (wantbh) {

/*           Compute the estimate of the relative error. */

#line 395 "SB03RD.f"
	    *ferr = dlamch_("Precision", (ftnlen)9) * dlanhs_("Frobenius", n, 
		    &a[a_offset], lda, &dwork[1], (ftnlen)9) / *sep;
#line 397 "SB03RD.f"
	}
#line 398 "SB03RD.f"
    }

#line 400 "SB03RD.f"
    dwork[1] = (doublereal) max(lwa,minwrk);

#line 402 "SB03RD.f"
    return 0;
/* *** Last line of SB03RD *** */
} /* sb03rd_ */

