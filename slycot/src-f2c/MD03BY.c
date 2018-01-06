#line 1 "MD03BY.f"
/* MD03BY.f -- translated by f2c (version 20100827).
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

#line 1 "MD03BY.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int md03by_(char *cond, integer *n, doublereal *r__, integer 
	*ldr, integer *ipvt, doublereal *diag, doublereal *qtb, doublereal *
	delta, doublereal *par, integer *rank, doublereal *x, doublereal *rx, 
	doublereal *tol, doublereal *dwork, integer *ldwork, integer *info, 
	ftnlen cond_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer j, l, n2;
    static doublereal fp, dum[3], parc;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal parl;
    static logical sing;
    static integer iter;
    static doublereal temp, paru;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int mb03od_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen);
    static logical econd;
    extern /* Subroutine */ int mb02yd_(char *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);
    static char condl[1];
    static logical ncond;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dwarf, dmino;
    static logical ucond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static doublereal gnorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dtrmv_(char *, char *, char *
	    , integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen, ftnlen), dtrsv_(char *, char *, char *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen,
	     ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal dxnorm;


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

/*     To determine a value for the parameter PAR such that if x solves */
/*     the system */

/*           A*x = b ,     sqrt(PAR)*D*x = 0 , */

/*     in the least squares sense, where A is an m-by-n matrix, D is an */
/*     n-by-n nonsingular diagonal matrix, and b is an m-vector, and if */
/*     DELTA is a positive number, DXNORM is the Euclidean norm of D*x, */
/*     then either PAR is zero and */

/*           ( DXNORM - DELTA ) .LE. 0.1*DELTA , */

/*     or PAR is positive and */

/*           ABS( DXNORM - DELTA ) .LE. 0.1*DELTA . */

/*     It is assumed that a QR factorization, with column pivoting, of A */
/*     is available, that is, A*P = Q*R, where P is a permutation matrix, */
/*     Q has orthogonal columns, and R is an upper triangular matrix */
/*     with diagonal elements of nonincreasing magnitude. */
/*     The routine needs the full upper triangle of R, the permutation */
/*     matrix P, and the first n components of Q'*b (' denotes the */
/*     transpose). On output, MD03BY also provides an upper triangular */
/*     matrix S such that */

/*           P'*(A'*A + PAR*D*D)*P = S'*S . */

/*     Matrix S is used in the solution process. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of the matrices R and S */
/*             should be estimated, as follows: */
/*             = 'E' :  use incremental condition estimation for R and S; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of R and S for zero values; */
/*             = 'U' :  use the rank already stored in RANK (for R). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix R.  N >= 0. */

/*     R       (input/output) DOUBLE PRECISION array, dimension (LDR, N) */
/*             On entry, the leading N-by-N upper triangular part of this */
/*             array must contain the upper triangular matrix R. */
/*             On exit, the full upper triangle is unaltered, and the */
/*             strict lower triangle contains the strict upper triangle */
/*             (transposed) of the upper triangular matrix S. */

/*     LDR     INTEGER */
/*             The leading dimension of array R.  LDR >= MAX(1,N). */

/*     IPVT    (input) INTEGER array, dimension (N) */
/*             This array must define the permutation matrix P such that */
/*             A*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     DIAG    (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the diagonal elements of the */
/*             matrix D.  DIAG(I) <> 0, I = 1,...,N. */

/*     QTB     (input) DOUBLE PRECISION array, dimension (N) */
/*             This array must contain the first n elements of the */
/*             vector Q'*b. */

/*     DELTA   (input) DOUBLE PRECISION */
/*             An upper bound on the Euclidean norm of D*x.  DELTA > 0. */

/*     PAR     (input/output) DOUBLE PRECISION */
/*             On entry, PAR must contain an initial estimate of the */
/*             Levenberg-Marquardt parameter.  PAR >= 0. */
/*             On exit, it contains the final estimate of this parameter. */

/*     RANK    (input or output) INTEGER */
/*             On entry, if COND = 'U', this parameter must contain the */
/*             (numerical) rank of the matrix R. */
/*             On exit, this parameter contains the numerical rank of */
/*             the matrix S. */

/*     X       (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the least squares solution of the */
/*             system A*x = b, sqrt(PAR)*D*x = 0. */

/*     RX      (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the matrix-vector product -R*P'*x. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             rank of the matrices R and S. If the user sets TOL > 0, */
/*             then the given value of TOL is used as a lower bound for */
/*             the reciprocal condition number;  a (sub)matrix whose */
/*             estimated condition number is less than 1/TOL is */
/*             considered to be of full rank.  If the user sets TOL <= 0, */
/*             then an implicitly computed, default tolerance, defined by */
/*             TOLDEF = N*EPS,  is used instead, where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH). */
/*             This parameter is not relevant if COND = 'U' or 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, the first N elements of this array contain the */
/*             diagonal elements of the upper triangular matrix S. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 4*N, if COND =  'E'; */
/*             LDWORK >= 2*N, if COND <> 'E'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The algorithm computes the Gauss-Newton direction. A least squares */
/*     solution is found if the Jacobian is rank deficient. If the Gauss- */
/*     Newton direction is not acceptable, then an iterative algorithm */
/*     obtains improved lower and upper bounds for the parameter PAR. */
/*     Only a few iterations are generally needed for convergence of the */
/*     algorithm. If, however, the limit of ITMAX = 10 iterations is */
/*     reached, then the output PAR will contain the best value obtained */
/*     so far. If the Gauss-Newton step is acceptable, it is stored in x, */
/*     and PAR is set to zero, hence S = R. */

/*     REFERENCES */

/*     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E. */
/*         User's Guide for MINPACK-1. */
/*         Applied Math. Division, Argonne National Laboratory, Argonne, */
/*         Illinois, Report ANL-80-74, 1980. */

/*     NUMERICAL ASPECTS */
/*                               2 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     FURTHER COMMENTS */

/*     This routine is a LAPACK-based modification of LMPAR from the */
/*     MINPACK package [1], and with optional condition estimation. */
/*     The option COND = 'U' is useful when dealing with several */
/*     right-hand side vectors, but RANK should be reset. */
/*     If COND = 'E', but the matrix S is guaranteed to be nonsingular */
/*     and well conditioned relative to TOL, i.e., rank(R) = N, and */
/*     min(DIAG) > 0, then its condition is not estimated. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2005. */

/*     KEYWORDS */

/*     Linear system of equations, matrix operations, plane rotations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 228 "MD03BY.f"
    /* Parameter adjustments */
#line 228 "MD03BY.f"
    r_dim1 = *ldr;
#line 228 "MD03BY.f"
    r_offset = 1 + r_dim1;
#line 228 "MD03BY.f"
    r__ -= r_offset;
#line 228 "MD03BY.f"
    --ipvt;
#line 228 "MD03BY.f"
    --diag;
#line 228 "MD03BY.f"
    --qtb;
#line 228 "MD03BY.f"
    --x;
#line 228 "MD03BY.f"
    --rx;
#line 228 "MD03BY.f"
    --dwork;
#line 228 "MD03BY.f"

#line 228 "MD03BY.f"
    /* Function Body */
#line 228 "MD03BY.f"
    econd = lsame_(cond, "E", (ftnlen)1, (ftnlen)1);
#line 229 "MD03BY.f"
    ncond = lsame_(cond, "N", (ftnlen)1, (ftnlen)1);
#line 230 "MD03BY.f"
    ucond = lsame_(cond, "U", (ftnlen)1, (ftnlen)1);
#line 231 "MD03BY.f"
    *info = 0;
#line 232 "MD03BY.f"
    if (! (econd || ncond || ucond)) {
#line 233 "MD03BY.f"
	*info = -1;
#line 234 "MD03BY.f"
    } else if (*n < 0) {
#line 235 "MD03BY.f"
	*info = -2;
#line 236 "MD03BY.f"
    } else if (*ldr < max(1,*n)) {
#line 237 "MD03BY.f"
	*info = -4;
#line 238 "MD03BY.f"
    } else if (*delta <= 0.) {
#line 239 "MD03BY.f"
	*info = -8;
#line 240 "MD03BY.f"
    } else if (*par < 0.) {
#line 241 "MD03BY.f"
	*info = -9;
#line 242 "MD03BY.f"
    } else if (ucond && (*rank < 0 || *rank > *n)) {
#line 243 "MD03BY.f"
	*info = -10;
#line 244 "MD03BY.f"
    } else if (*ldwork < *n << 1 || econd && *ldwork < *n << 2) {
#line 245 "MD03BY.f"
	*info = -15;
#line 246 "MD03BY.f"
    } else if (*n > 0) {
#line 247 "MD03BY.f"
	dmino = diag[1];
#line 248 "MD03BY.f"
	sing = FALSE_;

#line 250 "MD03BY.f"
	i__1 = *n;
#line 250 "MD03BY.f"
	for (j = 1; j <= i__1; ++j) {
#line 251 "MD03BY.f"
	    if (diag[j] < dmino) {
#line 251 "MD03BY.f"
		dmino = diag[j];
#line 251 "MD03BY.f"
	    }
#line 253 "MD03BY.f"
	    sing = sing || diag[j] == 0.;
#line 254 "MD03BY.f"
/* L10: */
#line 254 "MD03BY.f"
	}

#line 256 "MD03BY.f"
	if (sing) {
#line 256 "MD03BY.f"
	    *info = -6;
#line 256 "MD03BY.f"
	}
#line 258 "MD03BY.f"
    }

/*     Return if there are illegal arguments. */

#line 262 "MD03BY.f"
    if (*info != 0) {
#line 263 "MD03BY.f"
	i__1 = -(*info);
#line 263 "MD03BY.f"
	xerbla_("MD03BY", &i__1, (ftnlen)6);
#line 264 "MD03BY.f"
	return 0;
#line 265 "MD03BY.f"
    }

/*     Quick return if possible. */

#line 269 "MD03BY.f"
    if (*n == 0) {
#line 270 "MD03BY.f"
	*par = 0.;
#line 271 "MD03BY.f"
	*rank = 0;
#line 272 "MD03BY.f"
	return 0;
#line 273 "MD03BY.f"
    }

/*     DWARF is the smallest positive magnitude. */

#line 277 "MD03BY.f"
    dwarf = dlamch_("Underflow", (ftnlen)9);
#line 278 "MD03BY.f"
    n2 = *n;

/*     Estimate the rank of R, if required. */

#line 282 "MD03BY.f"
    if (econd) {
#line 283 "MD03BY.f"
	n2 = *n << 1;
#line 284 "MD03BY.f"
	temp = *tol;
#line 285 "MD03BY.f"
	if (temp <= 0.) {

/*           Use the default tolerance in rank determination. */

#line 289 "MD03BY.f"
	    temp = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
#line 290 "MD03BY.f"
	}

/*        Estimate the reciprocal condition number of R and set the rank. */
/*        Workspace: 2*N. */

#line 295 "MD03BY.f"
	mb03od_("No QR", n, n, &r__[r_offset], ldr, &ipvt[1], &temp, &c_b10, &
		dwork[1], rank, dum, &dwork[1], ldwork, info, (ftnlen)5);

#line 298 "MD03BY.f"
    } else if (ncond) {
#line 299 "MD03BY.f"
	j = 1;

#line 301 "MD03BY.f"
L20:
#line 302 "MD03BY.f"
	if (r__[j + j * r_dim1] != 0.) {
#line 303 "MD03BY.f"
	    ++j;
#line 304 "MD03BY.f"
	    if (j <= *n) {
#line 304 "MD03BY.f"
		goto L20;
#line 304 "MD03BY.f"
	    }
#line 306 "MD03BY.f"
	}

#line 308 "MD03BY.f"
	*rank = j - 1;
#line 309 "MD03BY.f"
    }

/*     Compute and store in x the Gauss-Newton direction. If the */
/*     Jacobian is rank-deficient, obtain a least squares solution. */
/*     The array RX is used as workspace. */

#line 315 "MD03BY.f"
    dcopy_(rank, &qtb[1], &c__1, &rx[1], &c__1);
#line 316 "MD03BY.f"
    dum[0] = 0.;
#line 317 "MD03BY.f"
    if (*rank < *n) {
#line 317 "MD03BY.f"
	i__1 = *n - *rank;
#line 317 "MD03BY.f"
	dcopy_(&i__1, dum, &c__0, &rx[*rank + 1], &c__1);
#line 317 "MD03BY.f"
    }
#line 319 "MD03BY.f"
    dtrsv_("Upper", "No transpose", "Non unit", rank, &r__[r_offset], ldr, &
	    rx[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);

#line 322 "MD03BY.f"
    i__1 = *n;
#line 322 "MD03BY.f"
    for (j = 1; j <= i__1; ++j) {
#line 323 "MD03BY.f"
	l = ipvt[j];
#line 324 "MD03BY.f"
	x[l] = rx[j];
#line 325 "MD03BY.f"
/* L30: */
#line 325 "MD03BY.f"
    }

/*     Initialize the iteration counter. */
/*     Evaluate the function at the origin, and test */
/*     for acceptance of the Gauss-Newton direction. */

#line 331 "MD03BY.f"
    iter = 0;

#line 333 "MD03BY.f"
    i__1 = *n;
#line 333 "MD03BY.f"
    for (j = 1; j <= i__1; ++j) {
#line 334 "MD03BY.f"
	dwork[j] = diag[j] * x[j];
#line 335 "MD03BY.f"
/* L40: */
#line 335 "MD03BY.f"
    }

#line 337 "MD03BY.f"
    dxnorm = dnrm2_(n, &dwork[1], &c__1);
#line 338 "MD03BY.f"
    fp = dxnorm - *delta;
#line 339 "MD03BY.f"
    if (fp > *delta * .1) {

/*        Set an appropriate option for estimating the condition of */
/*        the matrix S. */

#line 344 "MD03BY.f"
	if (ucond) {
#line 345 "MD03BY.f"
	    if (*ldwork >= *n << 2) {
#line 346 "MD03BY.f"
		*(unsigned char *)condl = 'E';
#line 347 "MD03BY.f"
		toldef = (doublereal) (*n) * dlamch_("Epsilon", (ftnlen)7);
#line 348 "MD03BY.f"
	    } else {
#line 349 "MD03BY.f"
		*(unsigned char *)condl = 'N';
#line 350 "MD03BY.f"
		toldef = *tol;
#line 351 "MD03BY.f"
	    }
#line 352 "MD03BY.f"
	} else {
#line 353 "MD03BY.f"
	    *(unsigned char *)condl = *(unsigned char *)cond;
#line 354 "MD03BY.f"
	    toldef = *tol;
#line 355 "MD03BY.f"
	}

/*        If the Jacobian is not rank deficient, the Newton */
/*        step provides a lower bound, PARL, for the zero of */
/*        the function. Otherwise set this bound to zero. */

#line 361 "MD03BY.f"
	if (*rank == *n) {

#line 363 "MD03BY.f"
	    i__1 = *n;
#line 363 "MD03BY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 364 "MD03BY.f"
		l = ipvt[j];
#line 365 "MD03BY.f"
		rx[j] = diag[l] * (dwork[l] / dxnorm);
#line 366 "MD03BY.f"
/* L50: */
#line 366 "MD03BY.f"
	    }

#line 368 "MD03BY.f"
	    dtrsv_("Upper", "Transpose", "Non unit", n, &r__[r_offset], ldr, &
		    rx[1], &c__1, (ftnlen)5, (ftnlen)9, (ftnlen)8);
#line 370 "MD03BY.f"
	    temp = dnrm2_(n, &rx[1], &c__1);
#line 371 "MD03BY.f"
	    parl = fp / *delta / temp / temp;

/*           For efficiency, use CONDL = 'U', if possible. */

#line 375 "MD03BY.f"
	    if (! lsame_(condl, "U", (ftnlen)1, (ftnlen)1) && dmino > 0.) {
#line 375 "MD03BY.f"
		*(unsigned char *)condl = 'U';
#line 375 "MD03BY.f"
	    }
#line 377 "MD03BY.f"
	} else {
#line 378 "MD03BY.f"
	    parl = 0.;
#line 379 "MD03BY.f"
	}

/*        Calculate an upper bound, PARU, for the zero of the function. */

#line 383 "MD03BY.f"
	i__1 = *n;
#line 383 "MD03BY.f"
	for (j = 1; j <= i__1; ++j) {
#line 384 "MD03BY.f"
	    l = ipvt[j];
#line 385 "MD03BY.f"
	    rx[j] = ddot_(&j, &r__[j * r_dim1 + 1], &c__1, &qtb[1], &c__1) / 
		    diag[l];
#line 386 "MD03BY.f"
/* L60: */
#line 386 "MD03BY.f"
	}

#line 388 "MD03BY.f"
	gnorm = dnrm2_(n, &rx[1], &c__1);
#line 389 "MD03BY.f"
	paru = gnorm / *delta;
#line 390 "MD03BY.f"
	if (paru == 0.) {
#line 390 "MD03BY.f"
	    paru = dwarf / min(*delta,.1) / .001;
#line 390 "MD03BY.f"
	}

/*        If the input PAR lies outside of the interval (PARL,PARU), */
/*        set PAR to the closer endpoint. */

#line 396 "MD03BY.f"
	*par = max(*par,parl);
#line 397 "MD03BY.f"
	*par = min(*par,paru);
#line 398 "MD03BY.f"
	if (*par == 0.) {
#line 398 "MD03BY.f"
	    *par = gnorm / dxnorm;
#line 398 "MD03BY.f"
	}

/*        Beginning of an iteration. */

#line 403 "MD03BY.f"
L70:
#line 404 "MD03BY.f"
	++iter;

/*           Evaluate the function at the current value of PAR. */

#line 408 "MD03BY.f"
	if (*par == 0.) {
/* Computing MAX */
#line 408 "MD03BY.f"
	    d__1 = dwarf, d__2 = paru * .001;
#line 408 "MD03BY.f"
	    *par = max(d__1,d__2);
#line 408 "MD03BY.f"
	}
#line 410 "MD03BY.f"
	temp = sqrt(*par);

#line 412 "MD03BY.f"
	i__1 = *n;
#line 412 "MD03BY.f"
	for (j = 1; j <= i__1; ++j) {
#line 413 "MD03BY.f"
	    rx[j] = temp * diag[j];
#line 414 "MD03BY.f"
/* L80: */
#line 414 "MD03BY.f"
	}

/*           Solve the system A*x = b , sqrt(PAR)*D*x = 0 , in a least */
/*           square sense. The first N elements of DWORK contain the */
/*           diagonal elements of the upper triangular matrix S, and */
/*           the next N elements contain the vector z, so that x = P*z. */
/*           The vector z is preserved if COND = 'E'. */
/*           Workspace:   4*N, if CONDL =  'E'; */
/*                        2*N, if CONDL <> 'E'. */

#line 424 "MD03BY.f"
	mb02yd_(condl, n, &r__[r_offset], ldr, &ipvt[1], &rx[1], &qtb[1], 
		rank, &x[1], &toldef, &dwork[1], ldwork, info, (ftnlen)1);

#line 427 "MD03BY.f"
	i__1 = *n;
#line 427 "MD03BY.f"
	for (j = 1; j <= i__1; ++j) {
#line 428 "MD03BY.f"
	    dwork[n2 + j] = diag[j] * x[j];
#line 429 "MD03BY.f"
/* L90: */
#line 429 "MD03BY.f"
	}

#line 431 "MD03BY.f"
	dxnorm = dnrm2_(n, &dwork[n2 + 1], &c__1);
#line 432 "MD03BY.f"
	temp = fp;
#line 433 "MD03BY.f"
	fp = dxnorm - *delta;

/*           If the function is small enough, accept the current value */
/*           of PAR. Also test for the exceptional cases where PARL */
/*           is zero or the number of iterations has reached ITMAX. */

#line 439 "MD03BY.f"
	if (abs(fp) > *delta * .1 && (parl != 0. || fp > temp || temp >= 0.) 
		&& iter < 10) {

/*              Compute the Newton correction. */

#line 445 "MD03BY.f"
	    i__1 = *rank;
#line 445 "MD03BY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 446 "MD03BY.f"
		l = ipvt[j];
#line 447 "MD03BY.f"
		rx[j] = diag[l] * (dwork[n2 + l] / dxnorm);
#line 448 "MD03BY.f"
/* L100: */
#line 448 "MD03BY.f"
	    }

#line 450 "MD03BY.f"
	    if (*rank < *n) {
#line 450 "MD03BY.f"
		i__1 = *n - *rank;
#line 450 "MD03BY.f"
		dcopy_(&i__1, dum, &c__0, &rx[*rank + 1], &c__1);
#line 450 "MD03BY.f"
	    }
#line 452 "MD03BY.f"
	    i__1 = *ldr + 1;
#line 452 "MD03BY.f"
	    dswap_(n, &r__[r_offset], &i__1, &dwork[1], &c__1);
#line 453 "MD03BY.f"
	    dtrsv_("Lower", "No transpose", "Non Unit", rank, &r__[r_offset], 
		    ldr, &rx[1], &c__1, (ftnlen)5, (ftnlen)12, (ftnlen)8);
#line 455 "MD03BY.f"
	    i__1 = *ldr + 1;
#line 455 "MD03BY.f"
	    dswap_(n, &r__[r_offset], &i__1, &dwork[1], &c__1);
#line 456 "MD03BY.f"
	    temp = dnrm2_(rank, &rx[1], &c__1);
#line 457 "MD03BY.f"
	    parc = fp / *delta / temp / temp;

/*              Depending on the sign of the function, update PARL */
/*              or PARU. */

#line 462 "MD03BY.f"
	    if (fp > 0.) {
#line 463 "MD03BY.f"
		parl = max(parl,*par);
#line 464 "MD03BY.f"
	    } else if (fp < 0.) {
#line 465 "MD03BY.f"
		paru = min(paru,*par);
#line 466 "MD03BY.f"
	    }

/*              Compute an improved estimate for PAR. */

/* Computing MAX */
#line 470 "MD03BY.f"
	    d__1 = parl, d__2 = *par + parc;
#line 470 "MD03BY.f"
	    *par = max(d__1,d__2);

/*              End of an iteration. */

#line 474 "MD03BY.f"
	    goto L70;
#line 475 "MD03BY.f"
	}
#line 476 "MD03BY.f"
    }

/*     Compute -R*P'*x = -R*z. */

#line 480 "MD03BY.f"
    if (econd && iter > 0) {

#line 482 "MD03BY.f"
	i__1 = *n;
#line 482 "MD03BY.f"
	for (j = 1; j <= i__1; ++j) {
#line 483 "MD03BY.f"
	    rx[j] = -dwork[*n + j];
#line 484 "MD03BY.f"
/* L110: */
#line 484 "MD03BY.f"
	}

#line 486 "MD03BY.f"
	dtrmv_("Upper", "NoTranspose", "NonUnit", n, &r__[r_offset], ldr, &rx[
		1], &c__1, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 488 "MD03BY.f"
    } else {

#line 490 "MD03BY.f"
	i__1 = *n;
#line 490 "MD03BY.f"
	for (j = 1; j <= i__1; ++j) {
#line 491 "MD03BY.f"
	    rx[j] = 0.;
#line 492 "MD03BY.f"
	    l = ipvt[j];
#line 493 "MD03BY.f"
	    d__1 = -x[l];
#line 493 "MD03BY.f"
	    daxpy_(&j, &d__1, &r__[j * r_dim1 + 1], &c__1, &rx[1], &c__1);
#line 494 "MD03BY.f"
/* L120: */
#line 494 "MD03BY.f"
	}

#line 496 "MD03BY.f"
    }

/*     Termination. If PAR = 0, set S. */

#line 500 "MD03BY.f"
    if (iter == 0) {
#line 501 "MD03BY.f"
	*par = 0.;

#line 503 "MD03BY.f"
	i__1 = *n - 1;
#line 503 "MD03BY.f"
	for (j = 1; j <= i__1; ++j) {
#line 504 "MD03BY.f"
	    dwork[j] = r__[j + j * r_dim1];
#line 505 "MD03BY.f"
	    i__2 = *n - j;
#line 505 "MD03BY.f"
	    dcopy_(&i__2, &r__[j + (j + 1) * r_dim1], ldr, &r__[j + 1 + j * 
		    r_dim1], &c__1);
#line 506 "MD03BY.f"
/* L130: */
#line 506 "MD03BY.f"
	}

#line 508 "MD03BY.f"
	dwork[*n] = r__[*n + *n * r_dim1];
#line 509 "MD03BY.f"
    }

#line 511 "MD03BY.f"
    return 0;

/* *** Last line of MD03BY *** */
} /* md03by_ */

