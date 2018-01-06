#line 1 "SB02MU.f"
/* SB02MU.f -- translated by f2c (version 20100827).
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

#line 1 "SB02MU.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b61 = 1.;

/* Subroutine */ int sb02mu_(char *dico, char *hinv, char *uplo, integer *n, 
	doublereal *a, integer *lda, doublereal *g, integer *ldg, doublereal *
	q, integer *ldq, doublereal *s, integer *lds, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen dico_len, 
	ftnlen hinv_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, q_dim1, q_offset, s_dim1, 
	    s_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, n2, nj, np1;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    static doublereal rcond, anorm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical lhinv, luplo;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), dlacpy_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgetri_(integer *, doublereal *, integer *, 
	    integer *, doublereal *, integer *, integer *), dgetrs_(char *, 
	    integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);
    static integer maxwrk;


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

/*     To construct the 2n-by-2n Hamiltonian or symplectic matrix S */
/*     associated to the linear-quadratic optimization problem, used to */
/*     solve the continuous- or discrete-time algebraic Riccati equation, */
/*     respectively. */

/*     For a continuous-time problem, S is defined by */

/*             (  A  -G ) */
/*         S = (        ),                                       (1) */
/*             ( -Q  -A') */

/*     and for a discrete-time problem by */

/*                 -1       -1 */
/*             (  A        A  *G     ) */
/*         S = (   -1           -1   ),                          (2) */
/*             ( QA     A' + Q*A  *G ) */

/*     or */

/*                       -T         -T */
/*             (  A + G*A  *Q   -G*A   ) */
/*         S = (      -T            -T ),                        (3) */
/*             (    -A  *Q         A   ) */

/*     where A, G, and Q are N-by-N matrices, with G and Q symmetric. */
/*     Matrix A must be nonsingular in the discrete-time case. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system as follows: */
/*             = 'C':  Continuous-time system; */
/*             = 'D':  Discrete-time system. */

/*     HINV    CHARACTER*1 */
/*             If DICO = 'D', specifies which of the matrices (2) or (3) */
/*             is constructed, as follows: */
/*             = 'D':  The matrix S in (2) is constructed; */
/*             = 'I':  The (inverse) matrix S in (3) is constructed. */
/*             HINV is not referenced if DICO = 'C'. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the matrices G and Q is */
/*             stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, G, and Q.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, if DICO = 'D', and INFO = 0, the leading N-by-N */
/*                                                     -1 */
/*             part of this array contains the matrix A  . */
/*             Otherwise, the array A is unchanged on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     G       (input) DOUBLE PRECISION array, dimension (LDG,N) */
/*             The leading N-by-N upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             must contain the upper triangular part or lower triangular */
/*             part, respectively, of the symmetric matrix G. The stricly */
/*             lower triangular part (if UPLO = 'U') or stricly upper */
/*             triangular part (if UPLO = 'L') is not referenced. */

/*     LDG     INTEGER */
/*             The leading dimension of array G.  LDG >= MAX(1,N). */

/*     Q       (input) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             The leading N-by-N upper triangular part (if UPLO = 'U') */
/*             or lower triangular part (if UPLO = 'L') of this array */
/*             must contain the upper triangular part or lower triangular */
/*             part, respectively, of the symmetric matrix Q. The stricly */
/*             lower triangular part (if UPLO = 'U') or stricly upper */
/*             triangular part (if UPLO = 'L') is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q.  LDQ >= MAX(1,N). */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,2*N) */
/*             If INFO = 0, the leading 2N-by-2N part of this array */
/*             contains the Hamiltonian or symplectic matrix of the */
/*             problem. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,2*N). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK; if DICO = 'D', DWORK(2) returns the reciprocal */
/*             condition number of the given matrix  A. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1          if DICO = 'C'; */
/*             LDWORK >= MAX(2,4*N) if DICO = 'D'. */
/*             For optimum performance LDWORK should be larger, if */
/*             DICO = 'D'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  if the leading i-by-i (1 <= i <= N) upper triangular */
/*                   submatrix of A is singular in discrete-time case; */
/*             = N+1:  if matrix A is numerically singular in discrete- */
/*                   time case. */

/*     METHOD */

/*     For a continuous-time problem, the 2n-by-2n Hamiltonian matrix (1) */
/*     is constructed. */
/*     For a discrete-time problem, the 2n-by-2n symplectic matrix (2) or */
/*     (3) - the inverse of the matrix in (2) - is constructed. */

/*     NUMERICAL ASPECTS */

/*     The discrete-time case needs the inverse of the matrix A, hence */
/*     the routine should not be used when A is ill-conditioned. */
/*                               3 */
/*     The algorithm requires 0(n ) floating point operations in the */
/*     discrete-time case. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, closed loop system, continuous-time */
/*     system, discrete-time system, optimal regulator, Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 203 "SB02MU.f"
    /* Parameter adjustments */
#line 203 "SB02MU.f"
    a_dim1 = *lda;
#line 203 "SB02MU.f"
    a_offset = 1 + a_dim1;
#line 203 "SB02MU.f"
    a -= a_offset;
#line 203 "SB02MU.f"
    g_dim1 = *ldg;
#line 203 "SB02MU.f"
    g_offset = 1 + g_dim1;
#line 203 "SB02MU.f"
    g -= g_offset;
#line 203 "SB02MU.f"
    q_dim1 = *ldq;
#line 203 "SB02MU.f"
    q_offset = 1 + q_dim1;
#line 203 "SB02MU.f"
    q -= q_offset;
#line 203 "SB02MU.f"
    s_dim1 = *lds;
#line 203 "SB02MU.f"
    s_offset = 1 + s_dim1;
#line 203 "SB02MU.f"
    s -= s_offset;
#line 203 "SB02MU.f"
    --iwork;
#line 203 "SB02MU.f"
    --dwork;
#line 203 "SB02MU.f"

#line 203 "SB02MU.f"
    /* Function Body */
#line 203 "SB02MU.f"
    *info = 0;
#line 204 "SB02MU.f"
    n2 = *n + *n;
#line 205 "SB02MU.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 206 "SB02MU.f"
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 207 "SB02MU.f"
    if (discr) {
#line 208 "SB02MU.f"
	lhinv = lsame_(hinv, "D", (ftnlen)1, (ftnlen)1);
#line 209 "SB02MU.f"
    } else {
#line 210 "SB02MU.f"
	lhinv = FALSE_;
#line 211 "SB02MU.f"
    }

/*     Test the input scalar arguments. */

#line 215 "SB02MU.f"
    if (! discr && ! lsame_(dico, "C", (ftnlen)1, (ftnlen)1)) {
#line 216 "SB02MU.f"
	*info = -1;
#line 217 "SB02MU.f"
    } else if (discr) {
#line 218 "SB02MU.f"
	if (! lhinv && ! lsame_(hinv, "I", (ftnlen)1, (ftnlen)1)) {
#line 218 "SB02MU.f"
	    *info = -2;
#line 218 "SB02MU.f"
	}
#line 220 "SB02MU.f"
    }
#line 221 "SB02MU.f"
    if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 222 "SB02MU.f"
	*info = -3;
#line 223 "SB02MU.f"
    } else if (*n < 0) {
#line 224 "SB02MU.f"
	*info = -4;
#line 225 "SB02MU.f"
    } else if (*lda < max(1,*n)) {
#line 226 "SB02MU.f"
	*info = -6;
#line 227 "SB02MU.f"
    } else if (*ldg < max(1,*n)) {
#line 228 "SB02MU.f"
	*info = -8;
#line 229 "SB02MU.f"
    } else if (*ldq < max(1,*n)) {
#line 230 "SB02MU.f"
	*info = -10;
#line 231 "SB02MU.f"
    } else if (*lds < max(1,n2)) {
#line 232 "SB02MU.f"
	*info = -12;
#line 233 "SB02MU.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 233 "SB02MU.f"
	i__1 = 2, i__2 = *n << 2;
#line 233 "SB02MU.f"
	if (*ldwork < 1 || discr && *ldwork < max(i__1,i__2)) {
#line 235 "SB02MU.f"
	    *info = -15;
#line 236 "SB02MU.f"
	}
#line 236 "SB02MU.f"
    }

#line 238 "SB02MU.f"
    if (*info != 0) {

/*        Error return. */

#line 242 "SB02MU.f"
	i__1 = -(*info);
#line 242 "SB02MU.f"
	xerbla_("SB02MU", &i__1, (ftnlen)6);
#line 243 "SB02MU.f"
	return 0;
#line 244 "SB02MU.f"
    }

/*     Quick return if possible. */

#line 248 "SB02MU.f"
    if (*n == 0) {
#line 249 "SB02MU.f"
	dwork[1] = 1.;
#line 250 "SB02MU.f"
	if (discr) {
#line 250 "SB02MU.f"
	    dwork[2] = 1.;
#line 250 "SB02MU.f"
	}
#line 251 "SB02MU.f"
	return 0;
#line 252 "SB02MU.f"
    }

/*     The code tries to exploit data locality as much as possible. */

#line 256 "SB02MU.f"
    if (! lhinv) {
#line 257 "SB02MU.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &s[s_offset], lds, (ftnlen)4)
		;

/*        Construct Hamiltonian matrix in the continuous-time case, or */
/*        prepare symplectic matrix in (3) in the discrete-time case: */

/*        Construct full Q in S(N+1:2*N,1:N) and change the sign, and */
/*        construct full G in S(1:N,N+1:2*N) and change the sign. */

#line 265 "SB02MU.f"
	i__1 = *n;
#line 265 "SB02MU.f"
	for (j = 1; j <= i__1; ++j) {
#line 266 "SB02MU.f"
	    nj = *n + j;
#line 267 "SB02MU.f"
	    if (luplo) {

#line 269 "SB02MU.f"
		i__2 = j;
#line 269 "SB02MU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 270 "SB02MU.f"
		    s[*n + i__ + j * s_dim1] = -q[i__ + j * q_dim1];
#line 271 "SB02MU.f"
/* L20: */
#line 271 "SB02MU.f"
		}

#line 273 "SB02MU.f"
		i__2 = *n;
#line 273 "SB02MU.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 274 "SB02MU.f"
		    s[*n + i__ + j * s_dim1] = -q[j + i__ * q_dim1];
#line 275 "SB02MU.f"
/* L40: */
#line 275 "SB02MU.f"
		}

#line 277 "SB02MU.f"
		i__2 = j;
#line 277 "SB02MU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 278 "SB02MU.f"
		    s[i__ + nj * s_dim1] = -g[i__ + j * g_dim1];
#line 279 "SB02MU.f"
/* L60: */
#line 279 "SB02MU.f"
		}

#line 281 "SB02MU.f"
		i__2 = *n;
#line 281 "SB02MU.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 282 "SB02MU.f"
		    s[i__ + nj * s_dim1] = -g[j + i__ * g_dim1];
#line 283 "SB02MU.f"
/* L80: */
#line 283 "SB02MU.f"
		}

#line 285 "SB02MU.f"
	    } else {

#line 287 "SB02MU.f"
		i__2 = j - 1;
#line 287 "SB02MU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 288 "SB02MU.f"
		    s[*n + i__ + j * s_dim1] = -q[j + i__ * q_dim1];
#line 289 "SB02MU.f"
/* L100: */
#line 289 "SB02MU.f"
		}

#line 291 "SB02MU.f"
		i__2 = *n;
#line 291 "SB02MU.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 292 "SB02MU.f"
		    s[*n + i__ + j * s_dim1] = -q[i__ + j * q_dim1];
#line 293 "SB02MU.f"
/* L120: */
#line 293 "SB02MU.f"
		}

#line 295 "SB02MU.f"
		i__2 = j - 1;
#line 295 "SB02MU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 296 "SB02MU.f"
		    s[i__ + nj * s_dim1] = -g[j + i__ * g_dim1];
#line 297 "SB02MU.f"
/* L140: */
#line 297 "SB02MU.f"
		}

#line 299 "SB02MU.f"
		i__2 = *n;
#line 299 "SB02MU.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 300 "SB02MU.f"
		    s[i__ + nj * s_dim1] = -g[i__ + j * g_dim1];
#line 301 "SB02MU.f"
/* L180: */
#line 301 "SB02MU.f"
		}

#line 303 "SB02MU.f"
	    }
#line 304 "SB02MU.f"
/* L200: */
#line 304 "SB02MU.f"
	}

#line 306 "SB02MU.f"
	if (! discr) {

#line 308 "SB02MU.f"
	    i__1 = *n;
#line 308 "SB02MU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 309 "SB02MU.f"
		nj = *n + j;

#line 311 "SB02MU.f"
		i__2 = *n;
#line 311 "SB02MU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 312 "SB02MU.f"
		    s[*n + i__ + nj * s_dim1] = -a[j + i__ * a_dim1];
#line 313 "SB02MU.f"
/* L220: */
#line 313 "SB02MU.f"
		}

#line 315 "SB02MU.f"
/* L240: */
#line 315 "SB02MU.f"
	    }

#line 317 "SB02MU.f"
	    dwork[1] = 1.;
#line 318 "SB02MU.f"
	}
#line 319 "SB02MU.f"
    }

#line 321 "SB02MU.f"
    if (discr) {

/*        Construct the symplectic matrix (2) or (3) in the discrete-time */
/*        case. */

/*        Compute workspace. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*        minimal amount of workspace needed at that point in the code, */
/*        as well as the preferred amount for good performance. */
/*        NB refers to the optimal block size for the immediately */
/*        following subroutine, as returned by ILAENV.) */

/* Computing MAX */
#line 333 "SB02MU.f"
	i__1 = *n << 2, i__2 = *n * ilaenv_(&c__1, "DGETRI", " ", n, &c_n1, &
		c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 333 "SB02MU.f"
	maxwrk = max(i__1,i__2);
#line 335 "SB02MU.f"
	np1 = *n + 1;

#line 337 "SB02MU.f"
	if (lhinv) {

/*           Put  A'  in  S(N+1:2*N,N+1:2*N). */

#line 341 "SB02MU.f"
	    i__1 = *n;
#line 341 "SB02MU.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 342 "SB02MU.f"
		dcopy_(n, &a[i__ + a_dim1], lda, &s[np1 + (*n + i__) * s_dim1]
			, &c__1);
#line 343 "SB02MU.f"
/* L260: */
#line 343 "SB02MU.f"
	    }

#line 345 "SB02MU.f"
	}

/*        Compute the norm of the matrix A. */

#line 349 "SB02MU.f"
	anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)
		6);

/*        Compute the LU factorization of A. */

#line 353 "SB02MU.f"
	dgetrf_(n, n, &a[a_offset], lda, &iwork[1], info);

/*        Return if INFO is non-zero. */

#line 357 "SB02MU.f"
	if (*info > 0) {
#line 358 "SB02MU.f"
	    dwork[2] = 0.;
#line 359 "SB02MU.f"
	    return 0;
#line 360 "SB02MU.f"
	}

/*        Compute the reciprocal of the condition number of A. */
/*        Workspace: need 4*N. */

#line 365 "SB02MU.f"
	dgecon_("1-norm", n, &a[a_offset], lda, &anorm, &rcond, &dwork[1], &
		iwork[np1], info, (ftnlen)6);

/*        Return if the matrix is singular to working precision. */

#line 370 "SB02MU.f"
	if (rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 371 "SB02MU.f"
	    *info = *n + 1;
#line 372 "SB02MU.f"
	    dwork[2] = rcond;
#line 373 "SB02MU.f"
	    return 0;
#line 374 "SB02MU.f"
	}

#line 376 "SB02MU.f"
	if (lhinv) {

/*           Compute S in (2). */

/*           Construct full Q in S(N+1:2*N,1:N). */

#line 382 "SB02MU.f"
	    if (luplo) {
#line 383 "SB02MU.f"
		i__1 = *n - 1;
#line 383 "SB02MU.f"
		for (j = 1; j <= i__1; ++j) {
#line 384 "SB02MU.f"
		    dcopy_(&j, &q[j * q_dim1 + 1], &c__1, &s[np1 + j * s_dim1]
			    , &c__1);
#line 385 "SB02MU.f"
		    i__2 = *n - j;
#line 385 "SB02MU.f"
		    dcopy_(&i__2, &q[j + (j + 1) * q_dim1], ldq, &s[np1 + j + 
			    j * s_dim1], &c__1);
#line 386 "SB02MU.f"
/* L270: */
#line 386 "SB02MU.f"
		}
#line 387 "SB02MU.f"
		dcopy_(n, &q[*n * q_dim1 + 1], &c__1, &s[np1 + *n * s_dim1], &
			c__1);
#line 388 "SB02MU.f"
	    } else {
#line 389 "SB02MU.f"
		dcopy_(n, &q[q_dim1 + 1], &c__1, &s[np1 + s_dim1], &c__1);
#line 390 "SB02MU.f"
		i__1 = *n;
#line 390 "SB02MU.f"
		for (j = 2; j <= i__1; ++j) {
#line 391 "SB02MU.f"
		    i__2 = j - 1;
#line 391 "SB02MU.f"
		    dcopy_(&i__2, &q[j + q_dim1], ldq, &s[np1 + j * s_dim1], &
			    c__1);
#line 392 "SB02MU.f"
		    i__2 = *n - j + 1;
#line 392 "SB02MU.f"
		    dcopy_(&i__2, &q[j + j * q_dim1], &c__1, &s[*n + j + j * 
			    s_dim1], &c__1);
#line 393 "SB02MU.f"
/* L280: */
#line 393 "SB02MU.f"
		}
#line 394 "SB02MU.f"
	    }

/*           Compute the solution matrix  X  of the system  X*A = Q  by */
/*                                                                    -1 */
/*           solving  A'*X' = Q and transposing the result to get  Q*A  . */

#line 400 "SB02MU.f"
	    dgetrs_("Transpose", n, n, &a[a_offset], lda, &iwork[1], &s[np1 + 
		    s_dim1], lds, info, (ftnlen)9);

#line 403 "SB02MU.f"
	    i__1 = *n - 1;
#line 403 "SB02MU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 404 "SB02MU.f"
		i__2 = *n - j;
#line 404 "SB02MU.f"
		dswap_(&i__2, &s[np1 + j + j * s_dim1], &c__1, &s[*n + j + (j 
			+ 1) * s_dim1], lds);
#line 405 "SB02MU.f"
/* L300: */
#line 405 "SB02MU.f"
	    }

/*           Construct full G in S(1:N,N+1:2*N). */

#line 409 "SB02MU.f"
	    if (luplo) {
#line 410 "SB02MU.f"
		i__1 = *n - 1;
#line 410 "SB02MU.f"
		for (j = 1; j <= i__1; ++j) {
#line 411 "SB02MU.f"
		    dcopy_(&j, &g[j * g_dim1 + 1], &c__1, &s[(*n + j) * 
			    s_dim1 + 1], &c__1);
#line 412 "SB02MU.f"
		    i__2 = *n - j;
#line 412 "SB02MU.f"
		    dcopy_(&i__2, &g[j + (j + 1) * g_dim1], ldg, &s[j + 1 + (*
			    n + j) * s_dim1], &c__1);
#line 413 "SB02MU.f"
/* L310: */
#line 413 "SB02MU.f"
		}
#line 414 "SB02MU.f"
		dcopy_(n, &g[*n * g_dim1 + 1], &c__1, &s[n2 * s_dim1 + 1], &
			c__1);
#line 415 "SB02MU.f"
	    } else {
#line 416 "SB02MU.f"
		dcopy_(n, &g[g_dim1 + 1], &c__1, &s[np1 * s_dim1 + 1], &c__1);
#line 417 "SB02MU.f"
		i__1 = *n;
#line 417 "SB02MU.f"
		for (j = 2; j <= i__1; ++j) {
#line 418 "SB02MU.f"
		    i__2 = j - 1;
#line 418 "SB02MU.f"
		    dcopy_(&i__2, &g[j + g_dim1], ldg, &s[(*n + j) * s_dim1 + 
			    1], &c__1);
#line 419 "SB02MU.f"
		    i__2 = *n - j + 1;
#line 419 "SB02MU.f"
		    dcopy_(&i__2, &g[j + j * g_dim1], &c__1, &s[j + (*n + j) *
			     s_dim1], &c__1);
#line 420 "SB02MU.f"
/* L320: */
#line 420 "SB02MU.f"
		}
#line 421 "SB02MU.f"
	    }
/*                            -1 */
/*           Compute  A' + Q*A  *G  in  S(N+1:2N,N+1:2N). */

#line 425 "SB02MU.f"
	    dgemm_("No transpose", "No transpose", n, n, n, &c_b61, &s[np1 + 
		    s_dim1], lds, &s[np1 * s_dim1 + 1], lds, &c_b61, &s[np1 + 
		    np1 * s_dim1], lds, (ftnlen)12, (ftnlen)12);

/*           Compute the solution matrix  Y  of the system  A*Y = G. */

#line 431 "SB02MU.f"
	    dgetrs_("No transpose", n, n, &a[a_offset], lda, &iwork[1], &s[
		    np1 * s_dim1 + 1], lds, info, (ftnlen)12);

/*           Compute the inverse of  A  in situ. */
/*           Workspace: need N;  prefer N*NB. */

#line 437 "SB02MU.f"
	    dgetri_(n, &a[a_offset], lda, &iwork[1], &dwork[1], ldwork, info);
/*                  -1 */
/*           Copy  A    in  S(1:N,1:N). */

#line 441 "SB02MU.f"
	    dlacpy_("Full", n, n, &a[a_offset], lda, &s[s_offset], lds, (
		    ftnlen)4);

#line 443 "SB02MU.f"
	} else {

/*           Compute S in (3) using the already prepared part. */

/*           Compute the solution matrix  X'  of the system  A*X' = -G */
/*                                                       -T */
/*           and transpose the result to obtain  X = -G*A  . */

#line 451 "SB02MU.f"
	    dgetrs_("No transpose", n, n, &a[a_offset], lda, &iwork[1], &s[
		    np1 * s_dim1 + 1], lds, info, (ftnlen)12);

#line 454 "SB02MU.f"
	    i__1 = *n - 1;
#line 454 "SB02MU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 455 "SB02MU.f"
		i__2 = *n - j;
#line 455 "SB02MU.f"
		dswap_(&i__2, &s[j + 1 + (*n + j) * s_dim1], &c__1, &s[j + (
			np1 + j) * s_dim1], lds);
#line 456 "SB02MU.f"
/* L340: */
#line 456 "SB02MU.f"
	    }
/*                           -T */
/*           Compute  A + G*A  *Q  in  S(1:N,1:N). */

#line 460 "SB02MU.f"
	    dgemm_("No transpose", "No transpose", n, n, n, &c_b61, &s[np1 * 
		    s_dim1 + 1], lds, &s[np1 + s_dim1], lds, &c_b61, &s[
		    s_offset], lds, (ftnlen)12, (ftnlen)12);

/*           Compute the solution matrix  Y  of the system  A'*Y = -Q. */

#line 465 "SB02MU.f"
	    dgetrs_("Transpose", n, n, &a[a_offset], lda, &iwork[1], &s[np1 + 
		    s_dim1], lds, info, (ftnlen)9);

/*           Compute the inverse of  A  in situ. */
/*           Workspace: need N;  prefer N*NB. */

#line 471 "SB02MU.f"
	    dgetri_(n, &a[a_offset], lda, &iwork[1], &dwork[1], ldwork, info);
/*                  -T */
/*           Copy  A    in  S(N+1:2N,N+1:2N). */

#line 475 "SB02MU.f"
	    i__1 = *n;
#line 475 "SB02MU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 476 "SB02MU.f"
		dcopy_(n, &a[j + a_dim1], lda, &s[np1 + (*n + j) * s_dim1], &
			c__1);
#line 477 "SB02MU.f"
/* L360: */
#line 477 "SB02MU.f"
	    }

#line 479 "SB02MU.f"
	}
#line 480 "SB02MU.f"
	dwork[1] = (doublereal) maxwrk;
#line 481 "SB02MU.f"
	dwork[2] = rcond;
#line 482 "SB02MU.f"
    }

/* *** Last line of SB02MU *** */
#line 485 "SB02MU.f"
    return 0;
} /* sb02mu_ */

