#line 1 "SB02RU.f"
/* SB02RU.f -- translated by f2c (version 20100827).
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

#line 1 "SB02RU.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b33 = 0.;
static doublereal c_b34 = 1.;
static doublereal c_b57 = -1.;

/* Subroutine */ int sb02ru_(char *dico, char *hinv, char *trana, char *uplo, 
	integer *n, doublereal *a, integer *lda, doublereal *g, integer *ldg, 
	doublereal *q, integer *ldq, doublereal *s, integer *lds, integer *
	iwork, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	dico_len, ftnlen hinv_len, ftnlen trana_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, q_dim1, q_offset, s_dim1, 
	    s_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, n2, nj, np1;
    static doublereal temp;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ma02ed_(char *, integer *, doublereal *, integer *, ftnlen), 
	    mb02pd_(char *, char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, char *, doublereal 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, doublereal *
	    , integer *, ftnlen, ftnlen, ftnlen), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char equed[1];
    static logical discr;
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical lhinv, luplo;
    static doublereal rconda;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static char tranat[1];
    static logical notrna;
    static doublereal pivotg;


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

/*             ( op(A)   -G    ) */
/*         S = (               ),                                     (1) */
/*             (  -Q   -op(A)' ) */

/*     and for a discrete-time problem by */

/*                     -1              -1 */
/*             (  op(A)           op(A)  *G       ) */
/*         S = (        -1                   -1   ),                  (2) */
/*             ( Q*op(A)     op(A)' + Q*op(A)  *G ) */

/*     or */
/*                              -T             -T */
/*             ( op(A) + G*op(A)  *Q   -G*op(A)   ) */
/*         S = (           -T                 -T  ),                  (3) */
/*             (     -op(A)  *Q          op(A)    ) */

/*     where op(A) = A or A' (A**T), A, G, and Q are n-by-n matrices, */
/*     with G and Q symmetric. Matrix A must be nonsingular in the */
/*     discrete-time case. */

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

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used, as follows: */
/*             = 'N':  op(A) = A    (No transpose); */
/*             = 'T':  op(A) = A**T (Transpose); */
/*             = 'C':  op(A) = A**T (Conjugate transpose = Transpose). */

/*     UPLO    CHARACTER*1 */
/*             Specifies which triangle of the matrices G and Q is */
/*             stored, as follows: */
/*             = 'U':  Upper triangle is stored; */
/*             = 'L':  Lower triangle is stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, G, and Q.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, the leading N-by-N upper triangular part (if */
/*             UPLO = 'U') or lower triangular part (if UPLO = 'L') of */
/*             this array must contain the upper triangular part or lower */
/*             triangular part, respectively, of the symmetric matrix G. */
/*             On exit, if DICO = 'D', the leading N-by-N part of this */
/*             array contains the symmetric matrix G fully stored. */
/*             If DICO = 'C', this array is not modified on exit, and the */
/*             strictly lower triangular part (if UPLO = 'U') or strictly */
/*             upper triangular part (if UPLO = 'L') is not referenced. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, the leading N-by-N upper triangular part (if */
/*             UPLO = 'U') or lower triangular part (if UPLO = 'L') of */
/*             this array must contain the upper triangular part or lower */
/*             triangular part, respectively, of the symmetric matrix Q. */
/*             On exit, if DICO = 'D', the leading N-by-N part of this */
/*             array contains the symmetric matrix Q fully stored. */
/*             If DICO = 'C', this array is not modified on exit, and the */
/*             strictly lower triangular part (if UPLO = 'U') or strictly */
/*             upper triangular part (if UPLO = 'L') is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= MAX(1,N). */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,2*N) */
/*             If INFO = 0, the leading 2N-by-2N part of this array */
/*             contains the Hamiltonian or symplectic matrix of the */
/*             problem. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S.  LDS >= MAX(1,2*N). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), where */
/*             LIWORK >= 0,   if DICO = 'C'; */
/*             LIWORK >= 2*N, if DICO = 'D'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if DICO = 'D', DWORK(1) returns the reciprocal */
/*             condition number  RCOND  of the given matrix  A,  and */
/*             DWORK(2) returns the reciprocal pivot growth factor */
/*             norm(A)/norm(U) (see SLICOT Library routine MB02PD). */
/*             If DWORK(2) is much less than 1, then the computed  S */
/*             and  RCOND  could be unreliable. If 0 < INFO <= N, then */
/*             DWORK(2) contains the reciprocal pivot growth factor for */
/*             the leading INFO columns of  A. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 0,          if DICO = 'C'; */
/*             LDWORK >= MAX(2,6*N), if DICO = 'D'. */

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

/*     FURTHER COMMENTS */

/*     This routine is a functionally extended and with improved accuracy */
/*     version of the SLICOT Library routine SB02MU. Transposed problems */
/*     can be dealt with as well. The LU factorization of  op(A)  (with */
/*     no equilibration) and iterative refinement are used for solving */
/*     the various linear algebraic systems involved. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 1999. */

/*     REVISIONS */

/*     - */

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

#line 223 "SB02RU.f"
    /* Parameter adjustments */
#line 223 "SB02RU.f"
    a_dim1 = *lda;
#line 223 "SB02RU.f"
    a_offset = 1 + a_dim1;
#line 223 "SB02RU.f"
    a -= a_offset;
#line 223 "SB02RU.f"
    g_dim1 = *ldg;
#line 223 "SB02RU.f"
    g_offset = 1 + g_dim1;
#line 223 "SB02RU.f"
    g -= g_offset;
#line 223 "SB02RU.f"
    q_dim1 = *ldq;
#line 223 "SB02RU.f"
    q_offset = 1 + q_dim1;
#line 223 "SB02RU.f"
    q -= q_offset;
#line 223 "SB02RU.f"
    s_dim1 = *lds;
#line 223 "SB02RU.f"
    s_offset = 1 + s_dim1;
#line 223 "SB02RU.f"
    s -= s_offset;
#line 223 "SB02RU.f"
    --iwork;
#line 223 "SB02RU.f"
    --dwork;
#line 223 "SB02RU.f"

#line 223 "SB02RU.f"
    /* Function Body */
#line 223 "SB02RU.f"
    n2 = *n + *n;
#line 224 "SB02RU.f"
    *info = 0;
#line 225 "SB02RU.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 226 "SB02RU.f"
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 227 "SB02RU.f"
    notrna = lsame_(trana, "N", (ftnlen)1, (ftnlen)1);
#line 228 "SB02RU.f"
    if (discr) {
#line 228 "SB02RU.f"
	lhinv = lsame_(hinv, "D", (ftnlen)1, (ftnlen)1);
#line 228 "SB02RU.f"
    }

/*     Test the input scalar arguments. */

#line 233 "SB02RU.f"
    if (! discr && ! lsame_(dico, "C", (ftnlen)1, (ftnlen)1)) {
#line 234 "SB02RU.f"
	*info = -1;
#line 235 "SB02RU.f"
    } else if (discr) {
#line 236 "SB02RU.f"
	if (! lhinv && ! lsame_(hinv, "I", (ftnlen)1, (ftnlen)1)) {
#line 236 "SB02RU.f"
	    *info = -2;
#line 236 "SB02RU.f"
	}
#line 238 "SB02RU.f"
    } else if (*info == 0) {
#line 239 "SB02RU.f"
	if (! notrna && ! lsame_(trana, "T", (ftnlen)1, (ftnlen)1) && ! 
		lsame_(trana, "C", (ftnlen)1, (ftnlen)1)) {
#line 241 "SB02RU.f"
	    *info = -3;
#line 242 "SB02RU.f"
	} else if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 243 "SB02RU.f"
	    *info = -4;
#line 244 "SB02RU.f"
	} else if (*n < 0) {
#line 245 "SB02RU.f"
	    *info = -5;
#line 246 "SB02RU.f"
	} else if (*lda < max(1,*n)) {
#line 247 "SB02RU.f"
	    *info = -7;
#line 248 "SB02RU.f"
	} else if (*ldg < max(1,*n)) {
#line 249 "SB02RU.f"
	    *info = -9;
#line 250 "SB02RU.f"
	} else if (*ldq < max(1,*n)) {
#line 251 "SB02RU.f"
	    *info = -11;
#line 252 "SB02RU.f"
	} else if (*lds < max(1,n2)) {
#line 253 "SB02RU.f"
	    *info = -13;
#line 254 "SB02RU.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
#line 254 "SB02RU.f"
	    i__1 = 2, i__2 = *n * 6;
#line 254 "SB02RU.f"
	    if (*ldwork < 0 || discr && *ldwork < max(i__1,i__2)) {
#line 256 "SB02RU.f"
		*info = -16;
#line 257 "SB02RU.f"
	    }
#line 257 "SB02RU.f"
	}
#line 258 "SB02RU.f"
    }

#line 260 "SB02RU.f"
    if (*info != 0) {

/*        Error return. */

#line 264 "SB02RU.f"
	i__1 = -(*info);
#line 264 "SB02RU.f"
	xerbla_("SB02RU", &i__1, (ftnlen)6);
#line 265 "SB02RU.f"
	return 0;
#line 266 "SB02RU.f"
    }

/*     Quick return if possible. */

#line 270 "SB02RU.f"
    if (*n == 0) {
#line 271 "SB02RU.f"
	if (discr) {
#line 272 "SB02RU.f"
	    dwork[1] = 1.;
#line 273 "SB02RU.f"
	    dwork[2] = 1.;
#line 274 "SB02RU.f"
	}
#line 275 "SB02RU.f"
	return 0;
#line 276 "SB02RU.f"
    }

/*     The code tries to exploit data locality as much as possible, */
/*     assuming that LDS is greater than LDA, LDQ, and/or LDG. */

#line 281 "SB02RU.f"
    if (! discr) {

/*        Continuous-time case: Construct Hamiltonian matrix column-wise. */

/*        Copy op(A) in S(1:N,1:N), and construct full Q */
/*        in S(N+1:2*N,1:N) and change the sign. */

#line 288 "SB02RU.f"
	i__1 = *n;
#line 288 "SB02RU.f"
	for (j = 1; j <= i__1; ++j) {
#line 289 "SB02RU.f"
	    if (notrna) {
#line 290 "SB02RU.f"
		dcopy_(n, &a[j * a_dim1 + 1], &c__1, &s[j * s_dim1 + 1], &
			c__1);
#line 291 "SB02RU.f"
	    } else {
#line 292 "SB02RU.f"
		dcopy_(n, &a[j + a_dim1], lda, &s[j * s_dim1 + 1], &c__1);
#line 293 "SB02RU.f"
	    }

#line 295 "SB02RU.f"
	    if (luplo) {

#line 297 "SB02RU.f"
		i__2 = j;
#line 297 "SB02RU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 298 "SB02RU.f"
		    s[*n + i__ + j * s_dim1] = -q[i__ + j * q_dim1];
#line 299 "SB02RU.f"
/* L20: */
#line 299 "SB02RU.f"
		}

#line 301 "SB02RU.f"
		i__2 = *n;
#line 301 "SB02RU.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 302 "SB02RU.f"
		    s[*n + i__ + j * s_dim1] = -q[j + i__ * q_dim1];
#line 303 "SB02RU.f"
/* L40: */
#line 303 "SB02RU.f"
		}

#line 305 "SB02RU.f"
	    } else {

#line 307 "SB02RU.f"
		i__2 = j - 1;
#line 307 "SB02RU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 308 "SB02RU.f"
		    s[*n + i__ + j * s_dim1] = -q[j + i__ * q_dim1];
#line 309 "SB02RU.f"
/* L60: */
#line 309 "SB02RU.f"
		}

#line 311 "SB02RU.f"
		i__2 = *n;
#line 311 "SB02RU.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 312 "SB02RU.f"
		    s[*n + i__ + j * s_dim1] = -q[i__ + j * q_dim1];
#line 313 "SB02RU.f"
/* L80: */
#line 313 "SB02RU.f"
		}

#line 315 "SB02RU.f"
	    }
#line 316 "SB02RU.f"
/* L100: */
#line 316 "SB02RU.f"
	}

/*        Construct full G in S(1:N,N+1:2*N) and change the sign, and */
/*        construct -op(A)' in S(N+1:2*N,N+1:2*N). */

#line 321 "SB02RU.f"
	i__1 = *n;
#line 321 "SB02RU.f"
	for (j = 1; j <= i__1; ++j) {
#line 322 "SB02RU.f"
	    nj = *n + j;
#line 323 "SB02RU.f"
	    if (luplo) {

#line 325 "SB02RU.f"
		i__2 = j;
#line 325 "SB02RU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 326 "SB02RU.f"
		    s[i__ + nj * s_dim1] = -g[i__ + j * g_dim1];
#line 327 "SB02RU.f"
/* L120: */
#line 327 "SB02RU.f"
		}

#line 329 "SB02RU.f"
		i__2 = *n;
#line 329 "SB02RU.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 330 "SB02RU.f"
		    s[i__ + nj * s_dim1] = -g[j + i__ * g_dim1];
#line 331 "SB02RU.f"
/* L140: */
#line 331 "SB02RU.f"
		}

#line 333 "SB02RU.f"
	    } else {

#line 335 "SB02RU.f"
		i__2 = j - 1;
#line 335 "SB02RU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 336 "SB02RU.f"
		    s[i__ + nj * s_dim1] = -g[j + i__ * g_dim1];
#line 337 "SB02RU.f"
/* L160: */
#line 337 "SB02RU.f"
		}

#line 339 "SB02RU.f"
		i__2 = *n;
#line 339 "SB02RU.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 340 "SB02RU.f"
		    s[i__ + nj * s_dim1] = -g[i__ + j * g_dim1];
#line 341 "SB02RU.f"
/* L180: */
#line 341 "SB02RU.f"
		}

#line 343 "SB02RU.f"
	    }

#line 345 "SB02RU.f"
	    if (notrna) {

#line 347 "SB02RU.f"
		i__2 = *n;
#line 347 "SB02RU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 348 "SB02RU.f"
		    s[*n + i__ + nj * s_dim1] = -a[j + i__ * a_dim1];
#line 349 "SB02RU.f"
/* L200: */
#line 349 "SB02RU.f"
		}

#line 351 "SB02RU.f"
	    } else {

#line 353 "SB02RU.f"
		i__2 = *n;
#line 353 "SB02RU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 354 "SB02RU.f"
		    s[*n + i__ + nj * s_dim1] = -a[i__ + j * a_dim1];
#line 355 "SB02RU.f"
/* L220: */
#line 355 "SB02RU.f"
		}

#line 357 "SB02RU.f"
	    }
#line 358 "SB02RU.f"
/* L240: */
#line 358 "SB02RU.f"
	}

#line 360 "SB02RU.f"
    } else {

/*        Discrete-time case: Construct the symplectic matrix (2) or (3). */

/*        Fill in the remaining triangles of the symmetric matrices Q */
/*        and G. */

#line 367 "SB02RU.f"
	ma02ed_(uplo, n, &q[q_offset], ldq, (ftnlen)1);
#line 368 "SB02RU.f"
	ma02ed_(uplo, n, &g[g_offset], ldg, (ftnlen)1);

/*        Prepare the construction of S in (2) or (3). */

#line 372 "SB02RU.f"
	np1 = *n + 1;
#line 373 "SB02RU.f"
	if (notrna) {
#line 374 "SB02RU.f"
	    *(unsigned char *)tranat = 'T';
#line 375 "SB02RU.f"
	} else {
#line 376 "SB02RU.f"
	    *(unsigned char *)tranat = 'N';
#line 377 "SB02RU.f"
	}

/*        Solve  op(A)'*X = Q  in  S(N+1:2*N,1:N),  using the LU */
/*        factorization of  op(A),  obtained in  S(1:N,1:N),  and */
/*        iterative refinement. No equilibration of  A  is used. */
/*        Workspace:  6*N. */

#line 384 "SB02RU.f"
	mb02pd_("No equilibration", tranat, n, n, &a[a_offset], lda, &s[
		s_offset], lds, &iwork[1], equed, &dwork[1], &dwork[1], &q[
		q_offset], ldq, &s[np1 + s_dim1], lds, &rcond, &dwork[1], &
		dwork[np1], &iwork[np1], &dwork[n2 + 1], info, (ftnlen)16, (
		ftnlen)1, (ftnlen)1);

/*        Return if the matrix is exactly singular or singular to */
/*        working precision. */

#line 392 "SB02RU.f"
	if (*info > 0) {
#line 393 "SB02RU.f"
	    dwork[1] = rcond;
#line 394 "SB02RU.f"
	    dwork[2] = dwork[n2 + 1];
#line 395 "SB02RU.f"
	    return 0;
#line 396 "SB02RU.f"
	}

#line 398 "SB02RU.f"
	rconda = rcond;
#line 399 "SB02RU.f"
	pivotg = dwork[n2 + 1];

#line 401 "SB02RU.f"
	if (lhinv) {

/*           Complete the construction of S in (2). */

/*           Transpose  X  in-situ. */

#line 407 "SB02RU.f"
	    i__1 = *n - 1;
#line 407 "SB02RU.f"
	    for (j = 1; j <= i__1; ++j) {
#line 408 "SB02RU.f"
		i__2 = *n - j;
#line 408 "SB02RU.f"
		dswap_(&i__2, &s[np1 + j + j * s_dim1], &c__1, &s[*n + j + (j 
			+ 1) * s_dim1], lds);
#line 409 "SB02RU.f"
/* L260: */
#line 409 "SB02RU.f"
	    }

/*           Solve  op(A)*X = I_n  in  S(N+1:2*N,N+1:2*N),  using the LU */
/*           factorization of  op(A),  computed in  S(1:N,1:N),  and */
/*           iterative refinement. */

#line 415 "SB02RU.f"
	    dlaset_("Full", n, n, &c_b33, &c_b34, &s[np1 * s_dim1 + 1], lds, (
		    ftnlen)4);
#line 416 "SB02RU.f"
	    mb02pd_("Factored", trana, n, n, &a[a_offset], lda, &s[s_offset], 
		    lds, &iwork[1], equed, &dwork[1], &dwork[1], &s[np1 * 
		    s_dim1 + 1], lds, &s[np1 + np1 * s_dim1], lds, &rcond, &
		    dwork[1], &dwork[np1], &iwork[np1], &dwork[n2 + 1], info, 
		    (ftnlen)8, (ftnlen)1, (ftnlen)1);

/*           Solve  op(A)*X = G  in  S(1:N,N+1:2*N),  using the LU */
/*           factorization of  op(A),  computed in  S(1:N,1:N),  and */
/*           iterative refinement. */

#line 425 "SB02RU.f"
	    mb02pd_("Factored", trana, n, n, &a[a_offset], lda, &s[s_offset], 
		    lds, &iwork[1], equed, &dwork[1], &dwork[1], &g[g_offset],
		     ldg, &s[np1 * s_dim1 + 1], lds, &rcond, &dwork[1], &
		    dwork[np1], &iwork[np1], &dwork[n2 + 1], info, (ftnlen)8, 
		    (ftnlen)1, (ftnlen)1);

/*                      -1 */
/*           Copy  op(A)    from  S(N+1:2*N,N+1:2*N)  in  S(1:N,1:N). */

#line 433 "SB02RU.f"
	    dlacpy_("Full", n, n, &s[np1 + np1 * s_dim1], lds, &s[s_offset], 
		    lds, (ftnlen)4);

/*                                    -1 */
/*           Compute  op(A)' + Q*op(A)  *G  in  S(N+1:2*N,N+1:2*N). */

#line 438 "SB02RU.f"
	    if (notrna) {
#line 439 "SB02RU.f"
		ma02ad_("Full", n, n, &a[a_offset], lda, &s[np1 + np1 * 
			s_dim1], lds, (ftnlen)4);
#line 440 "SB02RU.f"
	    } else {
#line 441 "SB02RU.f"
		dlacpy_("Full", n, n, &a[a_offset], lda, &s[np1 + np1 * 
			s_dim1], lds, (ftnlen)4);
#line 442 "SB02RU.f"
	    }
#line 443 "SB02RU.f"
	    dgemm_("No transpose", "No transpose", n, n, n, &c_b34, &q[
		    q_offset], ldq, &s[np1 * s_dim1 + 1], lds, &c_b34, &s[np1 
		    + np1 * s_dim1], lds, (ftnlen)12, (ftnlen)12);

#line 446 "SB02RU.f"
	} else {

/*           Complete the construction of S in (3). */

/*           Change the sign of  X. */

#line 452 "SB02RU.f"
	    i__1 = *n;
#line 452 "SB02RU.f"
	    for (j = 1; j <= i__1; ++j) {

#line 454 "SB02RU.f"
		i__2 = n2;
#line 454 "SB02RU.f"
		for (i__ = np1; i__ <= i__2; ++i__) {
#line 455 "SB02RU.f"
		    s[i__ + j * s_dim1] = -s[i__ + j * s_dim1];
#line 456 "SB02RU.f"
/* L280: */
#line 456 "SB02RU.f"
		}

#line 458 "SB02RU.f"
/* L300: */
#line 458 "SB02RU.f"
	    }

/*           Solve  op(A)'*X = I_n  in  S(N+1:2*N,N+1:2*N),  using the LU */
/*           factorization of  op(A),  computed in  S(1:N,1:N),  and */
/*           iterative refinement. */

#line 464 "SB02RU.f"
	    dlaset_("Full", n, n, &c_b33, &c_b34, &s[np1 * s_dim1 + 1], lds, (
		    ftnlen)4);
#line 465 "SB02RU.f"
	    mb02pd_("Factored", tranat, n, n, &a[a_offset], lda, &s[s_offset],
		     lds, &iwork[1], equed, &dwork[1], &dwork[1], &s[np1 * 
		    s_dim1 + 1], lds, &s[np1 + np1 * s_dim1], lds, &rcond, &
		    dwork[1], &dwork[np1], &iwork[np1], &dwork[n2 + 1], info, 
		    (ftnlen)8, (ftnlen)1, (ftnlen)1);

/*           Solve  op(A)*X' = -G  in  S(1:N,N+1:2*N),  using the LU */
/*           factorization of  op(A),  obtained in  S(1:N,1:N),  and */
/*           iterative refinement. */

#line 474 "SB02RU.f"
	    mb02pd_("Factored", trana, n, n, &a[a_offset], lda, &s[s_offset], 
		    lds, &iwork[1], equed, &dwork[1], &dwork[1], &g[g_offset],
		     ldg, &s[np1 * s_dim1 + 1], lds, &rcond, &dwork[1], &
		    dwork[np1], &iwork[np1], &dwork[n2 + 1], info, (ftnlen)8, 
		    (ftnlen)1, (ftnlen)1);

/*           Change the sign of  X  and transpose it in-situ. */

#line 481 "SB02RU.f"
	    i__1 = n2;
#line 481 "SB02RU.f"
	    for (j = np1; j <= i__1; ++j) {

#line 483 "SB02RU.f"
		i__2 = *n;
#line 483 "SB02RU.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 484 "SB02RU.f"
		    temp = -s[i__ + j * s_dim1];
#line 485 "SB02RU.f"
		    s[i__ + j * s_dim1] = -s[j - *n + (i__ + *n) * s_dim1];
#line 486 "SB02RU.f"
		    s[j - *n + (i__ + *n) * s_dim1] = temp;
#line 487 "SB02RU.f"
/* L320: */
#line 487 "SB02RU.f"
		}

#line 489 "SB02RU.f"
/* L340: */
#line 489 "SB02RU.f"
	    }
/*                                   -T */
/*           Compute  op(A) + G*op(A)  *Q  in  S(1:N,1:N). */

#line 493 "SB02RU.f"
	    if (notrna) {
#line 494 "SB02RU.f"
		dlacpy_("Full", n, n, &a[a_offset], lda, &s[s_offset], lds, (
			ftnlen)4);
#line 495 "SB02RU.f"
	    } else {
#line 496 "SB02RU.f"
		ma02ad_("Full", n, n, &a[a_offset], lda, &s[s_offset], lds, (
			ftnlen)4);
#line 497 "SB02RU.f"
	    }
#line 498 "SB02RU.f"
	    dgemm_("No transpose", "No transpose", n, n, n, &c_b57, &g[
		    g_offset], ldg, &s[np1 + s_dim1], lds, &c_b34, &s[
		    s_offset], lds, (ftnlen)12, (ftnlen)12);

#line 501 "SB02RU.f"
	}
#line 502 "SB02RU.f"
	dwork[1] = rconda;
#line 503 "SB02RU.f"
	dwork[2] = pivotg;
#line 504 "SB02RU.f"
    }
#line 505 "SB02RU.f"
    return 0;

/* *** Last line of SB02RU *** */
} /* sb02ru_ */

