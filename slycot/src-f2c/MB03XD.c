#line 1 "MB03XD.f"
/* MB03XD.f -- translated by f2c (version 20100827).
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

#line 1 "MB03XD.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b30 = 1.;
static doublereal c_b31 = -1.;
static doublereal c_b56 = 0.;
static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int mb03xd_(char *balanc, char *job, char *jobu, char *jobv, 
	integer *n, doublereal *a, integer *lda, doublereal *qg, integer *
	ldqg, doublereal *t, integer *ldt, doublereal *u1, integer *ldu1, 
	doublereal *u2, integer *ldu2, doublereal *v1, integer *ldv1, 
	doublereal *v2, integer *ldv2, doublereal *wr, doublereal *wi, 
	integer *ilo, doublereal *scale, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen balanc_len, ftnlen job_len, ftnlen jobu_len, 
	ftnlen jobv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, t_dim1, t_offset, u1_dim1, 
	    u1_offset, u2_dim1, u2_offset, v1_dim1, v1_offset, v2_dim1, 
	    v2_offset, i__1, i__2, i__3, i__4, i__5, i__6, i__7;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, l, pq, pz;
    static doublereal eps;
    static integer pdw, ilo1, ierr, pcsl;
    static doublereal hnrm, temp;
    static integer pcsr;
    extern /* Subroutine */ int ma01ad_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), mb04dd_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);
    extern doublereal ma02id_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb04qb_(char *, char *, char *, char *, char *
	    , integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen), dscal_(
	    integer *, doublereal *, doublereal *, integer *), mb04tb_(char *,
	     char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     doublereal *, integer *, integer *, ftnlen, ftnlen), dgemm_(char 
	    *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen);
    static integer pbeta;
    static logical lscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char uchar[1], vchar[1];
    extern /* Subroutine */ int mb03xp_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, ftnlen, ftnlen, ftnlen);
    static doublereal tempi;
    static logical lperm, wantg;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer ptaul;
    static doublereal tempr;
    static integer ptaur;
    static logical wants, wantu, wantv;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cscale;
    static logical scaleh;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlacpy_(char *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static doublereal bignum;
    static integer wrkmin;
    static doublereal smlnum;
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

/*     To compute the eigenvalues of a Hamiltonian matrix, */

/*                   [  A   G  ]         T        T */
/*             H  =  [       T ],   G = G,   Q = Q,                  (1) */
/*                   [  Q  -A  ] */

/*     where A, G and Q are real n-by-n matrices. */

/*     Due to the structure of H all eigenvalues appear in pairs */
/*     (lambda,-lambda). This routine computes the eigenvalues of H */
/*     using an algorithm based on the symplectic URV and the periodic */
/*     Schur decompositions as described in [1], */

/*           T       [  T   G  ] */
/*          U H V =  [       T ],                                    (2) */
/*                   [  0  -S  ] */

/*     where U and V are 2n-by-2n orthogonal symplectic matrices, */
/*     S is in real Schur form and T is upper triangular. */

/*     The algorithm is backward stable and preserves the eigenvalue */
/*     pairings in finite precision arithmetic. */

/*     Optionally, a symplectic balancing transformation to improve the */
/*     conditioning of eigenvalues is computed (see MB04DD). In this */
/*     case, the matrix H in decomposition (2) must be replaced by the */
/*     balanced matrix. */

/*     The SLICOT Library routine MB03ZD can be used to compute invariant */
/*     subspaces of H from the output of this routine. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     BALANC  CHARACTER*1 */
/*             Indicates how H should be diagonally scaled and/or */
/*             permuted to reduce its norm. */
/*             = 'N': Do not diagonally scale or permute; */
/*             = 'P': Perform symplectic permutations to make the matrix */
/*                    closer to Hamiltonian Schur form. Do not diagonally */
/*                    scale; */
/*             = 'S': Diagonally scale the matrix, i.e., replace A, G and */
/*                    Q by D*A*D**(-1), D*G*D and D**(-1)*Q*D**(-1) where */
/*                    D is a diagonal matrix chosen to make the rows and */
/*                    columns of H more equal in norm. Do not permute; */
/*             = 'B': Both diagonally scale and permute A, G and Q. */
/*             Permuting does not change the norm of H, but scaling does. */

/*     JOB     CHARACTER*1 */
/*             Indicates whether the user wishes to compute the full */
/*             decomposition (2) or the eigenvalues only, as follows: */
/*             = 'E': compute the eigenvalues only; */
/*             = 'S': compute matrices T and S of (2); */
/*             = 'G': compute matrices T, S and G of (2). */

/*     JOBU    CHARACTER*1 */
/*             Indicates whether or not the user wishes to compute the */
/*             orthogonal symplectic matrix U of (2) as follows: */
/*             = 'N': the matrix U is not computed; */
/*             = 'U': the matrix U is computed. */

/*     JOBV    CHARACTER*1 */
/*             Indicates whether or not the user wishes to compute the */
/*             orthogonal symplectic matrix V of (2) as follows: */
/*             = 'N': the matrix V is not computed; */
/*             = 'V': the matrix V is computed. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, this array is overwritten. If JOB = 'S' or */
/*             JOB = 'G', the leading N-by-N part of this array contains */
/*             the matrix S in real Schur form of decomposition (2). */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*                            (LDQG,N+1) */
/*             On entry, the leading N-by-N+1 part of this array must */
/*             contain in columns 1:N the lower triangular part of the */
/*             matrix Q and in columns 2:N+1 the upper triangular part */
/*             of the matrix G. */
/*             On exit, this array is overwritten. If JOB = 'G', the */
/*             leading N-by-N+1 part of this array contains in columns */
/*             2:N+1 the matrix G of decomposition (2). */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= max(1,N). */

/*     T       (output) DOUBLE PRECISION array, dimension (LDT,N) */
/*             On exit, if JOB = 'S' or JOB = 'G', the leading N-by-N */
/*             part of this array contains the upper triangular matrix T */
/*             of the decomposition (2). Otherwise, this array is used as */
/*             workspace. */

/*     LDT     INTEGER */
/*             The leading dimension of the array T.  LDT >= MAX(1,N). */

/*     U1      (output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On exit, if JOBU = 'U', the leading N-by-N part of this */
/*             array contains the (1,1) block of the orthogonal */
/*             symplectic matrix U of decomposition (2). */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1.  LDU1 >= 1. */
/*             LDU1 >= N,    if JOBU = 'U'. */

/*     U2      (output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On exit, if JOBU = 'U', the leading N-by-N part of this */
/*             array contains the (2,1) block of the orthogonal */
/*             symplectic matrix U of decomposition (2). */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2.  LDU2 >= 1. */
/*             LDU2 >= N,    if JOBU = 'U'. */

/*     V1      (output) DOUBLE PRECISION array, dimension (LDV1,N) */
/*             On exit, if JOBV = 'V', the leading N-by-N part of this */
/*             array contains the (1,1) block of the orthogonal */
/*             symplectic matrix V of decomposition (2). */

/*     LDV1    INTEGER */
/*             The leading dimension of the array V1.  LDV1 >= 1. */
/*             LDV1 >= N,    if JOBV = 'V'. */

/*     V2      (output) DOUBLE PRECISION array, dimension (LDV2,N) */
/*             On exit, if JOBV = 'V', the leading N-by-N part of this */
/*             array contains the (2,1) block of the orthogonal */
/*             symplectic matrix V of decomposition (2). */

/*     LDV2    INTEGER */
/*             The leading dimension of the array V2.  LDV2 >= 1. */
/*             LDV2 >= N,    if JOBV = 'V'. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             On exit, the leading N elements of WR and WI contain the */
/*             real and imaginary parts, respectively, of N eigenvalues */
/*             that have nonpositive real part. Complex conjugate pairs */
/*             of eigenvalues with real part not equal to zero will */
/*             appear consecutively with the eigenvalue having the */
/*             positive imaginary part first. For complex conjugate pairs */
/*             of eigenvalues on the imaginary axis only the eigenvalue */
/*             having nonnegative imaginary part will be returned. */

/*     ILO     (output) INTEGER */
/*             ILO is an integer value determined when H was balanced. */
/*             The balanced A(i,j) = 0 if I > J and J = 1,...,ILO-1. */
/*             The balanced Q(i,j) = 0 if J = 1,...,ILO-1 or */
/*             I = 1,...,ILO-1. */

/*     SCALE   (output) DOUBLE PRECISION array, dimension (N) */
/*             On exit, if SCALE = 'S', the leading N elements of this */
/*             array contain details of the permutation and scaling */
/*             factors applied when balancing H, see MB04DD. */
/*             This array is not referenced if BALANC = 'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -25,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  (input) INTEGER */
/*             The dimension of the array DWORK. LDWORK >= max( 1, 8*N ). */
/*             Moreover: */
/*             If JOB = 'E' or 'S' and JOBU = 'N' and JOBV = 'N', */
/*                LDWORK >= 7*N+N*N. */
/*             If JOB = 'G' and JOBU = 'N' and JOBV = 'N', */
/*                LDWORK >= max( 7*N+N*N, 2*N+3*N*N ). */
/*             If JOB = 'G' and JOBU = 'U' and JOBV = 'N', */
/*                LDWORK >= 7*N+2*N*N. */
/*             If JOB = 'G' and JOBU = 'N' and JOBV = 'V', */
/*                LDWORK >= 7*N+2*N*N. */
/*             If JOB = 'G' and JOBU = 'U' and JOBV = 'V', */
/*                LDWORK >= 7*N+N*N. */
/*             For good performance, LDWORK must generally be larger. */

/*     Error Indicator */

/*     INFO     (output) INTEGER */
/*              = 0:  successful exit; */
/*              < 0:  if INFO = -i, the i-th argument had an illegal */
/*                    value; */
/*              > 0:  if INFO = i, the periodic QR algorithm failed to */
/*                    compute all the eigenvalues, elements i+1:N of WR */
/*                    and WI contain eigenvalues which have converged. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A numerically stable, structure preserving method for */
/*         computing the eigenvalues of real Hamiltonian or symplectic */
/*         pencils. */
/*         Numer. Math., Vol. 78(3), pp. 329-358, 1998. */

/*     [2] Benner, P., Mehrmann, V., and Xu, H. */
/*         A new method for computing the stable invariant subspace of a */
/*         real Hamiltonian matrix,  J. Comput. Appl. Math., vol. 86, */
/*         pp. 17-43, 1997. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHAESU). */

/*     KEYWORDS */

/*     Eigenvalues, invariant subspace, Hamiltonian matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 284 "MB03XD.f"
    /* Parameter adjustments */
#line 284 "MB03XD.f"
    a_dim1 = *lda;
#line 284 "MB03XD.f"
    a_offset = 1 + a_dim1;
#line 284 "MB03XD.f"
    a -= a_offset;
#line 284 "MB03XD.f"
    qg_dim1 = *ldqg;
#line 284 "MB03XD.f"
    qg_offset = 1 + qg_dim1;
#line 284 "MB03XD.f"
    qg -= qg_offset;
#line 284 "MB03XD.f"
    t_dim1 = *ldt;
#line 284 "MB03XD.f"
    t_offset = 1 + t_dim1;
#line 284 "MB03XD.f"
    t -= t_offset;
#line 284 "MB03XD.f"
    u1_dim1 = *ldu1;
#line 284 "MB03XD.f"
    u1_offset = 1 + u1_dim1;
#line 284 "MB03XD.f"
    u1 -= u1_offset;
#line 284 "MB03XD.f"
    u2_dim1 = *ldu2;
#line 284 "MB03XD.f"
    u2_offset = 1 + u2_dim1;
#line 284 "MB03XD.f"
    u2 -= u2_offset;
#line 284 "MB03XD.f"
    v1_dim1 = *ldv1;
#line 284 "MB03XD.f"
    v1_offset = 1 + v1_dim1;
#line 284 "MB03XD.f"
    v1 -= v1_offset;
#line 284 "MB03XD.f"
    v2_dim1 = *ldv2;
#line 284 "MB03XD.f"
    v2_offset = 1 + v2_dim1;
#line 284 "MB03XD.f"
    v2 -= v2_offset;
#line 284 "MB03XD.f"
    --wr;
#line 284 "MB03XD.f"
    --wi;
#line 284 "MB03XD.f"
    --scale;
#line 284 "MB03XD.f"
    --dwork;
#line 284 "MB03XD.f"

#line 284 "MB03XD.f"
    /* Function Body */
#line 284 "MB03XD.f"
    *info = 0;
#line 285 "MB03XD.f"
    lperm = lsame_(balanc, "P", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "B", (
	    ftnlen)1, (ftnlen)1);
#line 286 "MB03XD.f"
    lscal = lsame_(balanc, "S", (ftnlen)1, (ftnlen)1) || lsame_(balanc, "B", (
	    ftnlen)1, (ftnlen)1);
#line 287 "MB03XD.f"
    wants = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "G", (
	    ftnlen)1, (ftnlen)1);
#line 288 "MB03XD.f"
    wantg = lsame_(job, "G", (ftnlen)1, (ftnlen)1);
#line 289 "MB03XD.f"
    wantu = lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);
#line 290 "MB03XD.f"
    wantv = lsame_(jobv, "V", (ftnlen)1, (ftnlen)1);

#line 292 "MB03XD.f"
    if (wantg) {
#line 293 "MB03XD.f"
	if (wantu) {
#line 294 "MB03XD.f"
	    if (wantv) {
/* Computing MAX */
#line 295 "MB03XD.f"
		i__1 = 1, i__2 = *n * 7 + *n * *n;
#line 295 "MB03XD.f"
		wrkmin = max(i__1,i__2);
#line 296 "MB03XD.f"
	    } else {
/* Computing MAX */
#line 297 "MB03XD.f"
		i__1 = 1, i__2 = *n * 7 + (*n << 1) * *n;
#line 297 "MB03XD.f"
		wrkmin = max(i__1,i__2);
#line 298 "MB03XD.f"
	    }
#line 299 "MB03XD.f"
	} else {
#line 300 "MB03XD.f"
	    if (wantv) {
/* Computing MAX */
#line 301 "MB03XD.f"
		i__1 = 1, i__2 = *n * 7 + (*n << 1) * *n;
#line 301 "MB03XD.f"
		wrkmin = max(i__1,i__2);
#line 302 "MB03XD.f"
	    } else {
/* Computing MAX */
#line 303 "MB03XD.f"
		i__1 = 1, i__2 = *n * 7 + *n * *n, i__1 = max(i__1,i__2), 
			i__2 = (*n << 1) + *n * 3 * *n;
#line 303 "MB03XD.f"
		wrkmin = max(i__1,i__2);
#line 304 "MB03XD.f"
	    }
#line 305 "MB03XD.f"
	}
#line 306 "MB03XD.f"
    } else {
#line 307 "MB03XD.f"
	if (wantu) {
#line 308 "MB03XD.f"
	    if (wantv) {
/* Computing MAX */
#line 309 "MB03XD.f"
		i__1 = 1, i__2 = *n << 3;
#line 309 "MB03XD.f"
		wrkmin = max(i__1,i__2);
#line 310 "MB03XD.f"
	    } else {
/* Computing MAX */
#line 311 "MB03XD.f"
		i__1 = 1, i__2 = *n << 3;
#line 311 "MB03XD.f"
		wrkmin = max(i__1,i__2);
#line 312 "MB03XD.f"
	    }
#line 313 "MB03XD.f"
	} else {
#line 314 "MB03XD.f"
	    if (wantv) {
/* Computing MAX */
#line 315 "MB03XD.f"
		i__1 = 1, i__2 = *n << 3;
#line 315 "MB03XD.f"
		wrkmin = max(i__1,i__2);
#line 316 "MB03XD.f"
	    } else {
/* Computing MAX */
#line 317 "MB03XD.f"
		i__1 = 1, i__2 = *n * 7 + *n * *n;
#line 317 "MB03XD.f"
		wrkmin = max(i__1,i__2);
#line 318 "MB03XD.f"
	    }
#line 319 "MB03XD.f"
	}
#line 320 "MB03XD.f"
    }

#line 322 "MB03XD.f"
    wrkopt = wrkmin;

/*     Test the scalar input parameters. */

#line 326 "MB03XD.f"
    if (! lperm && ! lscal && ! lsame_(balanc, "N", (ftnlen)1, (ftnlen)1)) {
#line 328 "MB03XD.f"
	*info = -1;
#line 329 "MB03XD.f"
    } else if (! wants && ! lsame_(job, "E", (ftnlen)1, (ftnlen)1)) {
#line 330 "MB03XD.f"
	*info = -2;
#line 331 "MB03XD.f"
    } else if (! wantu && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
#line 332 "MB03XD.f"
	*info = -3;
#line 333 "MB03XD.f"
    } else if (! wantv && ! lsame_(jobv, "N", (ftnlen)1, (ftnlen)1)) {
#line 334 "MB03XD.f"
	*info = -4;
#line 335 "MB03XD.f"
    } else if (*n < 0) {
#line 336 "MB03XD.f"
	*info = -5;
#line 337 "MB03XD.f"
    } else if (*lda < max(1,*n)) {
#line 338 "MB03XD.f"
	*info = -7;
#line 339 "MB03XD.f"
    } else if (*ldqg < max(1,*n)) {
#line 340 "MB03XD.f"
	*info = -9;
#line 341 "MB03XD.f"
    } else if (*ldt < max(1,*n)) {
#line 342 "MB03XD.f"
	*info = -11;
#line 343 "MB03XD.f"
    } else if (*ldu1 < 1 || wantu && *ldu1 < *n) {
#line 344 "MB03XD.f"
	*info = -13;
#line 345 "MB03XD.f"
    } else if (*ldu2 < 1 || wantu && *ldu2 < *n) {
#line 346 "MB03XD.f"
	*info = -15;
#line 347 "MB03XD.f"
    } else if (*ldv1 < 1 || wantv && *ldv1 < *n) {
#line 348 "MB03XD.f"
	*info = -17;
#line 349 "MB03XD.f"
    } else if (*ldv2 < 1 || wantv && *ldv2 < *n) {
#line 350 "MB03XD.f"
	*info = -19;
#line 351 "MB03XD.f"
    } else if (*ldwork < wrkmin) {
#line 352 "MB03XD.f"
	*info = -25;
#line 353 "MB03XD.f"
	dwork[1] = (doublereal) wrkmin;
#line 354 "MB03XD.f"
    }

/*     Return if there were illegal values. */

#line 358 "MB03XD.f"
    if (*info != 0) {
#line 359 "MB03XD.f"
	i__1 = -(*info);
#line 359 "MB03XD.f"
	xerbla_("MB03XD", &i__1, (ftnlen)6);
#line 360 "MB03XD.f"
	return 0;
#line 361 "MB03XD.f"
    }

/*     Quick return if possible. */

#line 365 "MB03XD.f"
    *ilo = 0;
#line 366 "MB03XD.f"
    if (*n == 0) {
#line 366 "MB03XD.f"
	return 0;
#line 366 "MB03XD.f"
    }

#line 369 "MB03XD.f"
    eps = dlamch_("P", (ftnlen)1);
#line 370 "MB03XD.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 371 "MB03XD.f"
    bignum = 1. / smlnum;
#line 372 "MB03XD.f"
    dlabad_(&smlnum, &bignum);
#line 373 "MB03XD.f"
    smlnum = sqrt(smlnum) / eps;
#line 374 "MB03XD.f"
    bignum = 1. / smlnum;

/*     Scale H if maximal element is outside range [SMLNUM,BIGNUM]. */

#line 378 "MB03XD.f"
    hnrm = ma02id_("Hamiltonian", "MaxElement", n, &a[a_offset], lda, &qg[
	    qg_offset], ldqg, &dwork[1], (ftnlen)11, (ftnlen)10);
#line 380 "MB03XD.f"
    scaleh = FALSE_;
#line 381 "MB03XD.f"
    if (hnrm > 0. && hnrm < smlnum) {
#line 382 "MB03XD.f"
	scaleh = TRUE_;
#line 383 "MB03XD.f"
	cscale = smlnum;
#line 384 "MB03XD.f"
    } else if (hnrm > bignum) {
#line 385 "MB03XD.f"
	scaleh = TRUE_;
#line 386 "MB03XD.f"
	cscale = bignum;
#line 387 "MB03XD.f"
    }
#line 388 "MB03XD.f"
    if (scaleh) {
#line 389 "MB03XD.f"
	dlascl_("General", &c__0, &c__0, &hnrm, &cscale, n, n, &a[a_offset], 
		lda, &ierr, (ftnlen)7);
#line 390 "MB03XD.f"
	i__1 = *n + 1;
#line 390 "MB03XD.f"
	dlascl_("General", &c__0, &c__0, &hnrm, &cscale, n, &i__1, &qg[
		qg_offset], ldqg, &ierr, (ftnlen)7);
#line 392 "MB03XD.f"
    }

/*     Balance the matrix. */

#line 396 "MB03XD.f"
    mb04dd_(balanc, n, &a[a_offset], lda, &qg[qg_offset], ldqg, ilo, &scale[1]
	    , &ierr, (ftnlen)1);

/*     Copy A to T and multiply A by -1. */

#line 400 "MB03XD.f"
    dlacpy_("All", n, n, &a[a_offset], lda, &t[t_offset], ldt, (ftnlen)3);
#line 401 "MB03XD.f"
    dlascl_("General", &c__0, &c__0, &c_b30, &c_b31, n, n, &a[a_offset], lda, 
	    &ierr, (ftnlen)7);

/*     --------------------------------------------- */
/*     Step 1: Compute symplectic URV decomposition. */
/*     --------------------------------------------- */

#line 407 "MB03XD.f"
    pcsl = 1;
#line 408 "MB03XD.f"
    pcsr = pcsl + (*n << 1);
#line 409 "MB03XD.f"
    ptaul = pcsr + (*n << 1);
#line 410 "MB03XD.f"
    ptaur = ptaul + *n;
#line 411 "MB03XD.f"
    pdw = ptaur + *n;
#line 413 "MB03XD.f"
    if (! wantu && ! wantv) {

/*         Copy Q and Q' to workspace. */

#line 417 "MB03XD.f"
	pq = pdw;
#line 418 "MB03XD.f"
	pdw += *n * *n;
#line 419 "MB03XD.f"
	i__1 = *n;
#line 419 "MB03XD.f"
	for (j = 1; j <= i__1; ++j) {
#line 420 "MB03XD.f"
	    k = pq + (*n + 1) * (j - 1);
#line 421 "MB03XD.f"
	    l = k;
#line 422 "MB03XD.f"
	    dwork[k] = qg[j + j * qg_dim1];
#line 423 "MB03XD.f"
	    i__2 = *n;
#line 423 "MB03XD.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 424 "MB03XD.f"
		++k;
#line 425 "MB03XD.f"
		l += *n;
#line 426 "MB03XD.f"
		temp = qg[i__ + j * qg_dim1];
#line 427 "MB03XD.f"
		dwork[k] = temp;
#line 428 "MB03XD.f"
		dwork[l] = temp;
#line 429 "MB03XD.f"
/* L10: */
#line 429 "MB03XD.f"
	    }
#line 430 "MB03XD.f"
/* L20: */
#line 430 "MB03XD.f"
	}
#line 431 "MB03XD.f"
    } else if (wantu) {

/*         Copy Q and Q' to U2. */

#line 435 "MB03XD.f"
	i__1 = *n;
#line 435 "MB03XD.f"
	for (j = 1; j <= i__1; ++j) {
#line 436 "MB03XD.f"
	    u2[j + j * u2_dim1] = qg[j + j * qg_dim1];
#line 437 "MB03XD.f"
	    i__2 = *n;
#line 437 "MB03XD.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 438 "MB03XD.f"
		temp = qg[i__ + j * qg_dim1];
#line 439 "MB03XD.f"
		u2[i__ + j * u2_dim1] = temp;
#line 440 "MB03XD.f"
		u2[j + i__ * u2_dim1] = temp;
#line 441 "MB03XD.f"
/* L30: */
#line 441 "MB03XD.f"
	    }
#line 442 "MB03XD.f"
/* L40: */
#line 442 "MB03XD.f"
	}
#line 443 "MB03XD.f"
    } else {

/*         Copy Q and Q' to V2. */

#line 447 "MB03XD.f"
	i__1 = *n;
#line 447 "MB03XD.f"
	for (j = 1; j <= i__1; ++j) {
#line 448 "MB03XD.f"
	    v2[j + j * v2_dim1] = qg[j + j * qg_dim1];
#line 449 "MB03XD.f"
	    i__2 = *n;
#line 449 "MB03XD.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 450 "MB03XD.f"
		temp = qg[i__ + j * qg_dim1];
#line 451 "MB03XD.f"
		v2[i__ + j * v2_dim1] = temp;
#line 452 "MB03XD.f"
		v2[j + i__ * v2_dim1] = temp;
#line 453 "MB03XD.f"
/* L50: */
#line 453 "MB03XD.f"
	    }
#line 454 "MB03XD.f"
/* L60: */
#line 454 "MB03XD.f"
	}
#line 455 "MB03XD.f"
    }

/*     Transpose G. */

#line 459 "MB03XD.f"
    i__1 = *n;
#line 459 "MB03XD.f"
    for (j = 1; j <= i__1; ++j) {
#line 460 "MB03XD.f"
	i__2 = *n;
#line 460 "MB03XD.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 461 "MB03XD.f"
	    qg[i__ + (j + 1) * qg_dim1] = qg[j + (i__ + 1) * qg_dim1];
#line 462 "MB03XD.f"
/* L70: */
#line 462 "MB03XD.f"
	}
#line 463 "MB03XD.f"
/* L80: */
#line 463 "MB03XD.f"
    }

#line 465 "MB03XD.f"
    if (! wantu && ! wantv) {
#line 466 "MB03XD.f"
	i__1 = *ldwork - pdw + 1;
#line 466 "MB03XD.f"
	mb04tb_("Not Transposed", "Transposed", n, ilo, &t[t_offset], ldt, &a[
		a_offset], lda, &qg[(qg_dim1 << 1) + 1], ldqg, &dwork[pq], n, 
		&dwork[pcsl], &dwork[pcsr], &dwork[ptaul], &dwork[ptaur], &
		dwork[pdw], &i__1, &ierr, (ftnlen)14, (ftnlen)10);
#line 470 "MB03XD.f"
    } else if (wantu) {
#line 471 "MB03XD.f"
	i__1 = *ldwork - pdw + 1;
#line 471 "MB03XD.f"
	mb04tb_("Not Transposed", "Transposed", n, ilo, &t[t_offset], ldt, &a[
		a_offset], lda, &qg[(qg_dim1 << 1) + 1], ldqg, &u2[u2_offset],
		 ldu2, &dwork[pcsl], &dwork[pcsr], &dwork[ptaul], &dwork[
		ptaur], &dwork[pdw], &i__1, &ierr, (ftnlen)14, (ftnlen)10);
#line 475 "MB03XD.f"
    } else {
#line 476 "MB03XD.f"
	i__1 = *ldwork - pdw + 1;
#line 476 "MB03XD.f"
	mb04tb_("Not Transposed", "Transposed", n, ilo, &t[t_offset], ldt, &a[
		a_offset], lda, &qg[(qg_dim1 << 1) + 1], ldqg, &v2[v2_offset],
		 ldv2, &dwork[pcsl], &dwork[pcsr], &dwork[ptaul], &dwork[
		ptaur], &dwork[pdw], &i__1, &ierr, (ftnlen)14, (ftnlen)10);
#line 480 "MB03XD.f"
    }
/* Computing MAX */
#line 481 "MB03XD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 481 "MB03XD.f"
    wrkopt = max(i__1,i__2);

#line 483 "MB03XD.f"
    if (wantu && ! wantv && ! wantg) {
#line 484 "MB03XD.f"
	if (*n > 1) {
#line 484 "MB03XD.f"
	    i__1 = *n - 1;
#line 484 "MB03XD.f"
	    i__2 = *n - 1;
#line 484 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &t[t_dim1 + 2], ldt, &qg[qg_dim1 + 
		    2], ldqg, (ftnlen)5);
#line 484 "MB03XD.f"
	}
#line 486 "MB03XD.f"
    } else if (! wantu && wantv && ! wantg) {
#line 487 "MB03XD.f"
	if (*n > 1) {
#line 488 "MB03XD.f"
	    i__1 = *n - 1;
#line 488 "MB03XD.f"
	    i__2 = *n - 1;
#line 488 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &qg[qg_dim1 + 
		    2], ldqg, (ftnlen)5);
#line 489 "MB03XD.f"
	    i__1 = *n - 1;
#line 489 "MB03XD.f"
	    i__2 = *n - 1;
#line 489 "MB03XD.f"
	    dlacpy_("Upper", &i__1, &i__2, &v2[(v2_dim1 << 1) + 1], ldv2, &qg[
		    (qg_dim1 << 1) + 1], ldqg, (ftnlen)5);
#line 491 "MB03XD.f"
	}
#line 492 "MB03XD.f"
    } else if (wantu && wantv && ! wantg) {
#line 493 "MB03XD.f"
	if (*n > 1) {
#line 494 "MB03XD.f"
	    i__1 = *n - 1;
#line 494 "MB03XD.f"
	    i__2 = *n - 1;
#line 494 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &t[t_dim1 + 2], ldt, &v2[v2_dim1 + 
		    2], ldv2, (ftnlen)5);
#line 495 "MB03XD.f"
	    i__1 = *n - 1;
#line 495 "MB03XD.f"
	    i__2 = *n - 1;
#line 495 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 2], lda, &qg[qg_dim1 + 
		    2], ldqg, (ftnlen)5);
#line 496 "MB03XD.f"
	}
#line 497 "MB03XD.f"
    } else if (wantu && ! wantv && wantg) {
#line 498 "MB03XD.f"
	if (*n > 1) {
#line 498 "MB03XD.f"
	    i__1 = *n - 1;
#line 498 "MB03XD.f"
	    i__2 = *n - 1;
#line 498 "MB03XD.f"
	    i__3 = *n - 1;
#line 498 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &t[t_dim1 + 2], ldt, &dwork[pdw + *
		    n * *n + *n], &i__3, (ftnlen)5);
#line 498 "MB03XD.f"
	}
#line 501 "MB03XD.f"
    } else if (! wantu && wantv && wantg) {
#line 502 "MB03XD.f"
	if (*n > 2) {
#line 502 "MB03XD.f"
	    i__1 = *n - 2;
#line 502 "MB03XD.f"
	    i__2 = *n - 2;
#line 502 "MB03XD.f"
	    i__3 = *n - 2;
#line 502 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 3], lda, &dwork[pdw + *
		    n * *n + *n], &i__3, (ftnlen)5);
#line 502 "MB03XD.f"
	}
#line 505 "MB03XD.f"
    } else if (wantu && wantv && wantg) {
#line 506 "MB03XD.f"
	if (*n > 1) {
#line 506 "MB03XD.f"
	    i__1 = *n - 1;
#line 506 "MB03XD.f"
	    i__2 = *n - 1;
#line 506 "MB03XD.f"
	    i__3 = *n - 1;
#line 506 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &t[t_dim1 + 2], ldt, &dwork[pdw + *
		    n], &i__3, (ftnlen)5);
#line 506 "MB03XD.f"
	}
#line 509 "MB03XD.f"
	if (*n > 2) {
#line 509 "MB03XD.f"
	    i__1 = *n - 2;
#line 509 "MB03XD.f"
	    i__2 = *n - 2;
#line 509 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &a[a_dim1 + 3], lda, &v2[v2_dim1 + 
		    3], ldv2, (ftnlen)5);
#line 509 "MB03XD.f"
	}
#line 511 "MB03XD.f"
    }

/*     ---------------------------------------------- */
/*     Step 2:  Compute periodic Schur decomposition. */
/*     ---------------------------------------------- */

#line 517 "MB03XD.f"
    if (*n > 2) {
#line 517 "MB03XD.f"
	i__1 = *n - 2;
#line 517 "MB03XD.f"
	i__2 = *n - 2;
#line 517 "MB03XD.f"
	dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &a[a_dim1 + 3], lda, (
		ftnlen)5);
#line 517 "MB03XD.f"
    }
#line 519 "MB03XD.f"
    if (*n > 1) {
#line 519 "MB03XD.f"
	i__1 = *n - 1;
#line 519 "MB03XD.f"
	i__2 = *n - 1;
#line 519 "MB03XD.f"
	dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], ldt, (
		ftnlen)5);
#line 519 "MB03XD.f"
    }
#line 521 "MB03XD.f"
    if (! wantu && ! wantv) {
#line 522 "MB03XD.f"
	pbeta = 1;
#line 523 "MB03XD.f"
    } else {
#line 524 "MB03XD.f"
	pbeta = pdw;
#line 525 "MB03XD.f"
    }

#line 527 "MB03XD.f"
    if (! wantg) {

/*        Workspace requirements: 2*N (8*N with U or V). */

#line 531 "MB03XD.f"
	pdw = pbeta + *n;
#line 532 "MB03XD.f"
	if (wantu) {
#line 533 "MB03XD.f"
	    *(unsigned char *)uchar = 'I';
#line 534 "MB03XD.f"
	} else {
#line 535 "MB03XD.f"
	    *(unsigned char *)uchar = 'N';
#line 536 "MB03XD.f"
	}
#line 537 "MB03XD.f"
	if (wantv) {
#line 538 "MB03XD.f"
	    *(unsigned char *)vchar = 'I';
#line 539 "MB03XD.f"
	} else {
#line 540 "MB03XD.f"
	    *(unsigned char *)vchar = 'N';
#line 541 "MB03XD.f"
	}
#line 542 "MB03XD.f"
	i__1 = *ldwork - pdw + 1;
#line 542 "MB03XD.f"
	mb03xp_(job, vchar, uchar, n, ilo, n, &a[a_offset], lda, &t[t_offset],
		 ldt, &v1[v1_offset], ldv1, &u1[u1_offset], ldu1, &wr[1], &wi[
		1], &dwork[pbeta], &dwork[pdw], &i__1, info, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
#line 545 "MB03XD.f"
	if (*info != 0) {
#line 545 "MB03XD.f"
	    goto L90;
#line 545 "MB03XD.f"
	}
/* Computing MAX */
#line 547 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 547 "MB03XD.f"
	wrkopt = max(i__1,i__2);

#line 549 "MB03XD.f"
    } else if (! wantu && ! wantv && wantg) {

/*        Workspace requirements: 3*N*N + 2*N. */

#line 553 "MB03XD.f"
	pq = pbeta + *n;
#line 554 "MB03XD.f"
	pz = pq + *n * *n;
#line 555 "MB03XD.f"
	pdw = pz + *n * *n;
#line 556 "MB03XD.f"
	i__1 = *ldwork - pdw + 1;
#line 556 "MB03XD.f"
	mb03xp_("Schur", "Init", "Init", n, ilo, n, &a[a_offset], lda, &t[
		t_offset], ldt, &dwork[pq], n, &dwork[pz], n, &wr[1], &wi[1], 
		&dwork[pbeta], &dwork[pdw], &i__1, info, (ftnlen)5, (ftnlen)4,
		 (ftnlen)4);
#line 559 "MB03XD.f"
	if (*info != 0) {
#line 559 "MB03XD.f"
	    goto L90;
#line 559 "MB03XD.f"
	}
/* Computing MAX */
#line 561 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 561 "MB03XD.f"
	wrkopt = max(i__1,i__2);
#line 562 "MB03XD.f"
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pz], n, &
		qg[(qg_dim1 << 1) + 1], ldqg, &c_b56, &dwork[pdw], n, (ftnlen)
		9, (ftnlen)12);
#line 564 "MB03XD.f"
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pdw], 
		n, &dwork[pq], n, &c_b56, &qg[(qg_dim1 << 1) + 1], ldqg, (
		ftnlen)12, (ftnlen)12);
#line 566 "MB03XD.f"
    } else if (wantu && ! wantv && wantg) {

/*        Workspace requirements: 2*N*N + 7*N. */

#line 570 "MB03XD.f"
	pq = pbeta + *n;
#line 571 "MB03XD.f"
	pdw = pq + *n * *n;
#line 572 "MB03XD.f"
	i__1 = *ldwork - pdw - (*n - 1) * (*n - 1) + 1;
#line 572 "MB03XD.f"
	mb03xp_("Schur", "Init", "Init", n, ilo, n, &a[a_offset], lda, &t[
		t_offset], ldt, &dwork[pq], n, &u1[u1_offset], ldu1, &wr[1], &
		wi[1], &dwork[pbeta], &dwork[pdw + (*n - 1) * (*n - 1)], &
		i__1, info, (ftnlen)5, (ftnlen)4, (ftnlen)4);
#line 576 "MB03XD.f"
	if (*info != 0) {
#line 576 "MB03XD.f"
	    goto L90;
#line 576 "MB03XD.f"
	}
/* Computing MAX */
#line 578 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + (*n - 1) * (*n - 1)] + 
		pdw + (*n - 1) * (*n - 1) - 1;
#line 578 "MB03XD.f"
	wrkopt = max(i__1,i__2);
#line 580 "MB03XD.f"
	if (*n > 1) {
#line 580 "MB03XD.f"
	    i__1 = *n - 1;
#line 580 "MB03XD.f"
	    i__2 = *n - 1;
#line 580 "MB03XD.f"
	    i__3 = *n - 1;
#line 580 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &dwork[pdw], &i__3, &t[t_dim1 + 2],
		     ldt, (ftnlen)5);
#line 580 "MB03XD.f"
	}
#line 583 "MB03XD.f"
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b30, &u1[u1_offset], 
		ldu1, &qg[(qg_dim1 << 1) + 1], ldqg, &c_b56, &dwork[pdw], n, (
		ftnlen)9, (ftnlen)12);
#line 585 "MB03XD.f"
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pdw], 
		n, &dwork[pq], n, &c_b56, &qg[(qg_dim1 << 1) + 1], ldqg, (
		ftnlen)12, (ftnlen)12);

#line 588 "MB03XD.f"
    } else if (! wantu && wantv && wantg) {

/*        Workspace requirements: 2*N*N + 7*N */

#line 592 "MB03XD.f"
	pz = pbeta + *n;
#line 593 "MB03XD.f"
	pdw = pz + *n * *n;
#line 594 "MB03XD.f"
	i__1 = *ldwork - pdw - (*n - 1) * (*n - 1) + 1;
#line 594 "MB03XD.f"
	mb03xp_("Schur", "Init", "Init", n, ilo, n, &a[a_offset], lda, &t[
		t_offset], ldt, &v1[v1_offset], ldv1, &dwork[pz], n, &wr[1], &
		wi[1], &dwork[pbeta], &dwork[pdw + (*n - 1) * (*n - 1)], &
		i__1, info, (ftnlen)5, (ftnlen)4, (ftnlen)4);
#line 598 "MB03XD.f"
	if (*info != 0) {
#line 598 "MB03XD.f"
	    goto L90;
#line 598 "MB03XD.f"
	}
/* Computing MAX */
#line 600 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + (*n - 1) * (*n - 1)] + 
		pdw + (*n - 1) * (*n - 1) - 1;
#line 600 "MB03XD.f"
	wrkopt = max(i__1,i__2);
#line 602 "MB03XD.f"
	if (*n > 2) {
#line 602 "MB03XD.f"
	    i__1 = *n - 2;
#line 602 "MB03XD.f"
	    i__2 = *n - 2;
#line 602 "MB03XD.f"
	    i__3 = *n - 2;
#line 602 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &dwork[pdw], &i__3, &a[a_dim1 + 3],
		     lda, (ftnlen)5);
#line 602 "MB03XD.f"
	}
#line 605 "MB03XD.f"
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pz], n, &
		qg[(qg_dim1 << 1) + 1], ldqg, &c_b56, &dwork[pdw], n, (ftnlen)
		9, (ftnlen)12);
#line 607 "MB03XD.f"
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pdw], 
		n, &v1[v1_offset], ldv1, &c_b56, &qg[(qg_dim1 << 1) + 1], 
		ldqg, (ftnlen)12, (ftnlen)12);

#line 610 "MB03XD.f"
    } else if (wantu && wantv && wantg) {

/*        Workspace requirements: N*N + 7*N. */

#line 614 "MB03XD.f"
	pdw = pbeta + *n;
#line 615 "MB03XD.f"
	i__1 = *ldwork - pdw - (*n - 1) * (*n - 1) + 1;
#line 615 "MB03XD.f"
	mb03xp_("Schur", "Init", "Init", n, ilo, n, &a[a_offset], lda, &t[
		t_offset], ldt, &v1[v1_offset], ldv1, &u1[u1_offset], ldu1, &
		wr[1], &wi[1], &dwork[pbeta], &dwork[pdw + (*n - 1) * (*n - 1)
		], &i__1, info, (ftnlen)5, (ftnlen)4, (ftnlen)4);
#line 619 "MB03XD.f"
	if (*info != 0) {
#line 619 "MB03XD.f"
	    goto L90;
#line 619 "MB03XD.f"
	}
/* Computing MAX */
#line 621 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + (*n - 1) * (*n - 1)] + 
		pdw + (*n - 1) * (*n - 1) - 1;
#line 621 "MB03XD.f"
	wrkopt = max(i__1,i__2);
#line 623 "MB03XD.f"
	if (*n > 1) {
#line 623 "MB03XD.f"
	    i__1 = *n - 1;
#line 623 "MB03XD.f"
	    i__2 = *n - 1;
#line 623 "MB03XD.f"
	    i__3 = *n - 1;
#line 623 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &dwork[pdw], &i__3, &t[t_dim1 + 2],
		     ldt, (ftnlen)5);
#line 623 "MB03XD.f"
	}
#line 626 "MB03XD.f"
	if (*n > 2) {
#line 626 "MB03XD.f"
	    i__1 = *n - 2;
#line 626 "MB03XD.f"
	    i__2 = *n - 2;
#line 626 "MB03XD.f"
	    dlacpy_("Lower", &i__1, &i__2, &v2[v2_dim1 + 3], ldv2, &a[a_dim1 
		    + 3], lda, (ftnlen)5);
#line 626 "MB03XD.f"
	}
#line 628 "MB03XD.f"
	dgemm_("Transpose", "No Transpose", n, n, n, &c_b30, &u1[u1_offset], 
		ldu1, &qg[(qg_dim1 << 1) + 1], ldqg, &c_b56, &dwork[pdw], n, (
		ftnlen)9, (ftnlen)12);
#line 630 "MB03XD.f"
	dgemm_("No Transpose", "No Transpose", n, n, n, &c_b30, &dwork[pdw], 
		n, &v1[v1_offset], ldv1, &c_b56, &qg[(qg_dim1 << 1) + 1], 
		ldqg, (ftnlen)12, (ftnlen)12);
#line 632 "MB03XD.f"
    }

#line 634 "MB03XD.f"
L90:

/*     Compute square roots of eigenvalues and rescale. */

#line 638 "MB03XD.f"
    i__1 = *n;
#line 638 "MB03XD.f"
    for (i__ = *info + 1; i__ <= i__1; ++i__) {
#line 639 "MB03XD.f"
	tempr = wr[i__];
#line 640 "MB03XD.f"
	tempi = wi[i__];
#line 641 "MB03XD.f"
	temp = dwork[pbeta + i__ - 1];
#line 642 "MB03XD.f"
	if (temp > 0.) {
#line 642 "MB03XD.f"
	    tempr = -tempr;
#line 642 "MB03XD.f"
	}
#line 644 "MB03XD.f"
	temp = abs(temp);
#line 645 "MB03XD.f"
	if (tempi == 0.) {
#line 646 "MB03XD.f"
	    if (tempr < 0.) {
#line 647 "MB03XD.f"
		wr[i__] = 0.;
#line 648 "MB03XD.f"
		wi[i__] = sqrt(temp) * sqrt(-tempr);
#line 649 "MB03XD.f"
	    } else {
#line 650 "MB03XD.f"
		wr[i__] = -sqrt(temp) * sqrt(tempr);
#line 651 "MB03XD.f"
		wi[i__] = 0.;
#line 652 "MB03XD.f"
	    }
#line 653 "MB03XD.f"
	} else {
#line 654 "MB03XD.f"
	    ma01ad_(&tempr, &tempi, &wr[i__], &wi[i__]);
#line 655 "MB03XD.f"
	    wr[i__] = -wr[i__] * sqrt(temp);
#line 656 "MB03XD.f"
	    if (temp > 0.) {
#line 657 "MB03XD.f"
		wi[i__] *= sqrt(temp);
#line 658 "MB03XD.f"
	    } else {
#line 659 "MB03XD.f"
		wi[i__] = 0.;
#line 660 "MB03XD.f"
	    }
#line 661 "MB03XD.f"
	}
#line 662 "MB03XD.f"
/* L100: */
#line 662 "MB03XD.f"
    }

#line 664 "MB03XD.f"
    if (scaleh) {

/*        Undo scaling. */

#line 668 "MB03XD.f"
	dlascl_("Hessenberg", &c__0, &c__0, &cscale, &hnrm, n, n, &a[a_offset]
		, lda, &ierr, (ftnlen)10);
#line 670 "MB03XD.f"
	dlascl_("Upper", &c__0, &c__0, &cscale, &hnrm, n, n, &t[t_offset], 
		ldt, &ierr, (ftnlen)5);
#line 671 "MB03XD.f"
	if (wantg) {
#line 671 "MB03XD.f"
	    dlascl_("General", &c__0, &c__0, &cscale, &hnrm, n, n, &qg[(
		    qg_dim1 << 1) + 1], ldqg, &ierr, (ftnlen)7);
#line 671 "MB03XD.f"
	}
#line 674 "MB03XD.f"
	dlascl_("General", &c__0, &c__0, &cscale, &hnrm, n, &c__1, &wr[1], n, 
		&ierr, (ftnlen)7);
#line 675 "MB03XD.f"
	dlascl_("General", &c__0, &c__0, &cscale, &hnrm, n, &c__1, &wi[1], n, 
		&ierr, (ftnlen)7);
#line 676 "MB03XD.f"
    }

#line 678 "MB03XD.f"
    if (*info != 0) {
#line 678 "MB03XD.f"
	return 0;
#line 678 "MB03XD.f"
    }

/*     ----------------------------------------------- */
/*     Step 3:  Compute orthogonal symplectic factors. */
/*     ----------------------------------------------- */

/*     Fix CSL and CSR for MB04QB. */

#line 687 "MB03XD.f"
    if (wantu) {
#line 687 "MB03XD.f"
	dscal_(n, &c_b31, &dwork[pcsl + 1], &c__2);
#line 687 "MB03XD.f"
    }
#line 689 "MB03XD.f"
    if (wantv) {
#line 689 "MB03XD.f"
	i__1 = *n - 1;
#line 689 "MB03XD.f"
	dscal_(&i__1, &c_b31, &dwork[pcsr + 1], &c__2);
#line 689 "MB03XD.f"
    }
/* Computing MIN */
#line 691 "MB03XD.f"
    i__1 = *n, i__2 = *ilo + 1;
#line 691 "MB03XD.f"
    ilo1 = min(i__1,i__2);

#line 693 "MB03XD.f"
    if (wantu && ! wantv && ! wantg) {

/*        Workspace requirements: 7*N. */

#line 697 "MB03XD.f"
	pdw = ptaur;
#line 698 "MB03XD.f"
	i__1 = *ldt + 1;
#line 698 "MB03XD.f"
	dcopy_(n, &t[t_dim1 + 1], &i__1, &dwork[pdw], &c__1);
#line 699 "MB03XD.f"
	dlacpy_("Lower", n, n, &u2[u2_offset], ldu2, &t[t_offset], ldt, (
		ftnlen)5);
#line 700 "MB03XD.f"
	dlaset_("All", n, n, &c_b56, &c_b56, &u2[u2_offset], ldu2, (ftnlen)3);
#line 701 "MB03XD.f"
	i__1 = *n - *ilo + 1;
#line 701 "MB03XD.f"
	i__2 = *n - *ilo + 1;
#line 701 "MB03XD.f"
	i__3 = *ldwork - pdw - *n + 1;
#line 701 "MB03XD.f"
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Columnwise", &i__1, n, &i__2, &qg[*ilo + *ilo * qg_dim1], 
		ldqg, &t[*ilo + *ilo * t_dim1], ldt, &u1[*ilo + u1_dim1], 
		ldu1, &u2[*ilo + u2_dim1], ldu2, &dwork[pcsl + (*ilo << 1) - 
		2], &dwork[ptaul + *ilo - 1], &dwork[pdw + *n], &i__3, &ierr, 
		(ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)10);
/* Computing MAX */
#line 707 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *n] + pdw + *n - 1;
#line 707 "MB03XD.f"
	wrkopt = max(i__1,i__2);
#line 708 "MB03XD.f"
	i__1 = *ldt + 1;
#line 708 "MB03XD.f"
	dcopy_(n, &dwork[pdw], &c__1, &t[t_dim1 + 1], &i__1);
#line 709 "MB03XD.f"
	if (*n > 1) {
#line 709 "MB03XD.f"
	    i__1 = *n - 1;
#line 709 "MB03XD.f"
	    i__2 = *n - 1;
#line 709 "MB03XD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], 
		    ldt, (ftnlen)5);
#line 709 "MB03XD.f"
	}

#line 712 "MB03XD.f"
    } else if (! wantu && wantv && ! wantg) {

/*        Workspace requirements: 7*N. */

#line 716 "MB03XD.f"
	pdw = ptaur + *n;
#line 717 "MB03XD.f"
	dlaset_("All", n, n, &c_b56, &c_b56, &v2[v2_offset], ldv2, (ftnlen)3);
/* Computing MAX */
#line 718 "MB03XD.f"
	i__2 = 0, i__3 = *n - *ilo;
#line 718 "MB03XD.f"
	i__1 = max(i__2,i__3);
/* Computing MAX */
#line 718 "MB03XD.f"
	i__5 = 0, i__6 = *n - *ilo;
#line 718 "MB03XD.f"
	i__4 = max(i__5,i__6);
#line 718 "MB03XD.f"
	i__7 = *ldwork - pdw + 1;
#line 718 "MB03XD.f"
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Rowwise", &i__1, n, &i__4, &qg[ilo1 + *ilo * qg_dim1], ldqg, 
		&qg[*ilo + ilo1 * qg_dim1], ldqg, &v1[ilo1 + v1_dim1], ldv1, &
		v2[ilo1 + v2_dim1], ldv2, &dwork[pcsr + (*ilo << 1) - 2], &
		dwork[ptaur + *ilo - 1], &dwork[pdw], &i__7, &ierr, (ftnlen)
		12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)7);
/* Computing MAX */
#line 724 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 724 "MB03XD.f"
	wrkopt = max(i__1,i__2);

#line 726 "MB03XD.f"
    } else if (wantu && wantv && ! wantg) {

/*        Workspace requirements: 8*N. */

#line 730 "MB03XD.f"
	pdw = ptaur + *n;
#line 731 "MB03XD.f"
	i__1 = *ldt + 1;
#line 731 "MB03XD.f"
	dcopy_(n, &t[t_dim1 + 1], &i__1, &dwork[pdw], &c__1);
#line 732 "MB03XD.f"
	dlacpy_("Lower", n, n, &v2[v2_offset], ldv2, &t[t_offset], ldt, (
		ftnlen)5);
#line 733 "MB03XD.f"
	dlaset_("All", n, n, &c_b56, &c_b56, &v2[v2_offset], ldv2, (ftnlen)3);
/* Computing MAX */
#line 734 "MB03XD.f"
	i__2 = 0, i__3 = *n - *ilo;
#line 734 "MB03XD.f"
	i__1 = max(i__2,i__3);
/* Computing MAX */
#line 734 "MB03XD.f"
	i__5 = 0, i__6 = *n - *ilo;
#line 734 "MB03XD.f"
	i__4 = max(i__5,i__6);
#line 734 "MB03XD.f"
	i__7 = *ldwork - pdw - *n + 1;
#line 734 "MB03XD.f"
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Rowwise", &i__1, n, &i__4, &qg[ilo1 + *ilo * qg_dim1], ldqg, 
		&u2[*ilo + ilo1 * u2_dim1], ldu2, &v1[ilo1 + v1_dim1], ldv1, &
		v2[ilo1 + v2_dim1], ldv2, &dwork[pcsr + (*ilo << 1) - 2], &
		dwork[ptaur + *ilo - 1], &dwork[pdw + *n], &i__7, &ierr, (
		ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)7);
/* Computing MAX */
#line 740 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *n] + pdw + *n - 1;
#line 740 "MB03XD.f"
	wrkopt = max(i__1,i__2);

#line 742 "MB03XD.f"
	dlacpy_("Lower", n, n, &u2[u2_offset], ldu2, &qg[qg_offset], ldqg, (
		ftnlen)5);
#line 743 "MB03XD.f"
	dlaset_("All", n, n, &c_b56, &c_b56, &u2[u2_offset], ldu2, (ftnlen)3);
#line 744 "MB03XD.f"
	i__1 = *n - *ilo + 1;
#line 744 "MB03XD.f"
	i__2 = *n - *ilo + 1;
#line 744 "MB03XD.f"
	i__3 = *ldwork - pdw - *n + 1;
#line 744 "MB03XD.f"
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Columnwise", &i__1, n, &i__2, &t[*ilo + *ilo * t_dim1], ldt, 
		&qg[*ilo + *ilo * qg_dim1], ldqg, &u1[*ilo + u1_dim1], ldu1, &
		u2[*ilo + u2_dim1], ldu2, &dwork[pcsl + (*ilo << 1) - 2], &
		dwork[ptaul + *ilo - 1], &dwork[pdw + *n], &i__3, &ierr, (
		ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)10);
/* Computing MAX */
#line 750 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *n] + pdw + *n - 1;
#line 750 "MB03XD.f"
	wrkopt = max(i__1,i__2);
#line 751 "MB03XD.f"
	i__1 = *ldt + 1;
#line 751 "MB03XD.f"
	dcopy_(n, &dwork[pdw], &c__1, &t[t_dim1 + 1], &i__1);
#line 752 "MB03XD.f"
	if (*n > 1) {
#line 752 "MB03XD.f"
	    i__1 = *n - 1;
#line 752 "MB03XD.f"
	    i__2 = *n - 1;
#line 752 "MB03XD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], 
		    ldt, (ftnlen)5);
#line 752 "MB03XD.f"
	}

#line 755 "MB03XD.f"
    } else if (wantu && ! wantv && wantg) {

/*        Workspace requirements: 6*N + N*N. */

#line 759 "MB03XD.f"
	pq = ptaur;
#line 760 "MB03XD.f"
	pdw = pq + *n * *n;
#line 761 "MB03XD.f"
	dlacpy_("Lower", n, n, &u2[u2_offset], ldu2, &dwork[pq], n, (ftnlen)5)
		;
#line 762 "MB03XD.f"
	dlaset_("All", n, n, &c_b56, &c_b56, &u2[u2_offset], ldu2, (ftnlen)3);
#line 763 "MB03XD.f"
	i__1 = *n - *ilo + 1;
#line 763 "MB03XD.f"
	i__2 = *n - *ilo + 1;
#line 763 "MB03XD.f"
	i__3 = *ldwork - pdw + 1;
#line 763 "MB03XD.f"
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Columnwise", &i__1, n, &i__2, &t[*ilo + *ilo * t_dim1], ldt, 
		&dwork[pq + (*ilo - 1) * (*n + 1)], n, &u1[*ilo + u1_dim1], 
		ldu1, &u2[*ilo + u2_dim1], ldu2, &dwork[pcsl + (*ilo << 1) - 
		2], &dwork[ptaul + *ilo - 1], &dwork[pdw], &i__3, &ierr, (
		ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)10);
/* Computing MAX */
#line 769 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 769 "MB03XD.f"
	wrkopt = max(i__1,i__2);
#line 770 "MB03XD.f"
	if (*n > 1) {
#line 770 "MB03XD.f"
	    i__1 = *n - 1;
#line 770 "MB03XD.f"
	    i__2 = *n - 1;
#line 770 "MB03XD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], 
		    ldt, (ftnlen)5);
#line 770 "MB03XD.f"
	}

#line 773 "MB03XD.f"
    } else if (! wantu && wantv && wantg) {

/*        Workspace requirements: 7*N + N*N. */

#line 777 "MB03XD.f"
	pq = ptaur + *n;
#line 778 "MB03XD.f"
	pdw = pq + *n * *n;
#line 779 "MB03XD.f"
	dlacpy_("Upper", n, n, &v2[v2_offset], ldv2, &dwork[pq], n, (ftnlen)5)
		;
#line 780 "MB03XD.f"
	dlaset_("All", n, n, &c_b56, &c_b56, &v2[v2_offset], ldv2, (ftnlen)3);
/* Computing MAX */
#line 781 "MB03XD.f"
	i__2 = 0, i__3 = *n - *ilo;
#line 781 "MB03XD.f"
	i__1 = max(i__2,i__3);
/* Computing MAX */
#line 781 "MB03XD.f"
	i__5 = 0, i__6 = *n - *ilo;
#line 781 "MB03XD.f"
	i__4 = max(i__5,i__6);
#line 781 "MB03XD.f"
	i__7 = *ldwork - pdw - *n + 1;
#line 781 "MB03XD.f"
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Rowwise", &i__1, n, &i__4, &a[ilo1 + *ilo * a_dim1], lda, &
		dwork[pq + *ilo * *n + *ilo - 1], n, &v1[ilo1 + v1_dim1], 
		ldv1, &v2[ilo1 + v2_dim1], ldv2, &dwork[pcsr + (*ilo << 1) - 
		2], &dwork[ptaur + *ilo - 1], &dwork[pdw + *n], &i__7, &ierr, 
		(ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)7);
/* Computing MAX */
#line 788 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw + *n] + pdw + *n - 1;
#line 788 "MB03XD.f"
	wrkopt = max(i__1,i__2);
#line 789 "MB03XD.f"
	if (*n > 2) {
#line 789 "MB03XD.f"
	    i__1 = *n - 2;
#line 789 "MB03XD.f"
	    i__2 = *n - 2;
#line 789 "MB03XD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &a[a_dim1 + 3], 
		    lda, (ftnlen)5);
#line 789 "MB03XD.f"
	}

#line 792 "MB03XD.f"
    } else if (wantu && wantv && wantg) {

/*        Workspace requirements: 6*N + N*N. */

#line 796 "MB03XD.f"
	pdw = ptaur + *n;
#line 797 "MB03XD.f"
	dlaset_("All", n, n, &c_b56, &c_b56, &v2[v2_offset], ldv2, (ftnlen)3);
/* Computing MAX */
#line 798 "MB03XD.f"
	i__2 = 0, i__3 = *n - *ilo;
#line 798 "MB03XD.f"
	i__1 = max(i__2,i__3);
/* Computing MAX */
#line 798 "MB03XD.f"
	i__5 = 0, i__6 = *n - *ilo;
#line 798 "MB03XD.f"
	i__4 = max(i__5,i__6);
#line 798 "MB03XD.f"
	i__7 = *ldwork - pdw + 1;
#line 798 "MB03XD.f"
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Rowwise", &i__1, n, &i__4, &a[ilo1 + *ilo * a_dim1], lda, &
		u2[*ilo + ilo1 * u2_dim1], ldu2, &v1[ilo1 + v1_dim1], ldv1, &
		v2[ilo1 + v2_dim1], ldv2, &dwork[pcsr + (*ilo << 1) - 2], &
		dwork[ptaur + *ilo - 1], &dwork[pdw], &i__7, &ierr, (ftnlen)
		12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)7);
/* Computing MAX */
#line 804 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 804 "MB03XD.f"
	wrkopt = max(i__1,i__2);

#line 806 "MB03XD.f"
	pq = ptaur;
#line 807 "MB03XD.f"
	pdw = pq + *n * *n;
#line 808 "MB03XD.f"
	dlacpy_("Lower", n, n, &u2[u2_offset], ldu2, &dwork[pq], n, (ftnlen)5)
		;
#line 809 "MB03XD.f"
	dlaset_("All", n, n, &c_b56, &c_b56, &u2[u2_offset], ldu2, (ftnlen)3);
#line 810 "MB03XD.f"
	i__1 = *n - *ilo + 1;
#line 810 "MB03XD.f"
	i__2 = *n - *ilo + 1;
#line 810 "MB03XD.f"
	i__3 = *ldwork - pdw + 1;
#line 810 "MB03XD.f"
	mb04qb_("No Transpose", "No Transpose", "No Transpose", "Columnwise", 
		"Columnwise", &i__1, n, &i__2, &t[*ilo + *ilo * t_dim1], ldt, 
		&dwork[pq + (*ilo - 1) * (*n + 1)], n, &u1[*ilo + u1_dim1], 
		ldu1, &u2[*ilo + u2_dim1], ldu2, &dwork[pcsl + (*ilo << 1) - 
		2], &dwork[ptaul + *ilo - 1], &dwork[pdw], &i__3, &ierr, (
		ftnlen)12, (ftnlen)12, (ftnlen)12, (ftnlen)10, (ftnlen)10);
/* Computing MAX */
#line 816 "MB03XD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[pdw] + pdw - 1;
#line 816 "MB03XD.f"
	wrkopt = max(i__1,i__2);
#line 817 "MB03XD.f"
	if (*n > 2) {
#line 817 "MB03XD.f"
	    i__1 = *n - 2;
#line 817 "MB03XD.f"
	    i__2 = *n - 2;
#line 817 "MB03XD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &a[a_dim1 + 3], 
		    lda, (ftnlen)5);
#line 817 "MB03XD.f"
	}
#line 819 "MB03XD.f"
	if (*n > 1) {
#line 819 "MB03XD.f"
	    i__1 = *n - 1;
#line 819 "MB03XD.f"
	    i__2 = *n - 1;
#line 819 "MB03XD.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b56, &c_b56, &t[t_dim1 + 2], 
		    ldt, (ftnlen)5);
#line 819 "MB03XD.f"
	}
#line 821 "MB03XD.f"
    }

#line 823 "MB03XD.f"
    dwork[1] = (doublereal) wrkopt;
#line 824 "MB03XD.f"
    return 0;
/* *** Last line of MB03XD *** */
} /* mb03xd_ */

