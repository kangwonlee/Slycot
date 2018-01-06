#line 1 "MB03XP.f"
/* MB03XP.f -- translated by f2c (version 20100827).
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

#line 1 "MB03XP.f"
/* Table of constant values */

static doublereal c_b12 = 0.;
static doublereal c_b13 = 1.;
static integer c__4 = 4;
static integer c__2 = 2;
static integer c__8 = 8;
static integer c__3 = 3;
static integer c__15 = 15;
static logical c_false = FALSE_;
static integer c__1 = 1;
static doublereal c_b89 = -1.;

/* Subroutine */ int mb03xp_(char *job, char *compq, char *compz, integer *n, 
	integer *ilo, integer *ihi, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *q, integer *ldq, doublereal *z__, 
	integer *ldz, doublereal *alphar, doublereal *alphai, doublereal *
	beta, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	job_len, ftnlen compq_len, ftnlen compz_len)
{
    /* System generated locals */
    address a__1[2];
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3[2], i__4, i__5;
    doublereal d__1, d__2;
    char ch__1[2];

    /* Builtin functions */
    /* Subroutine */ int s_cat(char *, char **, integer *, integer *, ftnlen);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal v[51];
    static integer i1, i2;
    static doublereal as[225]	/* was [15][15] */, bs[225]	/* was [15][
	    15] */;
    static integer kk, nh, nr, ns, nv, pv2, pv3, dum, itn, its;
    static doublereal ulp, tst;
    static integer maxb, ierr;
    static doublereal unfl, temp, ovfl, tauv, tauw;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static integer iseed[4];
    extern integer ue01md_(integer *, char *, char *, integer *, integer *, 
	    integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb03ya_(logical *, logical *, logical *, 
	    integer *, integer *, integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *), mb03yd_(logical *,
	     logical *, logical *, integer *, integer *, integer *, integer *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer itemp;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical initq;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical wantq;
    extern /* Subroutine */ int dtrmv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static logical initz, wantt, wantz;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dlarfx_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, ftnlen), 
	    dlarnv_(integer *, integer *, integer *, doublereal *);
    static doublereal smlnum;


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

/*     To compute the periodic Schur decomposition and the eigenvalues of */
/*     a product of matrices, H = A*B, with A upper Hessenberg and B */
/*     upper triangular without evaluating any part of the product. */
/*     Specifically, the matrices Q and Z are computed, so that */

/*          Q' * A * Z = S,    Z' * B * Q = T */

/*     where S is in real Schur form, and T is upper triangular. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates whether the user wishes to compute the full */
/*             Schur form or the eigenvalues only, as follows: */
/*             = 'E':  Compute the eigenvalues only; */
/*             = 'S':  compute the factors S and T of the full */
/*                     Schur form. */

/*     COMPQ   CHARACTER*1 */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Q as follows: */
/*             = 'N':  The matrix Q is not required; */
/*             = 'I':  Q is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Q is returned; */
/*             = 'V':  Q must contain an orthogonal matrix U on entry, */
/*                     and the product U*Q is returned. */

/*     COMPZ   CHARACTER*1 */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Z as follows: */
/*             = 'N':  The matrix Z is not required; */
/*             = 'I':  Z is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Z is returned; */
/*             = 'V':  Z must contain an orthogonal matrix U on entry, */
/*                     and the product U*Z is returned. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and B. N >= 0. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that the matrices A and B are already upper */
/*             triangular in rows and columns 1:ILO-1 and IHI+1:N. */
/*             The routine works primarily with the submatrices in rows */
/*             and columns ILO to IHI, but applies the transformations to */
/*             all the rows and columns of the matrices A and B, if */
/*             JOB = 'S'. */
/*             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array A must */
/*             contain the upper Hessenberg matrix A. */
/*             On exit, if JOB = 'S', the leading N-by-N part of this */
/*             array is upper quasi-triangular with any 2-by-2 diagonal */
/*             blocks corresponding to a pair of complex conjugated */
/*             eigenvalues. */
/*             If JOB = 'E', the diagonal elements and 2-by-2 diagonal */
/*             blocks of A will be correct, but the remaining parts of A */
/*             are unspecified on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,N) */
/*             On entry, the leading N-by-N part of this array B must */
/*             contain the upper triangular matrix B. */
/*             On exit, if JOB = 'S', the leading N-by-N part of this */
/*             array contains the transformed upper triangular matrix. */
/*             2-by-2 blocks in B corresponding to 2-by-2 blocks in A */
/*             will be reduced to positive diagonal form. (I.e., if */
/*             A(j+1,j) is non-zero, then B(j+1,j)=B(j,j+1)=0 and B(j,j) */
/*             and B(j+1,j+1) will be positive.) */
/*             If JOB = 'E', the elements corresponding to diagonal */
/*             elements and 2-by-2 diagonal blocks in A will be correct, */
/*             but the remaining parts of B are unspecified on exit. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             On entry, if COMPQ = 'V', then the leading N-by-N part of */
/*             this array must contain a matrix Q which is assumed to be */
/*             equal to the unit matrix except for the submatrix */
/*             Q(ILO:IHI,ILO:IHI). */
/*             If COMPQ = 'I', Q need not be set on entry. */
/*             On exit, if COMPQ = 'V' or COMPQ = 'I' the leading N-by-N */
/*             part of this array contains the transformation matrix */
/*             which produced the Schur form. */
/*             If COMPQ = 'N', Q is not referenced. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q.  LDQ >= 1. */
/*             If COMPQ <> 'N', LDQ >= MAX(1,N). */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N) */
/*             On entry, if COMPZ = 'V', then the leading N-by-N part of */
/*             this array must contain a matrix Z which is assumed to be */
/*             equal to the unit matrix except for the submatrix */
/*             Z(ILO:IHI,ILO:IHI). */
/*             If COMPZ = 'I', Z need not be set on entry. */
/*             On exit, if COMPZ = 'V' or COMPZ = 'I' the leading N-by-N */
/*             part of this array contains the transformation matrix */
/*             which produced the Schur form. */
/*             If COMPZ = 'N', Z is not referenced. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= 1. */
/*             If COMPZ <> 'N', LDZ >= MAX(1,N). */

/*     ALPHAR  (output) DOUBLE PRECISION array, dimension (N) */
/*     ALPHAI  (output) DOUBLE PRECISION array, dimension (N) */
/*     BETA    (output) DOUBLE PRECISION array, dimension (N) */
/*             The i-th (1 <= i <= N) computed eigenvalue is given by */
/*             BETA(I) * ( ALPHAR(I) + sqrt(-1)*ALPHAI(I) ). If two */
/*             eigenvalues are computed as a complex conjugate pair, */
/*             they are stored in consecutive elements of ALPHAR, ALPHAI */
/*             and BETA. If JOB = 'S', the eigenvalues are stored in the */
/*             same order as on the diagonales of the Schur forms of A */
/*             and B. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -19,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, then MB03XP failed to compute the Schur */
/*                   form in a total of 30*(IHI-ILO+1) iterations; */
/*                   elements 1:ilo-1 and i+1:n of ALPHAR, ALPHAI and */
/*                   BETA contain successfully computed eigenvalues. */

/*     METHOD */

/*     The implemented algorithm is a multi-shift version of the periodic */
/*     QR algorithm described in [1,3] with some minor modifications */
/*     proposed in [2]. */

/*     REFERENCES */

/*     [1] Bojanczyk, A.W., Golub, G.H., and Van Dooren, P. */
/*         The periodic Schur decomposition: Algorithms and applications. */
/*         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42, */
/*         1992. */

/*     [2] Kressner, D. */
/*         An efficient and reliable implementation of the periodic QZ */
/*         algorithm. Proc. of the IFAC Workshop on Periodic Control */
/*         Systems, pp. 187-192, 2001. */

/*     [3] Van Loan, C. */
/*         Generalized Singular Values with Algorithms and Applications. */
/*         Ph. D. Thesis, University of Michigan, 1973. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O(N**3) floating point operations and is */
/*     backward stable. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHGPQR). */

/*     KEYWORDS */

/*     Eigenvalue, eigenvalue decomposition, Hessenberg form, orthogonal */
/*     transformation, (periodic) Schur form */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 249 "MB03XP.f"
    /* Parameter adjustments */
#line 249 "MB03XP.f"
    a_dim1 = *lda;
#line 249 "MB03XP.f"
    a_offset = 1 + a_dim1;
#line 249 "MB03XP.f"
    a -= a_offset;
#line 249 "MB03XP.f"
    b_dim1 = *ldb;
#line 249 "MB03XP.f"
    b_offset = 1 + b_dim1;
#line 249 "MB03XP.f"
    b -= b_offset;
#line 249 "MB03XP.f"
    q_dim1 = *ldq;
#line 249 "MB03XP.f"
    q_offset = 1 + q_dim1;
#line 249 "MB03XP.f"
    q -= q_offset;
#line 249 "MB03XP.f"
    z_dim1 = *ldz;
#line 249 "MB03XP.f"
    z_offset = 1 + z_dim1;
#line 249 "MB03XP.f"
    z__ -= z_offset;
#line 249 "MB03XP.f"
    --alphar;
#line 249 "MB03XP.f"
    --alphai;
#line 249 "MB03XP.f"
    --beta;
#line 249 "MB03XP.f"
    --dwork;
#line 249 "MB03XP.f"

#line 249 "MB03XP.f"
    /* Function Body */
#line 249 "MB03XP.f"
    wantt = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 250 "MB03XP.f"
    initq = lsame_(compq, "I", (ftnlen)1, (ftnlen)1);
#line 251 "MB03XP.f"
    wantq = initq || lsame_(compq, "V", (ftnlen)1, (ftnlen)1);
#line 252 "MB03XP.f"
    initz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
#line 253 "MB03XP.f"
    wantz = initz || lsame_(compz, "V", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 257 "MB03XP.f"
    *info = 0;
#line 258 "MB03XP.f"
    if (! lsame_(job, "E", (ftnlen)1, (ftnlen)1) && ! wantt) {
#line 259 "MB03XP.f"
	*info = -1;
#line 260 "MB03XP.f"
    } else if (! lsame_(compq, "N", (ftnlen)1, (ftnlen)1) && ! wantq) {
#line 261 "MB03XP.f"
	*info = -2;
#line 262 "MB03XP.f"
    } else if (! lsame_(compz, "N", (ftnlen)1, (ftnlen)1) && ! wantz) {
#line 263 "MB03XP.f"
	*info = -3;
#line 264 "MB03XP.f"
    } else if (*n < 0) {
#line 265 "MB03XP.f"
	*info = -4;
#line 266 "MB03XP.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 267 "MB03XP.f"
	*info = -5;
#line 268 "MB03XP.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 269 "MB03XP.f"
	*info = -6;
#line 270 "MB03XP.f"
    } else if (*lda < max(1,*n)) {
#line 271 "MB03XP.f"
	*info = -8;
#line 272 "MB03XP.f"
    } else if (*ldb < max(1,*n)) {
#line 273 "MB03XP.f"
	*info = -10;
#line 274 "MB03XP.f"
    } else if (*ldq < 1 || wantq && *ldq < *n) {
#line 275 "MB03XP.f"
	*info = -12;
#line 276 "MB03XP.f"
    } else if (*ldz < 1 || wantz && *ldz < *n) {
#line 277 "MB03XP.f"
	*info = -14;
#line 278 "MB03XP.f"
    } else if (*ldwork < max(1,*n)) {
#line 279 "MB03XP.f"
	dwork[1] = (doublereal) max(1,*n);
#line 280 "MB03XP.f"
	*info = -19;
#line 281 "MB03XP.f"
    }

/*     Return if there were illegal values. */

#line 285 "MB03XP.f"
    if (*info != 0) {
#line 286 "MB03XP.f"
	i__1 = -(*info);
#line 286 "MB03XP.f"
	xerbla_("MB03XP", &i__1, (ftnlen)6);
#line 287 "MB03XP.f"
	return 0;
#line 288 "MB03XP.f"
    }

/*     Initialize Q and Z, if necessary. */

#line 292 "MB03XP.f"
    if (initq) {
#line 292 "MB03XP.f"
	dlaset_("All", n, n, &c_b12, &c_b13, &q[q_offset], ldq, (ftnlen)3);
#line 292 "MB03XP.f"
    }
#line 294 "MB03XP.f"
    if (initz) {
#line 294 "MB03XP.f"
	dlaset_("All", n, n, &c_b12, &c_b13, &z__[z_offset], ldz, (ftnlen)3);
#line 294 "MB03XP.f"
    }

/*     Store isolated eigenvalues and standardize B. */

/*     FOR I = [1:ILO-1, IHI+1:N] */
#line 300 "MB03XP.f"
    i__ = 1;
#line 301 "MB03XP.f"
L10:
#line 302 "MB03XP.f"
    if (i__ == *ilo) {
#line 303 "MB03XP.f"
	i__ = *ihi + 1;
#line 304 "MB03XP.f"
    }
#line 305 "MB03XP.f"
    if (i__ <= *n) {
#line 306 "MB03XP.f"
	if (b[i__ + i__ * b_dim1] < 0.) {
#line 307 "MB03XP.f"
	    if (wantt) {
#line 308 "MB03XP.f"
		i__1 = i__;
#line 308 "MB03XP.f"
		for (k = *ilo; k <= i__1; ++k) {
#line 309 "MB03XP.f"
		    b[k + i__ * b_dim1] = -b[k + i__ * b_dim1];
#line 310 "MB03XP.f"
/* L20: */
#line 310 "MB03XP.f"
		}
#line 311 "MB03XP.f"
		i__1 = *ihi;
#line 311 "MB03XP.f"
		for (k = i__; k <= i__1; ++k) {
#line 312 "MB03XP.f"
		    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1];
#line 313 "MB03XP.f"
/* L30: */
#line 313 "MB03XP.f"
		}
#line 314 "MB03XP.f"
	    } else {
#line 315 "MB03XP.f"
		b[i__ + i__ * b_dim1] = -b[i__ + i__ * b_dim1];
#line 316 "MB03XP.f"
		a[i__ + i__ * a_dim1] = -a[i__ + i__ * a_dim1];
#line 317 "MB03XP.f"
	    }
#line 318 "MB03XP.f"
	    if (wantq) {
#line 319 "MB03XP.f"
		i__1 = *ihi;
#line 319 "MB03XP.f"
		for (k = *ilo; k <= i__1; ++k) {
#line 320 "MB03XP.f"
		    q[k + i__ * q_dim1] = -q[k + i__ * q_dim1];
#line 321 "MB03XP.f"
/* L40: */
#line 321 "MB03XP.f"
		}
#line 322 "MB03XP.f"
	    }
#line 323 "MB03XP.f"
	}
#line 324 "MB03XP.f"
	alphar[i__] = a[i__ + i__ * a_dim1];
#line 325 "MB03XP.f"
	alphai[i__] = 0.;
#line 326 "MB03XP.f"
	beta[i__] = b[i__ + i__ * b_dim1];
#line 327 "MB03XP.f"
	++i__;
/*        END FOR */
#line 329 "MB03XP.f"
	goto L10;
#line 330 "MB03XP.f"
    }

/*     Quick return if possible. */

#line 334 "MB03XP.f"
    if (*n == 0 || *ilo == *ihi + 1) {
#line 335 "MB03XP.f"
	dwork[1] = 1.;
#line 336 "MB03XP.f"
	return 0;
#line 337 "MB03XP.f"
    }

/*     Set rows and coloms ILO to IHI of B (A) to zero below the first */
/*     (sub)diagonal. */

#line 342 "MB03XP.f"
    i__1 = *ihi - 2;
#line 342 "MB03XP.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 343 "MB03XP.f"
	i__2 = *n;
#line 343 "MB03XP.f"
	for (i__ = j + 2; i__ <= i__2; ++i__) {
#line 344 "MB03XP.f"
	    a[i__ + j * a_dim1] = 0.;
#line 345 "MB03XP.f"
/* L50: */
#line 345 "MB03XP.f"
	}
#line 346 "MB03XP.f"
/* L60: */
#line 346 "MB03XP.f"
    }
#line 347 "MB03XP.f"
    i__1 = *ihi - 1;
#line 347 "MB03XP.f"
    for (j = *ilo; j <= i__1; ++j) {
#line 348 "MB03XP.f"
	i__2 = *n;
#line 348 "MB03XP.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 349 "MB03XP.f"
	    b[i__ + j * b_dim1] = 0.;
#line 350 "MB03XP.f"
/* L70: */
#line 350 "MB03XP.f"
	}
#line 351 "MB03XP.f"
/* L80: */
#line 351 "MB03XP.f"
    }
#line 352 "MB03XP.f"
    nh = *ihi - *ilo + 1;

/*     Suboptimal choice of the number of shifts. */

#line 356 "MB03XP.f"
    if (wantq) {
/* Writing concatenation */
#line 357 "MB03XP.f"
	i__3[0] = 1, a__1[0] = job;
#line 357 "MB03XP.f"
	i__3[1] = 1, a__1[1] = compq;
#line 357 "MB03XP.f"
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 357 "MB03XP.f"
	ns = ue01md_(&c__4, "MB03XP", ch__1, n, ilo, ihi, (ftnlen)6, (ftnlen)
		2);
/* Writing concatenation */
#line 358 "MB03XP.f"
	i__3[0] = 1, a__1[0] = job;
#line 358 "MB03XP.f"
	i__3[1] = 1, a__1[1] = compq;
#line 358 "MB03XP.f"
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 358 "MB03XP.f"
	maxb = ue01md_(&c__8, "MB03XP", ch__1, n, ilo, ihi, (ftnlen)6, (
		ftnlen)2);
#line 359 "MB03XP.f"
    } else {
/* Writing concatenation */
#line 360 "MB03XP.f"
	i__3[0] = 1, a__1[0] = job;
#line 360 "MB03XP.f"
	i__3[1] = 1, a__1[1] = compz;
#line 360 "MB03XP.f"
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 360 "MB03XP.f"
	ns = ue01md_(&c__4, "MB03XP", ch__1, n, ilo, ihi, (ftnlen)6, (ftnlen)
		2);
/* Writing concatenation */
#line 361 "MB03XP.f"
	i__3[0] = 1, a__1[0] = job;
#line 361 "MB03XP.f"
	i__3[1] = 1, a__1[1] = compz;
#line 361 "MB03XP.f"
	s_cat(ch__1, a__1, i__3, &c__2, (ftnlen)2);
#line 361 "MB03XP.f"
	maxb = ue01md_(&c__8, "MB03XP", ch__1, n, ilo, ihi, (ftnlen)6, (
		ftnlen)2);
#line 362 "MB03XP.f"
    }

#line 364 "MB03XP.f"
    if (ns <= 2 || ns > nh || maxb >= nh) {

/*        Standard double-shift product QR. */

#line 368 "MB03XP.f"
	mb03yd_(&wantt, &wantq, &wantz, n, ilo, ihi, ilo, ihi, &a[a_offset], 
		lda, &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], 
		ldz, &alphar[1], &alphai[1], &beta[1], &dwork[1], ldwork, 
		info);
#line 371 "MB03XP.f"
	return 0;
#line 372 "MB03XP.f"
    }
#line 373 "MB03XP.f"
    maxb = max(3,maxb);
/* Computing MIN */
#line 374 "MB03XP.f"
    i__1 = min(ns,maxb);
#line 374 "MB03XP.f"
    ns = min(i__1,15);

/*     Set machine-dependent constants for the stopping criterion. */
/*     If max(norm(A),norm(B)) <= sqrt(OVFL), then overflow should not */
/*     occur. */

#line 380 "MB03XP.f"
    unfl = dlamch_("Safe minimum", (ftnlen)12);
#line 381 "MB03XP.f"
    ovfl = 1. / unfl;
#line 382 "MB03XP.f"
    dlabad_(&unfl, &ovfl);
#line 383 "MB03XP.f"
    ulp = dlamch_("Precision", (ftnlen)9);
#line 384 "MB03XP.f"
    smlnum = unfl * ((doublereal) nh / ulp);

/*     I1 and I2 are the indices of the first rows and last columns of */
/*     A and B to which transformations must be applied. */

#line 389 "MB03XP.f"
    if (wantt) {
#line 390 "MB03XP.f"
	i1 = 1;
#line 391 "MB03XP.f"
	i2 = *n;
#line 392 "MB03XP.f"
    }
#line 393 "MB03XP.f"
    iseed[0] = 1;
#line 394 "MB03XP.f"
    iseed[1] = 0;
#line 395 "MB03XP.f"
    iseed[2] = 0;
#line 396 "MB03XP.f"
    iseed[3] = 1;

/*     ITN is the maximal number of QR iterations. */

#line 400 "MB03XP.f"
    itn = nh * 30;
#line 401 "MB03XP.f"
    dum = 0;

/*     Main loop. Eigenvalues I+1:IHI have converged. Either L = ILO */
/*     or A(L,L-1) is negligible. */

#line 406 "MB03XP.f"
    i__ = *ihi;
#line 407 "MB03XP.f"
L90:
#line 408 "MB03XP.f"
    l = *ilo;
#line 409 "MB03XP.f"
    if (i__ < *ilo) {
#line 409 "MB03XP.f"
	goto L210;
#line 409 "MB03XP.f"
    }

#line 412 "MB03XP.f"
    i__1 = itn;
#line 412 "MB03XP.f"
    for (its = 0; its <= i__1; ++its) {
#line 413 "MB03XP.f"
	dum += (*ihi - *ilo) * (*ihi - *ilo);

/*        Look for deflations in A. */

#line 417 "MB03XP.f"
	i__2 = l + 1;
#line 417 "MB03XP.f"
	for (k = i__; k >= i__2; --k) {
#line 418 "MB03XP.f"
	    tst = (d__1 = a[k - 1 + (k - 1) * a_dim1], abs(d__1)) + (d__2 = a[
		    k + k * a_dim1], abs(d__2));
#line 419 "MB03XP.f"
	    if (tst == 0.) {
#line 419 "MB03XP.f"
		i__4 = i__ - l + 1;
#line 419 "MB03XP.f"
		tst = dlanhs_("1", &i__4, &a[l + l * a_dim1], lda, &dwork[1], 
			(ftnlen)1);
#line 419 "MB03XP.f"
	    }
/* Computing MAX */
#line 421 "MB03XP.f"
	    d__2 = ulp * tst;
#line 421 "MB03XP.f"
	    if ((d__1 = a[k + (k - 1) * a_dim1], abs(d__1)) <= max(d__2,
		    smlnum)) {
#line 421 "MB03XP.f"
		goto L110;
#line 421 "MB03XP.f"
	    }
#line 423 "MB03XP.f"
/* L100: */
#line 423 "MB03XP.f"
	}
#line 424 "MB03XP.f"
L110:

/*        Look for deflation in B if problem size is greater than 1. */

#line 428 "MB03XP.f"
	if (i__ - k >= 1) {
#line 429 "MB03XP.f"
	    i__2 = k;
#line 429 "MB03XP.f"
	    for (kk = i__; kk >= i__2; --kk) {
#line 430 "MB03XP.f"
		if (kk == i__) {
#line 431 "MB03XP.f"
		    tst = (d__1 = b[kk - 1 + kk * b_dim1], abs(d__1));
#line 432 "MB03XP.f"
		} else if (kk == k) {
#line 433 "MB03XP.f"
		    tst = (d__1 = b[kk + (kk + 1) * b_dim1], abs(d__1));
#line 434 "MB03XP.f"
		} else {
#line 435 "MB03XP.f"
		    tst = (d__1 = b[kk - 1 + kk * b_dim1], abs(d__1)) + (d__2 
			    = b[kk + (kk + 1) * b_dim1], abs(d__2));
#line 436 "MB03XP.f"
		}
#line 437 "MB03XP.f"
		if (tst == 0.) {
#line 437 "MB03XP.f"
		    i__4 = i__ - k + 1;
#line 437 "MB03XP.f"
		    tst = dlanhs_("1", &i__4, &b[k + k * b_dim1], ldb, &dwork[
			    1], (ftnlen)1);
#line 437 "MB03XP.f"
		}
/* Computing MAX */
#line 439 "MB03XP.f"
		d__2 = ulp * tst;
#line 439 "MB03XP.f"
		if ((d__1 = b[kk + kk * b_dim1], abs(d__1)) <= max(d__2,
			smlnum)) {
#line 439 "MB03XP.f"
		    goto L130;
#line 439 "MB03XP.f"
		}
#line 441 "MB03XP.f"
/* L120: */
#line 441 "MB03XP.f"
	    }
#line 442 "MB03XP.f"
	} else {
#line 443 "MB03XP.f"
	    kk = k - 1;
#line 444 "MB03XP.f"
	}
#line 445 "MB03XP.f"
L130:
#line 446 "MB03XP.f"
	if (kk >= k) {

/*           B has an element close to zero at position (KK,KK). */

#line 450 "MB03XP.f"
	    b[kk + kk * b_dim1] = 0.;
#line 451 "MB03XP.f"
	    mb03ya_(&wantt, &wantq, &wantz, n, &k, &i__, ilo, ihi, &kk, &a[
		    a_offset], lda, &b[b_offset], ldb, &q[q_offset], ldq, &
		    z__[z_offset], ldz, info);
#line 453 "MB03XP.f"
	    k = kk + 1;
#line 454 "MB03XP.f"
	}
#line 455 "MB03XP.f"
	l = k;
#line 456 "MB03XP.f"
	if (l > *ilo) {

/*           A(L,L-1) is negligible. */

#line 460 "MB03XP.f"
	    a[l + (l - 1) * a_dim1] = 0.;
#line 461 "MB03XP.f"
	}

/*        Exit from loop if a submatrix of order <= MAXB has split off. */

#line 465 "MB03XP.f"
	if (l >= i__ - maxb + 1) {
#line 465 "MB03XP.f"
	    goto L200;
#line 465 "MB03XP.f"
	}

/*        The active submatrices are now in rows and columns L:I. */

#line 470 "MB03XP.f"
	if (! wantt) {
#line 471 "MB03XP.f"
	    i1 = l;
#line 472 "MB03XP.f"
	    i2 = i__;
#line 473 "MB03XP.f"
	}
#line 474 "MB03XP.f"
	if (its == 10 || its == 20) {

/*           Exceptional shift. The first column of the shift polynomial */
/*           is a pseudo-random vector. */

#line 479 "MB03XP.f"
	    i__2 = ns + 1;
#line 479 "MB03XP.f"
	    dlarnv_(&c__3, iseed, &i__2, v);
#line 480 "MB03XP.f"
	} else {

/*           Use eigenvalues of trailing submatrix as shifts. */

#line 484 "MB03XP.f"
	    dlacpy_("Full", &ns, &ns, &a[i__ - ns + 1 + (i__ - ns + 1) * 
		    a_dim1], lda, as, &c__15, (ftnlen)4);
#line 486 "MB03XP.f"
	    dlacpy_("Full", &ns, &ns, &b[i__ - ns + 1 + (i__ - ns + 1) * 
		    b_dim1], ldb, bs, &c__15, (ftnlen)4);
#line 488 "MB03XP.f"
	    mb03yd_(&c_false, &c_false, &c_false, &ns, &c__1, &ns, &c__1, &ns,
		     as, &c__15, bs, &c__15, &q[q_offset], ldq, &z__[z_offset]
		    , ldz, &alphar[i__ - ns + 1], &alphai[i__ - ns + 1], &
		    beta[i__ - ns + 1], &dwork[1], ldwork, &ierr);
#line 492 "MB03XP.f"
	}

/*        Compute the nonzero elements of the first column of */
/*        (A*B-w(1)) (A*B-w(2)) .. (A*B-w(ns)). */

#line 497 "MB03XP.f"
	v[0] = 1.;
#line 498 "MB03XP.f"
	nv = 1;
/*        WHILE NV <= NS */
#line 500 "MB03XP.f"
L140:
#line 501 "MB03XP.f"
	if (nv <= ns) {
#line 502 "MB03XP.f"
	    if (nv == ns || as[nv + 1 + nv * 15 - 16] == 0.) {

/*              Real shift. */

#line 506 "MB03XP.f"
		v[nv] = 0.;
#line 507 "MB03XP.f"
		pv2 = nv + 2;
#line 508 "MB03XP.f"
		dcopy_(&nv, v, &c__1, &v[pv2 - 1], &c__1);
#line 509 "MB03XP.f"
		dtrmv_("Upper", "No transpose", "No unit diagonal", &nv, &b[l 
			+ l * b_dim1], ldb, &v[pv2 - 1], &c__1, (ftnlen)5, (
			ftnlen)12, (ftnlen)16);
#line 511 "MB03XP.f"
		dscal_(&nv, &bs[nv + nv * 15 - 16], v, &c__1);
#line 512 "MB03XP.f"
		i__2 = (nv << 1) + 1;
#line 512 "MB03XP.f"
		itemp = idamax_(&i__2, v, &c__1);
/* Computing MAX */
#line 513 "MB03XP.f"
		d__2 = (d__1 = v[itemp - 1], abs(d__1));
#line 513 "MB03XP.f"
		temp = 1. / max(d__2,smlnum);
#line 514 "MB03XP.f"
		i__2 = (nv << 1) + 1;
#line 514 "MB03XP.f"
		dscal_(&i__2, &temp, v, &c__1);
#line 515 "MB03XP.f"
		i__2 = nv + 1;
#line 515 "MB03XP.f"
		d__1 = -as[nv + nv * 15 - 16];
#line 515 "MB03XP.f"
		dgemv_("No transpose", &i__2, &nv, &c_b13, &a[l + l * a_dim1],
			 lda, &v[pv2 - 1], &c__1, &d__1, v, &c__1, (ftnlen)12)
			;
#line 517 "MB03XP.f"
		++nv;
#line 518 "MB03XP.f"
	    } else {

/*              Double shift using a product formulation of the shift */
/*              polynomial [2]. */

#line 523 "MB03XP.f"
		v[nv] = 0.;
#line 524 "MB03XP.f"
		v[nv + 1] = 0.;
#line 525 "MB03XP.f"
		pv2 = nv + 3;
#line 526 "MB03XP.f"
		pv3 = (nv << 1) + 5;
#line 527 "MB03XP.f"
		i__2 = nv + 2;
#line 527 "MB03XP.f"
		dcopy_(&i__2, v, &c__1, &v[pv2 - 1], &c__1);
#line 528 "MB03XP.f"
		i__2 = nv + 1;
#line 528 "MB03XP.f"
		dcopy_(&i__2, v, &c__1, &v[pv3 - 1], &c__1);
#line 529 "MB03XP.f"
		dscal_(&nv, &bs[nv + 1 + (nv + 1) * 15 - 16], &v[pv2 - 1], &
			c__1);
#line 530 "MB03XP.f"
		dtrmv_("Upper", "No transpose", "No unit diagonal", &nv, &b[l 
			+ l * b_dim1], ldb, &v[pv3 - 1], &c__1, (ftnlen)5, (
			ftnlen)12, (ftnlen)16);
#line 532 "MB03XP.f"
		i__2 = (nv << 1) + 3;
#line 532 "MB03XP.f"
		itemp = idamax_(&i__2, &v[pv2 - 1], &c__1);
/* Computing MAX */
#line 533 "MB03XP.f"
		d__2 = (d__1 = v[pv2 + itemp - 2], abs(d__1));
#line 533 "MB03XP.f"
		temp = 1. / max(d__2,smlnum);
#line 534 "MB03XP.f"
		i__2 = (nv << 1) + 3;
#line 534 "MB03XP.f"
		dscal_(&i__2, &temp, &v[pv2 - 1], &c__1);

#line 536 "MB03XP.f"
		dcopy_(&nv, &v[pv2 - 1], &c__1, v, &c__1);
#line 537 "MB03XP.f"
		i__2 = nv + 1;
#line 537 "MB03XP.f"
		dgemv_("No transpose", &i__2, &nv, &c_b89, &a[l + l * a_dim1],
			 lda, &v[pv3 - 1], &c__1, &as[nv + 1 + (nv + 1) * 15 
			- 16], &v[pv2 - 1], &c__1, (ftnlen)12);
#line 539 "MB03XP.f"
		dscal_(&nv, &as[nv + (nv + 1) * 15 - 16], v, &c__1);
#line 540 "MB03XP.f"
		i__2 = (nv << 1) + 3;
#line 540 "MB03XP.f"
		itemp = idamax_(&i__2, v, &c__1);
/* Computing MAX */
#line 541 "MB03XP.f"
		d__2 = (d__1 = v[itemp - 1], abs(d__1));
#line 541 "MB03XP.f"
		temp = 1. / max(d__2,smlnum);
#line 542 "MB03XP.f"
		i__2 = (nv << 1) + 3;
#line 542 "MB03XP.f"
		dscal_(&i__2, &temp, v, &c__1);

#line 544 "MB03XP.f"
		d__1 = -as[nv + 1 + nv * 15 - 16];
#line 544 "MB03XP.f"
		dscal_(&nv, &d__1, v, &c__1);
#line 545 "MB03XP.f"
		i__2 = nv + 1;
#line 545 "MB03XP.f"
		daxpy_(&i__2, &as[nv + nv * 15 - 16], &v[pv2 - 1], &c__1, v, &
			c__1);
#line 546 "MB03XP.f"
		i__2 = (nv << 1) + 3;
#line 546 "MB03XP.f"
		itemp = idamax_(&i__2, v, &c__1);
/* Computing MAX */
#line 547 "MB03XP.f"
		d__2 = (d__1 = v[itemp - 1], abs(d__1));
#line 547 "MB03XP.f"
		temp = 1. / max(d__2,smlnum);
#line 548 "MB03XP.f"
		i__2 = (nv << 1) + 3;
#line 548 "MB03XP.f"
		dscal_(&i__2, &temp, v, &c__1);

#line 550 "MB03XP.f"
		i__2 = nv + 1;
#line 550 "MB03XP.f"
		dscal_(&i__2, &bs[nv + nv * 15 - 16], v, &c__1);
#line 551 "MB03XP.f"
		i__2 = nv + 1;
#line 551 "MB03XP.f"
		dtrmv_("Upper", "No transpose", "No unit diagonal", &i__2, &b[
			l + l * b_dim1], ldb, &v[pv2 - 1], &c__1, (ftnlen)5, (
			ftnlen)12, (ftnlen)16);
#line 553 "MB03XP.f"
		i__2 = (nv << 1) + 3;
#line 553 "MB03XP.f"
		itemp = idamax_(&i__2, v, &c__1);
/* Computing MAX */
#line 554 "MB03XP.f"
		d__2 = (d__1 = v[itemp - 1], abs(d__1));
#line 554 "MB03XP.f"
		temp = 1. / max(d__2,smlnum);
#line 555 "MB03XP.f"
		i__2 = (nv << 1) + 3;
#line 555 "MB03XP.f"
		dscal_(&i__2, &temp, v, &c__1);

#line 557 "MB03XP.f"
		i__2 = nv + 2;
#line 557 "MB03XP.f"
		i__4 = nv + 1;
#line 557 "MB03XP.f"
		dgemv_("No transpose", &i__2, &i__4, &c_b89, &a[l + l * 
			a_dim1], lda, &v[pv2 - 1], &c__1, &c_b13, v, &c__1, (
			ftnlen)12);
#line 559 "MB03XP.f"
		nv += 2;
#line 560 "MB03XP.f"
	    }
#line 561 "MB03XP.f"
	    itemp = idamax_(&nv, v, &c__1);
#line 562 "MB03XP.f"
	    temp = (d__1 = v[itemp - 1], abs(d__1));
#line 563 "MB03XP.f"
	    if (temp == 0.) {
#line 564 "MB03XP.f"
		v[0] = 1.;
#line 565 "MB03XP.f"
		i__2 = nv;
#line 565 "MB03XP.f"
		for (k = 2; k <= i__2; ++k) {
#line 566 "MB03XP.f"
		    v[k - 1] = 0.;
#line 567 "MB03XP.f"
/* L150: */
#line 567 "MB03XP.f"
		}
#line 568 "MB03XP.f"
	    } else {
#line 569 "MB03XP.f"
		temp = max(temp,smlnum);
#line 570 "MB03XP.f"
		d__1 = 1. / temp;
#line 570 "MB03XP.f"
		dscal_(&nv, &d__1, v, &c__1);
#line 571 "MB03XP.f"
	    }
#line 572 "MB03XP.f"
	    goto L140;
/*        END WHILE */
#line 574 "MB03XP.f"
	}

/*        Multi-shift product QR step. */

#line 578 "MB03XP.f"
	pv2 = ns + 2;
#line 579 "MB03XP.f"
	i__2 = i__ - 1;
#line 579 "MB03XP.f"
	for (k = l; k <= i__2; ++k) {
/* Computing MIN */
#line 580 "MB03XP.f"
	    i__4 = ns + 1, i__5 = i__ - k + 1;
#line 580 "MB03XP.f"
	    nr = min(i__4,i__5);
#line 581 "MB03XP.f"
	    if (k > l) {
#line 581 "MB03XP.f"
		dcopy_(&nr, &a[k + (k - 1) * a_dim1], &c__1, v, &c__1);
#line 581 "MB03XP.f"
	    }
#line 583 "MB03XP.f"
	    dlarfg_(&nr, v, &v[1], &c__1, &tauv);
#line 584 "MB03XP.f"
	    if (k > l) {
#line 585 "MB03XP.f"
		a[k + (k - 1) * a_dim1] = v[0];
#line 586 "MB03XP.f"
		i__4 = i__;
#line 586 "MB03XP.f"
		for (kk = k + 1; kk <= i__4; ++kk) {
#line 587 "MB03XP.f"
		    a[kk + (k - 1) * a_dim1] = 0.;
#line 588 "MB03XP.f"
/* L160: */
#line 588 "MB03XP.f"
		}
#line 589 "MB03XP.f"
	    }

/*           Apply reflector V from the right to B in rows */
/*           I1:min(K+NS,I). */

#line 594 "MB03XP.f"
	    v[0] = 1.;
/* Computing MIN */
#line 595 "MB03XP.f"
	    i__5 = k + ns;
#line 595 "MB03XP.f"
	    i__4 = min(i__5,i__) - i1 + 1;
#line 595 "MB03XP.f"
	    dlarfx_("Right", &i__4, &nr, v, &tauv, &b[i1 + k * b_dim1], ldb, &
		    dwork[1], (ftnlen)5);

/*           Annihilate the introduced nonzeros in the K-th column. */

#line 600 "MB03XP.f"
	    dcopy_(&nr, &b[k + k * b_dim1], &c__1, &v[pv2 - 1], &c__1);
#line 601 "MB03XP.f"
	    dlarfg_(&nr, &v[pv2 - 1], &v[pv2], &c__1, &tauw);
#line 602 "MB03XP.f"
	    b[k + k * b_dim1] = v[pv2 - 1];
#line 603 "MB03XP.f"
	    i__4 = i__;
#line 603 "MB03XP.f"
	    for (kk = k + 1; kk <= i__4; ++kk) {
#line 604 "MB03XP.f"
		b[kk + k * b_dim1] = 0.;
#line 605 "MB03XP.f"
/* L170: */
#line 605 "MB03XP.f"
	    }
#line 606 "MB03XP.f"
	    v[pv2 - 1] = 1.;

/*           Apply reflector W from the left to transform the rows of the */
/*           matrix B in columns K+1:I2. */

#line 611 "MB03XP.f"
	    i__4 = i2 - k;
#line 611 "MB03XP.f"
	    dlarfx_("Left", &nr, &i__4, &v[pv2 - 1], &tauw, &b[k + (k + 1) * 
		    b_dim1], ldb, &dwork[1], (ftnlen)4);

/*           Apply reflector V from the left to transform the rows of the */
/*           matrix A in columns K:I2. */

#line 617 "MB03XP.f"
	    i__4 = i2 - k + 1;
#line 617 "MB03XP.f"
	    dlarfx_("Left", &nr, &i__4, v, &tauv, &a[k + k * a_dim1], lda, &
		    dwork[1], (ftnlen)4);

/*           Apply reflector W from the right to transform the columns of */
/*           the matrix A in rows I1:min(K+NS,I). */

/* Computing MIN */
#line 623 "MB03XP.f"
	    i__5 = k + ns + 1;
#line 623 "MB03XP.f"
	    i__4 = min(i__5,i__) - i1 + 1;
#line 623 "MB03XP.f"
	    dlarfx_("Right", &i__4, &nr, &v[pv2 - 1], &tauw, &a[i1 + k * 
		    a_dim1], lda, &dwork[1], (ftnlen)5);

/*           Accumulate transformations in the matrices Q and Z. */

#line 628 "MB03XP.f"
	    if (wantq) {
#line 628 "MB03XP.f"
		dlarfx_("Right", &nh, &nr, v, &tauv, &q[*ilo + k * q_dim1], 
			ldq, &dwork[1], (ftnlen)5);
#line 628 "MB03XP.f"
	    }
#line 631 "MB03XP.f"
	    if (wantz) {
#line 631 "MB03XP.f"
		dlarfx_("Right", &nh, &nr, &v[pv2 - 1], &tauw, &z__[*ilo + k *
			 z_dim1], ldz, &dwork[1], (ftnlen)5);
#line 631 "MB03XP.f"
	    }
#line 634 "MB03XP.f"
/* L180: */
#line 634 "MB03XP.f"
	}
#line 635 "MB03XP.f"
/* L190: */
#line 635 "MB03XP.f"
    }

/*     Failure to converge. */

#line 639 "MB03XP.f"
    *info = i__;
#line 640 "MB03XP.f"
    return 0;
#line 641 "MB03XP.f"
L200:

/*     Submatrix of order <= MAXB has split off. Use double-shift */
/*     periodic QR algorithm. */

#line 646 "MB03XP.f"
    mb03yd_(&wantt, &wantq, &wantz, n, &l, &i__, ilo, ihi, &a[a_offset], lda, 
	    &b[b_offset], ldb, &q[q_offset], ldq, &z__[z_offset], ldz, &
	    alphar[1], &alphai[1], &beta[1], &dwork[1], ldwork, info);
#line 649 "MB03XP.f"
    if (*info > 0) {
#line 649 "MB03XP.f"
	return 0;
#line 649 "MB03XP.f"
    }
#line 651 "MB03XP.f"
    itn -= its;
#line 652 "MB03XP.f"
    i__ = l - 1;
#line 653 "MB03XP.f"
    goto L90;

#line 655 "MB03XP.f"
L210:
#line 656 "MB03XP.f"
    dwork[1] = (doublereal) max(1,*n);
#line 657 "MB03XP.f"
    return 0;
/* *** Last line of MB03XP *** */
} /* mb03xp_ */

