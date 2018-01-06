#line 1 "MB04ZD.f"
/* MB04ZD.f -- translated by f2c (version 20100827).
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

#line 1 "MB04ZD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b13 = -1.;
static doublereal c_b15 = 0.;
static doublereal c_b18 = 1.;
static integer c__2 = 2;
static integer c__0 = 0;

/* Subroutine */ int mb04zd_(char *compu, integer *n, doublereal *a, integer *
	lda, doublereal *qg, integer *ldqg, doublereal *u, integer *ldu, 
	doublereal *dwork, integer *info, ftnlen compu_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, u_dim1, u_offset, i__1, 
	    i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal t[4]	/* was [2][2] */, x, y, tau;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal sine;
    static logical form;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dsyr2_(char 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen);
    static logical accum;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), daxpy_(integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *)
	    ;
    static doublereal dummy[1];
    extern /* Subroutine */ int dsymv_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dlarfg_(integer *, doublereal *,
	     doublereal *, integer *, doublereal *), dlartg_(doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), xerbla_(
	    char *, integer *, ftnlen);
    static doublereal cosine;
    extern /* Subroutine */ int dlarfx_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen);
    static logical forget;


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

/*     To transform a Hamiltonian matrix */

/*               ( A   G  ) */
/*           H = (      T )                                           (1) */
/*               ( Q  -A  ) */

/*     into a square-reduced Hamiltonian matrix */

/*                ( A'  G'  ) */
/*           H' = (       T )                                         (2) */
/*                ( Q' -A'  ) */
/*                                                                 T */
/*     by an orthogonal symplectic similarity transformation H' = U H U, */
/*     where */
/*               (  U1   U2 ) */
/*           U = (          ).                                        (3) */
/*               ( -U2   U1 ) */
/*                                                              T */
/*     The square-reduced Hamiltonian matrix satisfies Q'A' - A' Q' = 0, */
/*     and */

/*           2       T     2     ( A''   G''  ) */
/*         H'  :=  (U  H U)   =  (          T ). */
/*                               ( 0     A''  ) */

/*     In addition, A'' is upper Hessenberg and G'' is skew symmetric. */
/*     The square roots of the eigenvalues of A'' = A'*A' + G'*Q' are the */
/*     eigenvalues of H. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     COMPU   CHARACTER*1 */
/*             Indicates whether the orthogonal symplectic similarity */
/*             transformation matrix U in (3) is returned or */
/*             accumulated into an orthogonal symplectic matrix, or if */
/*             the transformation matrix is not required, as follows: */
/*             = 'N':         U is not required; */
/*             = 'I' or 'F':  on entry, U need not be set; */
/*                            on exit, U contains the orthogonal */
/*                            symplectic matrix U from (3); */
/*             = 'V' or 'A':  the orthogonal symplectic similarity */
/*                            transformations are accumulated into U; */
/*                            on input, U must contain an orthogonal */
/*                            symplectic matrix S; */
/*                            on exit, U contains S*U with U from (3). */
/*             See the description of U below for details. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, G, and Q.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On input, the leading N-by-N part of this array must */
/*             contain the upper left block A of the Hamiltonian matrix H */
/*             in (1). */
/*             On output, the leading N-by-N part of this array contains */
/*             the upper left block A' of the square-reduced Hamiltonian */
/*             matrix H' in (2). */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*             (LDQG,N+1) */
/*             On input, the leading N-by-N lower triangular part of this */
/*             array must contain the lower triangle of the lower left */
/*             symmetric block Q of the Hamiltonian matrix H in (1), and */
/*             the N-by-N upper triangular part of the submatrix in the */
/*             columns 2 to N+1 of this array must contain the upper */
/*             triangle of the upper right symmetric block G of H in (1). */
/*             So, if i >= j, then Q(i,j) = Q(j,i) is stored in QG(i,j) */
/*             and G(i,j) = G(j,i) is stored in QG(j,i+1). */
/*             On output, the leading N-by-N lower triangular part of */
/*             this array contains the lower triangle of the lower left */
/*             symmetric block Q', and the N-by-N upper triangular part */
/*             of the submatrix in the columns 2 to N+1 of this array */
/*             contains the upper triangle of the upper right symmetric */
/*             block G' of the square-reduced Hamiltonian matrix H' */
/*             in (2). */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= MAX(1,N). */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,2*N) */
/*             If COMPU = 'N', then this array is not referenced. */
/*             If COMPU = 'I' or 'F', then the input contents of this */
/*             array are not specified.  On output, the leading */
/*             N-by-(2*N) part of this array contains the first N rows */
/*             of the orthogonal symplectic matrix U in (3). */
/*             If COMPU = 'V' or 'A', then, on input, the leading */
/*             N-by-(2*N) part of this array must contain the first N */
/*             rows of an orthogonal symplectic matrix S. On output, the */
/*             leading N-by-(2*N) part of this array contains the first N */
/*             rows of the product S*U where U is the orthogonal */
/*             symplectic matrix from (3). */
/*             The storage scheme implied by (3) is used for orthogonal */
/*             symplectic matrices, i.e., only the first N rows are */
/*             stored, as they contain all relevant information. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= MAX(1,N), if COMPU <> 'N'; */
/*             LDU >= 1,        if COMPU =  'N'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (2*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, then the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The Hamiltonian matrix H is transformed into a square-reduced */
/*     Hamiltonian matrix H' using the implicit version of Van Loan's */
/*     method as proposed in [1,2,3]. */

/*     REFERENCES */

/*     [1] Van Loan, C. F. */
/*         A Symplectic Method for Approximating All the Eigenvalues of */
/*         a Hamiltonian Matrix. */
/*         Linear Algebra and its Applications, 61, pp. 233-251, 1984. */

/*     [2] Byers, R. */
/*         Hamiltonian and Symplectic Algorithms for the Algebraic */
/*         Riccati Equation. */
/*         Ph. D. Thesis, Cornell University, Ithaca, NY, January 1983. */

/*     [3] Benner, P., Byers, R., and Barth, E. */
/*         Fortran 77 Subroutines for Computing the Eigenvalues of */
/*         Hamiltonian Matrices. I: The Square-Reduced Method. */
/*         ACM Trans. Math. Software, 26, 1, pp. 49-77, 2000. */

/*     NUMERICAL ASPECTS */

/*     This algorithm requires approximately 20*N**3 flops for */
/*     transforming H into square-reduced form. If the transformations */
/*     are required, this adds another 8*N**3 flops. The method is */
/*     strongly backward stable in the sense that if H' and U are the */
/*     computed square-reduced Hamiltonian and computed orthogonal */
/*     symplectic similarity transformation, then there is an orthogonal */
/*     symplectic matrix T and a Hamiltonian matrix M such that */

/*                  H T  =  T M */

/*        || T - U ||   <=  c1 * eps */

/*        || H' - M ||  <=  c2 * eps * || H || */

/*     where c1, c2 are modest constants depending on the dimension N and */
/*     eps is the machine precision. */

/*     Eigenvalues computed by explicitly forming the upper Hessenberg */
/*     matrix  A'' = A'A' + G'Q', with A', G', and Q' as in (2), and */
/*     applying the Hessenberg QR iteration to A'' are exactly */
/*     eigenvalues of a perturbed Hamiltonian matrix H + E,  where */

/*        || E ||  <=  c3 * sqrt(eps) * || H ||, */

/*     and c3 is a modest constant depending on the dimension N and eps */
/*     is the machine precision.  Moreover, if the norm of H and an */
/*     eigenvalue lambda are of roughly the same magnitude, the computed */
/*     eigenvalue is essentially as accurate as the computed eigenvalue */
/*     from traditional methods.  See [1] or [2]. */

/*     CONTRIBUTOR */

/*     P. Benner, Universitaet Bremen, Germany, */
/*     R. Byers, University of Kansas, Lawrence, USA, and */
/*     E. Barth, Kalamazoo College, Kalamazoo, USA, */
/*     Aug. 1998, routine DHASRD. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998, SLICOT Library version. */

/*     REVISIONS */

/*     May 2001, A. Varga, German Aeropsce Center, DLR Oberpfaffenhofen. */
/*     May 2009, V. Sima, Research Institute for Informatics, Bucharest. */

/*     KEYWORDS */

/*     Orthogonal transformation, (square-reduced) Hamiltonian matrix, */
/*     symplectic similarity transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */


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

#line 252 "MB04ZD.f"
    /* Parameter adjustments */
#line 252 "MB04ZD.f"
    a_dim1 = *lda;
#line 252 "MB04ZD.f"
    a_offset = 1 + a_dim1;
#line 252 "MB04ZD.f"
    a -= a_offset;
#line 252 "MB04ZD.f"
    qg_dim1 = *ldqg;
#line 252 "MB04ZD.f"
    qg_offset = 1 + qg_dim1;
#line 252 "MB04ZD.f"
    qg -= qg_offset;
#line 252 "MB04ZD.f"
    u_dim1 = *ldu;
#line 252 "MB04ZD.f"
    u_offset = 1 + u_dim1;
#line 252 "MB04ZD.f"
    u -= u_offset;
#line 252 "MB04ZD.f"
    --dwork;
#line 252 "MB04ZD.f"

#line 252 "MB04ZD.f"
    /* Function Body */
#line 252 "MB04ZD.f"
    *info = 0;
#line 253 "MB04ZD.f"
    accum = lsame_(compu, "A", (ftnlen)1, (ftnlen)1) || lsame_(compu, "V", (
	    ftnlen)1, (ftnlen)1);
#line 254 "MB04ZD.f"
    form = lsame_(compu, "F", (ftnlen)1, (ftnlen)1) || lsame_(compu, "I", (
	    ftnlen)1, (ftnlen)1);
#line 255 "MB04ZD.f"
    forget = lsame_(compu, "N", (ftnlen)1, (ftnlen)1);

#line 257 "MB04ZD.f"
    if (! accum && ! form && ! forget) {
#line 258 "MB04ZD.f"
	*info = -1;
#line 259 "MB04ZD.f"
    } else if (*n < 0) {
#line 260 "MB04ZD.f"
	*info = -2;
#line 261 "MB04ZD.f"
    } else if (*lda < max(1,*n)) {
#line 262 "MB04ZD.f"
	*info = -4;
#line 263 "MB04ZD.f"
    } else if (*ldqg < max(1,*n)) {
#line 264 "MB04ZD.f"
	*info = -6;
#line 265 "MB04ZD.f"
    } else if (*ldu < 1 || ! forget && *ldu < max(1,*n)) {
#line 267 "MB04ZD.f"
	*info = -8;
#line 268 "MB04ZD.f"
    }

#line 270 "MB04ZD.f"
    if (*info != 0) {
#line 271 "MB04ZD.f"
	i__1 = -(*info);
#line 271 "MB04ZD.f"
	xerbla_("MB04ZD", &i__1, (ftnlen)6);
#line 272 "MB04ZD.f"
	return 0;
#line 273 "MB04ZD.f"
    }

/*     Quick return if possible. */

#line 277 "MB04ZD.f"
    if (*n == 0) {
#line 277 "MB04ZD.f"
	return 0;
#line 277 "MB04ZD.f"
    }

/*     Transform to square-reduced form. */

#line 282 "MB04ZD.f"
    i__1 = *n - 1;
#line 282 "MB04ZD.f"
    for (j = 1; j <= i__1; ++j) {
/*                         T */
/*        DWORK <- (Q*A - A *Q)(J+1:N,J). */

#line 286 "MB04ZD.f"
	i__2 = j - 1;
#line 286 "MB04ZD.f"
	dcopy_(&i__2, &qg[j + qg_dim1], ldqg, &dwork[*n + 1], &c__1);
#line 287 "MB04ZD.f"
	i__2 = *n - j + 1;
#line 287 "MB04ZD.f"
	dcopy_(&i__2, &qg[j + j * qg_dim1], &c__1, &dwork[*n + j], &c__1);
#line 288 "MB04ZD.f"
	i__2 = *n - j;
#line 288 "MB04ZD.f"
	dgemv_("Transpose", n, &i__2, &c_b13, &a[(j + 1) * a_dim1 + 1], lda, &
		dwork[*n + 1], &c__1, &c_b15, &dwork[j + 1], &c__1, (ftnlen)9)
		;
#line 290 "MB04ZD.f"
	i__2 = *n - j;
#line 290 "MB04ZD.f"
	dgemv_("NoTranspose", &i__2, &j, &c_b18, &qg[j + 1 + qg_dim1], ldqg, &
		a[j * a_dim1 + 1], &c__1, &c_b18, &dwork[j + 1], &c__1, (
		ftnlen)11);
#line 292 "MB04ZD.f"
	i__2 = *n - j;
#line 292 "MB04ZD.f"
	dsymv_("Lower", &i__2, &c_b18, &qg[j + 1 + (j + 1) * qg_dim1], ldqg, &
		a[j + 1 + j * a_dim1], &c__1, &c_b18, &dwork[j + 1], &c__1, (
		ftnlen)5);

/*        Symplectic reflection to zero (H*H)((N+J+2):2N,J). */

#line 297 "MB04ZD.f"
	i__2 = *n - j;
#line 297 "MB04ZD.f"
	dlarfg_(&i__2, &dwork[j + 1], &dwork[j + 2], &c__1, &tau);
#line 298 "MB04ZD.f"
	y = dwork[j + 1];
#line 299 "MB04ZD.f"
	dwork[j + 1] = 1.;

#line 301 "MB04ZD.f"
	i__2 = *n - j;
#line 301 "MB04ZD.f"
	dlarfx_("Left", &i__2, n, &dwork[j + 1], &tau, &a[j + 1 + a_dim1], 
		lda, &dwork[*n + 1], (ftnlen)4);
#line 303 "MB04ZD.f"
	i__2 = *n - j;
#line 303 "MB04ZD.f"
	dlarfx_("Right", n, &i__2, &dwork[j + 1], &tau, &a[(j + 1) * a_dim1 + 
		1], lda, &dwork[*n + 1], (ftnlen)5);

#line 306 "MB04ZD.f"
	i__2 = *n - j;
#line 306 "MB04ZD.f"
	dlarfx_("Left", &i__2, &j, &dwork[j + 1], &tau, &qg[j + 1 + qg_dim1], 
		ldqg, &dwork[*n + 1], (ftnlen)4);
#line 308 "MB04ZD.f"
	i__2 = *n - j;
#line 308 "MB04ZD.f"
	dsymv_("Lower", &i__2, &tau, &qg[j + 1 + (j + 1) * qg_dim1], ldqg, &
		dwork[j + 1], &c__1, &c_b15, &dwork[*n + j + 1], &c__1, (
		ftnlen)5);
#line 310 "MB04ZD.f"
	i__2 = *n - j;
#line 310 "MB04ZD.f"
	i__3 = *n - j;
#line 310 "MB04ZD.f"
	d__1 = -tau * ddot_(&i__3, &dwork[*n + j + 1], &c__1, &dwork[j + 1], &
		c__1) / 2.;
#line 310 "MB04ZD.f"
	daxpy_(&i__2, &d__1, &dwork[j + 1], &c__1, &dwork[*n + j + 1], &c__1);
#line 312 "MB04ZD.f"
	i__2 = *n - j;
#line 312 "MB04ZD.f"
	dsyr2_("Lower", &i__2, &c_b13, &dwork[j + 1], &c__1, &dwork[*n + j + 
		1], &c__1, &qg[j + 1 + (j + 1) * qg_dim1], ldqg, (ftnlen)5);

#line 315 "MB04ZD.f"
	i__2 = *n - j;
#line 315 "MB04ZD.f"
	dlarfx_("Right", &j, &i__2, &dwork[j + 1], &tau, &qg[(j + 2) * 
		qg_dim1 + 1], ldqg, &dwork[*n + 1], (ftnlen)5);
#line 317 "MB04ZD.f"
	i__2 = *n - j;
#line 317 "MB04ZD.f"
	dsymv_("Upper", &i__2, &tau, &qg[j + 1 + (j + 2) * qg_dim1], ldqg, &
		dwork[j + 1], &c__1, &c_b15, &dwork[*n + j + 1], &c__1, (
		ftnlen)5);
#line 319 "MB04ZD.f"
	i__2 = *n - j;
#line 319 "MB04ZD.f"
	i__3 = *n - j;
#line 319 "MB04ZD.f"
	d__1 = -tau * ddot_(&i__3, &dwork[*n + j + 1], &c__1, &dwork[j + 1], &
		c__1) / 2.;
#line 319 "MB04ZD.f"
	daxpy_(&i__2, &d__1, &dwork[j + 1], &c__1, &dwork[*n + j + 1], &c__1);
#line 321 "MB04ZD.f"
	i__2 = *n - j;
#line 321 "MB04ZD.f"
	dsyr2_("Upper", &i__2, &c_b13, &dwork[j + 1], &c__1, &dwork[*n + j + 
		1], &c__1, &qg[j + 1 + (j + 2) * qg_dim1], ldqg, (ftnlen)5);

#line 324 "MB04ZD.f"
	if (form) {

/*           Save reflection. */

#line 328 "MB04ZD.f"
	    i__2 = *n - j;
#line 328 "MB04ZD.f"
	    dcopy_(&i__2, &dwork[j + 1], &c__1, &u[j + 1 + j * u_dim1], &c__1)
		    ;
#line 329 "MB04ZD.f"
	    u[j + 1 + j * u_dim1] = tau;

#line 331 "MB04ZD.f"
	} else if (accum) {

/*           Accumulate reflection. */

#line 335 "MB04ZD.f"
	    i__2 = *n - j;
#line 335 "MB04ZD.f"
	    dlarfx_("Right", n, &i__2, &dwork[j + 1], &tau, &u[(j + 1) * 
		    u_dim1 + 1], ldu, &dwork[*n + 1], (ftnlen)5);
#line 337 "MB04ZD.f"
	    i__2 = *n - j;
#line 337 "MB04ZD.f"
	    dlarfx_("Right", n, &i__2, &dwork[j + 1], &tau, &u[(*n + j + 1) * 
		    u_dim1 + 1], ldu, &dwork[*n + 1], (ftnlen)5);
#line 339 "MB04ZD.f"
	}

/*        (X,Y) := ((J+1,J),(N+J+1,J)) component of H*H. */

#line 343 "MB04ZD.f"
	i__2 = *n - j;
#line 343 "MB04ZD.f"
	x = ddot_(&j, &qg[(j + 2) * qg_dim1 + 1], &c__1, &qg[j + qg_dim1], 
		ldqg) + ddot_(&i__2, &qg[j + 1 + (j + 2) * qg_dim1], ldqg, &
		qg[j + 1 + j * qg_dim1], &c__1) + ddot_(n, &a[j + 1 + a_dim1],
		 lda, &a[j * a_dim1 + 1], &c__1);

/*        Symplectic rotation to zero (H*H)(N+J+1,J). */

#line 349 "MB04ZD.f"
	dlartg_(&x, &y, &cosine, &sine, &temp);

#line 351 "MB04ZD.f"
	drot_(&j, &a[j + 1 + a_dim1], lda, &qg[j + 1 + qg_dim1], ldqg, &
		cosine, &sine);
#line 352 "MB04ZD.f"
	drot_(&j, &a[(j + 1) * a_dim1 + 1], &c__1, &qg[(j + 2) * qg_dim1 + 1],
		 &c__1, &cosine, &sine);
#line 353 "MB04ZD.f"
	if (j < *n - 1) {
#line 354 "MB04ZD.f"
	    i__2 = *n - j - 1;
#line 354 "MB04ZD.f"
	    drot_(&i__2, &a[j + 1 + (j + 2) * a_dim1], lda, &qg[j + 2 + (j + 
		    1) * qg_dim1], &c__1, &cosine, &sine);
#line 356 "MB04ZD.f"
	    i__2 = *n - j - 1;
#line 356 "MB04ZD.f"
	    drot_(&i__2, &a[j + 2 + (j + 1) * a_dim1], &c__1, &qg[j + 1 + (j 
		    + 3) * qg_dim1], ldqg, &cosine, &sine);
#line 358 "MB04ZD.f"
	}

#line 360 "MB04ZD.f"
	t[0] = a[j + 1 + (j + 1) * a_dim1];
#line 361 "MB04ZD.f"
	t[2] = qg[j + 1 + (j + 2) * qg_dim1];
#line 362 "MB04ZD.f"
	t[1] = qg[j + 1 + (j + 1) * qg_dim1];
#line 363 "MB04ZD.f"
	t[3] = -t[0];
#line 364 "MB04ZD.f"
	drot_(&c__2, t, &c__1, &t[2], &c__1, &cosine, &sine);
#line 365 "MB04ZD.f"
	drot_(&c__2, t, &c__2, &t[1], &c__2, &cosine, &sine);
#line 366 "MB04ZD.f"
	a[j + 1 + (j + 1) * a_dim1] = t[0];
#line 367 "MB04ZD.f"
	qg[j + 1 + (j + 2) * qg_dim1] = t[2];
#line 368 "MB04ZD.f"
	qg[j + 1 + (j + 1) * qg_dim1] = t[1];

#line 370 "MB04ZD.f"
	if (form) {

/*           Save rotation. */

#line 374 "MB04ZD.f"
	    u[j + j * u_dim1] = cosine;
#line 375 "MB04ZD.f"
	    u[j + (*n + j) * u_dim1] = sine;

#line 377 "MB04ZD.f"
	} else if (accum) {

/*           Accumulate rotation. */

#line 381 "MB04ZD.f"
	    drot_(n, &u[(j + 1) * u_dim1 + 1], &c__1, &u[(*n + j + 1) * 
		    u_dim1 + 1], &c__1, &cosine, &sine);
#line 382 "MB04ZD.f"
	}

/*        DWORK := (A*A  + G*Q)(J+1:N,J). */

#line 386 "MB04ZD.f"
	i__2 = *n - j;
#line 386 "MB04ZD.f"
	dgemv_("NoTranspose", &i__2, n, &c_b18, &a[j + 1 + a_dim1], lda, &a[j 
		* a_dim1 + 1], &c__1, &c_b15, &dwork[j + 1], &c__1, (ftnlen)
		11);
#line 388 "MB04ZD.f"
	i__2 = *n - j;
#line 388 "MB04ZD.f"
	dgemv_("Transpose", &j, &i__2, &c_b18, &qg[(j + 2) * qg_dim1 + 1], 
		ldqg, &qg[j + qg_dim1], ldqg, &c_b18, &dwork[j + 1], &c__1, (
		ftnlen)9);
#line 390 "MB04ZD.f"
	i__2 = *n - j;
#line 390 "MB04ZD.f"
	dsymv_("Upper", &i__2, &c_b18, &qg[j + 1 + (j + 2) * qg_dim1], ldqg, &
		qg[j + 1 + j * qg_dim1], &c__1, &c_b18, &dwork[j + 1], &c__1, 
		(ftnlen)5);

/*        Symplectic reflection to zero (H*H)(J+2:N,J). */

#line 395 "MB04ZD.f"
	i__2 = *n - j;
#line 395 "MB04ZD.f"
	dlarfg_(&i__2, &dwork[j + 1], &dwork[j + 2], &c__1, &tau);
#line 396 "MB04ZD.f"
	dwork[j + 1] = 1.;

#line 398 "MB04ZD.f"
	i__2 = *n - j;
#line 398 "MB04ZD.f"
	dlarfx_("Left", &i__2, n, &dwork[j + 1], &tau, &a[j + 1 + a_dim1], 
		lda, &dwork[*n + 1], (ftnlen)4);
#line 400 "MB04ZD.f"
	i__2 = *n - j;
#line 400 "MB04ZD.f"
	dlarfx_("Right", n, &i__2, &dwork[j + 1], &tau, &a[(j + 1) * a_dim1 + 
		1], lda, &dwork[*n + 1], (ftnlen)5);

#line 403 "MB04ZD.f"
	i__2 = *n - j;
#line 403 "MB04ZD.f"
	dlarfx_("Left", &i__2, &j, &dwork[j + 1], &tau, &qg[j + 1 + qg_dim1], 
		ldqg, &dwork[*n + 1], (ftnlen)4);
#line 405 "MB04ZD.f"
	i__2 = *n - j;
#line 405 "MB04ZD.f"
	dsymv_("Lower", &i__2, &tau, &qg[j + 1 + (j + 1) * qg_dim1], ldqg, &
		dwork[j + 1], &c__1, &c_b15, &dwork[*n + j + 1], &c__1, (
		ftnlen)5);
#line 407 "MB04ZD.f"
	i__2 = *n - j;
#line 407 "MB04ZD.f"
	i__3 = *n - j;
#line 407 "MB04ZD.f"
	d__1 = -tau * ddot_(&i__3, &dwork[*n + j + 1], &c__1, &dwork[j + 1], &
		c__1) / 2.;
#line 407 "MB04ZD.f"
	daxpy_(&i__2, &d__1, &dwork[j + 1], &c__1, &dwork[*n + j + 1], &c__1);
#line 409 "MB04ZD.f"
	i__2 = *n - j;
#line 409 "MB04ZD.f"
	dsyr2_("Lower", &i__2, &c_b13, &dwork[j + 1], &c__1, &dwork[*n + j + 
		1], &c__1, &qg[j + 1 + (j + 1) * qg_dim1], ldqg, (ftnlen)5);

#line 412 "MB04ZD.f"
	i__2 = *n - j;
#line 412 "MB04ZD.f"
	dlarfx_("Right", &j, &i__2, &dwork[j + 1], &tau, &qg[(j + 2) * 
		qg_dim1 + 1], ldqg, &dwork[*n + 1], (ftnlen)5);
#line 414 "MB04ZD.f"
	i__2 = *n - j;
#line 414 "MB04ZD.f"
	dsymv_("Upper", &i__2, &tau, &qg[j + 1 + (j + 2) * qg_dim1], ldqg, &
		dwork[j + 1], &c__1, &c_b15, &dwork[*n + j + 1], &c__1, (
		ftnlen)5);
#line 416 "MB04ZD.f"
	i__2 = *n - j;
#line 416 "MB04ZD.f"
	i__3 = *n - j;
#line 416 "MB04ZD.f"
	d__1 = -tau * ddot_(&i__3, &dwork[*n + j + 1], &c__1, &dwork[j + 1], &
		c__1) / 2.;
#line 416 "MB04ZD.f"
	daxpy_(&i__2, &d__1, &dwork[j + 1], &c__1, &dwork[*n + j + 1], &c__1);
#line 418 "MB04ZD.f"
	i__2 = *n - j;
#line 418 "MB04ZD.f"
	dsyr2_("Upper", &i__2, &c_b13, &dwork[j + 1], &c__1, &dwork[*n + j + 
		1], &c__1, &qg[j + 1 + (j + 2) * qg_dim1], ldqg, (ftnlen)5);

#line 421 "MB04ZD.f"
	if (form) {

/*           Save reflection. */

#line 425 "MB04ZD.f"
	    i__2 = *n - j;
#line 425 "MB04ZD.f"
	    dcopy_(&i__2, &dwork[j + 1], &c__1, &u[j + 1 + (*n + j) * u_dim1],
		     &c__1);
#line 426 "MB04ZD.f"
	    u[j + 1 + (*n + j) * u_dim1] = tau;

#line 428 "MB04ZD.f"
	} else if (accum) {

/*           Accumulate reflection. */

#line 432 "MB04ZD.f"
	    i__2 = *n - j;
#line 432 "MB04ZD.f"
	    dlarfx_("Right", n, &i__2, &dwork[j + 1], &tau, &u[(j + 1) * 
		    u_dim1 + 1], ldu, &dwork[*n + 1], (ftnlen)5);
#line 434 "MB04ZD.f"
	    i__2 = *n - j;
#line 434 "MB04ZD.f"
	    dlarfx_("Right", n, &i__2, &dwork[j + 1], &tau, &u[(*n + j + 1) * 
		    u_dim1 + 1], ldu, &dwork[*n + 1], (ftnlen)5);
#line 436 "MB04ZD.f"
	}

#line 438 "MB04ZD.f"
/* L10: */
#line 438 "MB04ZD.f"
    }

#line 440 "MB04ZD.f"
    if (form) {
#line 441 "MB04ZD.f"
	dummy[0] = 0.;

/*        Form S by accumulating transformations. */

#line 445 "MB04ZD.f"
	for (j = *n - 1; j >= 1; --j) {

/*           Initialize (J+1)st column of S. */

#line 449 "MB04ZD.f"
	    dcopy_(n, dummy, &c__0, &u[(j + 1) * u_dim1 + 1], &c__1);
#line 450 "MB04ZD.f"
	    u[j + 1 + (j + 1) * u_dim1] = 1.;
#line 451 "MB04ZD.f"
	    dcopy_(n, dummy, &c__0, &u[(*n + j + 1) * u_dim1 + 1], &c__1);

/*           Second reflection. */

#line 455 "MB04ZD.f"
	    tau = u[j + 1 + (*n + j) * u_dim1];
#line 456 "MB04ZD.f"
	    u[j + 1 + (*n + j) * u_dim1] = 1.;
#line 457 "MB04ZD.f"
	    i__1 = *n - j;
#line 457 "MB04ZD.f"
	    i__2 = *n - j;
#line 457 "MB04ZD.f"
	    dlarfx_("Left", &i__1, &i__2, &u[j + 1 + (*n + j) * u_dim1], &tau,
		     &u[j + 1 + (j + 1) * u_dim1], ldu, &dwork[*n + 1], (
		    ftnlen)4);
#line 459 "MB04ZD.f"
	    i__1 = *n - j;
#line 459 "MB04ZD.f"
	    i__2 = *n - j;
#line 459 "MB04ZD.f"
	    dlarfx_("Left", &i__1, &i__2, &u[j + 1 + (*n + j) * u_dim1], &tau,
		     &u[j + 1 + (*n + j + 1) * u_dim1], ldu, &dwork[*n + 1], (
		    ftnlen)4);

/*           Rotation. */

#line 464 "MB04ZD.f"
	    i__1 = *n - j;
#line 464 "MB04ZD.f"
	    drot_(&i__1, &u[j + 1 + (j + 1) * u_dim1], ldu, &u[j + 1 + (*n + 
		    j + 1) * u_dim1], ldu, &u[j + j * u_dim1], &u[j + (*n + j)
		     * u_dim1]);

/*           First reflection. */

#line 469 "MB04ZD.f"
	    tau = u[j + 1 + j * u_dim1];
#line 470 "MB04ZD.f"
	    u[j + 1 + j * u_dim1] = 1.;
#line 471 "MB04ZD.f"
	    i__1 = *n - j;
#line 471 "MB04ZD.f"
	    i__2 = *n - j;
#line 471 "MB04ZD.f"
	    dlarfx_("Left", &i__1, &i__2, &u[j + 1 + j * u_dim1], &tau, &u[j 
		    + 1 + (j + 1) * u_dim1], ldu, &dwork[*n + 1], (ftnlen)4);
#line 473 "MB04ZD.f"
	    i__1 = *n - j;
#line 473 "MB04ZD.f"
	    i__2 = *n - j;
#line 473 "MB04ZD.f"
	    dlarfx_("Left", &i__1, &i__2, &u[j + 1 + j * u_dim1], &tau, &u[j 
		    + 1 + (*n + j + 1) * u_dim1], ldu, &dwork[*n + 1], (
		    ftnlen)4);
#line 475 "MB04ZD.f"
/* L20: */
#line 475 "MB04ZD.f"
	}

/*        The first column is the first column of identity. */

#line 479 "MB04ZD.f"
	dcopy_(n, dummy, &c__0, &u[u_offset], &c__1);
#line 480 "MB04ZD.f"
	u[u_dim1 + 1] = 1.;
#line 481 "MB04ZD.f"
	dcopy_(n, dummy, &c__0, &u[(*n + 1) * u_dim1 + 1], &c__1);
#line 482 "MB04ZD.f"
    }

#line 484 "MB04ZD.f"
    return 0;
/*     *** Last line of MB04ZD *** */
} /* mb04zd_ */

