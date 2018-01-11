#line 1 "MB04DY.f"
/* MB04DY.f -- translated by f2c (version 20100827).
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

#line 1 "MB04DY.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b19 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb04dy_(char *jobscl, integer *n, doublereal *a, integer 
	*lda, doublereal *qg, integer *ldqg, doublereal *d__, doublereal *
	dwork, integer *info, ftnlen jobscl_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal y;
    static integer ihi;
    static doublereal ofl;
    static integer ilo;
    static doublereal ufl, eps, rho, tau, base, anrm;
    static logical none;
    static integer ierr;
    static doublereal gnrm;
    static logical norm;
    static doublereal qnrm;
    static logical symp;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal sfmin, sfmax;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebal_(
	    char *, integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), xerbla_(char *, integer *, ftnlen);
    extern doublereal dlansy_(char *, char *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen, ftnlen);


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

/*     To perform a symplectic scaling on the Hamiltonian matrix */

/*              ( A    G  ) */
/*          H = (       T ),                                          (1) */
/*              ( Q   -A  ) */

/*     i.e., perform either the symplectic scaling transformation */

/*                                   -1 */
/*                 ( A'   G'  )   ( D   0 ) ( A   G  ) ( D  0   ) */
/*          H' <-- (        T ) = (       ) (      T ) (     -1 ),    (2) */
/*                 ( Q'  -A'  )   ( 0   D ) ( Q  -A  ) ( 0  D   ) */

/*     where D is a diagonal scaling matrix, or the symplectic norm */
/*     scaling transformation */

/*                  ( A''   G''  )    1  (   A   G/tau ) */
/*          H'' <-- (          T ) = --- (           T ),             (3) */
/*                  ( Q''  -A''  )   tau ( tau Q   -A  ) */

/*     where tau is a real scalar.  Note that if tau is not equal to 1, */
/*     then (3) is NOT a similarity transformation.  The eigenvalues */
/*     of H are then tau times the eigenvalues of H''. */

/*     For symplectic scaling (2), D is chosen to give the rows and */
/*     columns of A' approximately equal 1-norms and to give Q' and G' */
/*     approximately equal norms.  (See METHOD below for details.) For */
/*     norm scaling, tau = MAX(1, ||A||, ||G||, ||Q||) where ||.|| */
/*     denotes the 1-norm (column sum norm). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBSCL  CHARACTER*1 */
/*             Indicates which scaling strategy is used, as follows: */
/*             = 'S'       :  do the symplectic scaling (2); */
/*             = '1' or 'O':  do the 1-norm scaling (3); */
/*             = 'N'       :  do nothing; set INFO and return. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A, G, and Q.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On input, if JOBSCL <> 'N', the leading N-by-N part of */
/*             this array must contain the upper left block A of the */
/*             Hamiltonian matrix H in (1). */
/*             On output, if JOBSCL <> 'N', the leading N-by-N part of */
/*             this array contains the leading N-by-N part of the scaled */
/*             Hamiltonian matrix H' in (2) or H'' in (3), depending on */
/*             the setting of JOBSCL. */
/*             If JOBSCL = 'N', this array is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,N), if JOBSCL <> 'N'; */
/*             LDA >= 1,        if JOBSCL =  'N'. */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*             (LDQG,N+1) */
/*             On input, if JOBSCL <> 'N', the leading N-by-N lower */
/*             triangular part of this array must contain the lower */
/*             triangle of the lower left symmetric block Q of the */
/*             Hamiltonian matrix H in (1), and the N-by-N upper */
/*             triangular part of the submatrix in the columns 2 to N+1 */
/*             of this array must contain the upper triangle of the upper */
/*             right symmetric block G of H in (1). */
/*             So, if i >= j, then Q(i,j) = Q(j,i) is stored in QG(i,j) */
/*             and G(i,j) = G(j,i) is stored in QG(j,i+1). */
/*             On output, if JOBSCL <> 'N', the leading N-by-N lower */
/*             triangular part of this array contains the lower triangle */
/*             of the lower left symmetric block Q' or Q'', and the */
/*             N-by-N upper triangular part of the submatrix in the */
/*             columns 2 to N+1 of this array contains the upper triangle */
/*             of the upper right symmetric block G' or G'' of the scaled */
/*             Hamiltonian matrix H' in (2) or H'' in (3), depending on */
/*             the setting of JOBSCL. */
/*             If JOBSCL = 'N', this array is not referenced. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG. */
/*             LDQG >= MAX(1,N), if JOBSCL <> 'N'; */
/*             LDQG >= 1,        if JOBSCL =  'N'. */

/*     D       (output) DOUBLE PRECISION array, dimension (nd) */
/*             If JOBSCL = 'S', then nd = N and D contains the diagonal */
/*             elements of the diagonal scaling matrix in (2). */
/*             If JOBSCL = '1' or 'O', then nd = 1 and D(1) is set to tau */
/*             from (3). In this case, no other elements of D are */
/*             referenced. */
/*             If JOBSCL = 'N', this array is not referenced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */
/*             If JOBSCL = 'N', this array is not referenced. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, then the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     1. Symplectic scaling (JOBSCL = 'S'): */

/*     First, LAPACK subroutine DGEBAL is used to equilibrate the 1-norms */
/*     of the rows and columns of A using a diagonal scaling matrix D_A. */
/*     Then, H is similarily transformed by the symplectic diagonal */
/*     matrix D1 = diag(D_A,D_A**(-1)).  Next, the off-diagonal blocks of */
/*     the resulting Hamiltonian matrix are equilibrated in the 1-norm */
/*     using the symplectic diagonal matrix D2 of the form */

/*                 ( I/rho    0   ) */
/*            D2 = (              ) */
/*                 (   0    rho*I ) */

/*     where rho is a real scalar. Thus, in (2), D = D1*D2. */

/*     2. Norm scaling (JOBSCL = '1' or 'O'): */

/*     The norm of the matrices A and G of (1) is reduced by setting */
/*     A := A/tau  and  G := G/(tau**2) where tau is the power of the */
/*     base of the arithmetic closest to MAX(1, ||A||, ||G||, ||Q||) and */
/*     ||.|| denotes the 1-norm. */

/*     REFERENCES */

/*     [1] Benner, P., Byers, R., and Barth, E. */
/*         Fortran 77 Subroutines for Computing the Eigenvalues of */
/*         Hamiltonian Matrices. I: The Square-Reduced Method. */
/*         ACM Trans. Math. Software, 26, 1, pp. 49-77, 2000. */

/*     NUMERICAL ASPECTS */

/*     For symplectic scaling, the complexity of the used algorithms is */
/*     hard to estimate and depends upon how well the rows and columns of */
/*     A in (1) are equilibrated.  In one sweep, each row/column of A is */
/*     scaled once, i.e., the cost of one sweep is N**2 multiplications. */
/*     Usually, 3-6 sweeps are enough to equilibrate the norms of the */
/*     rows and columns of a matrix.  Roundoff errors are possible as */
/*     LAPACK routine DGEBAL does NOT use powers of the machine base for */
/*     scaling. The second stage (equilibrating ||G|| and ||Q||) requires */
/*     N**2 multiplications. */
/*     For norm scaling, 3*N**2 + O(N) multiplications are required and */
/*     NO rounding errors occur as all multiplications are performed with */
/*     powers of the machine base. */

/*     CONTRIBUTOR */

/*     P. Benner, Universitaet Bremen, Germany, and */
/*     R. Byers, University of Kansas, Lawrence, USA. */
/*     Aug. 1998, routine DHABL. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Romania, */
/*     Oct. 1998, SLICOT Library version. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2009. */

/*     KEYWORDS */

/*     Balancing, Hamiltonian matrix, norms, symplectic similarity */
/*     transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */

/*     .. Scalar Arguments .. */
/*    .. */
/*    .. Array Arguments .. */
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

#line 224 "MB04DY.f"
    /* Parameter adjustments */
#line 224 "MB04DY.f"
    a_dim1 = *lda;
#line 224 "MB04DY.f"
    a_offset = 1 + a_dim1;
#line 224 "MB04DY.f"
    a -= a_offset;
#line 224 "MB04DY.f"
    qg_dim1 = *ldqg;
#line 224 "MB04DY.f"
    qg_offset = 1 + qg_dim1;
#line 224 "MB04DY.f"
    qg -= qg_offset;
#line 224 "MB04DY.f"
    --d__;
#line 224 "MB04DY.f"
    --dwork;
#line 224 "MB04DY.f"

#line 224 "MB04DY.f"
    /* Function Body */
#line 224 "MB04DY.f"
    *info = 0;
#line 225 "MB04DY.f"
    symp = lsame_(jobscl, "S", (ftnlen)1, (ftnlen)1);
#line 226 "MB04DY.f"
    norm = lsame_(jobscl, "1", (ftnlen)1, (ftnlen)1) || lsame_(jobscl, "O", (
	    ftnlen)1, (ftnlen)1);
#line 227 "MB04DY.f"
    none = lsame_(jobscl, "N", (ftnlen)1, (ftnlen)1);

#line 229 "MB04DY.f"
    if (! symp && ! norm && ! none) {
#line 230 "MB04DY.f"
	*info = -1;
#line 231 "MB04DY.f"
    } else if (*n < 0) {
#line 232 "MB04DY.f"
	*info = -2;
#line 233 "MB04DY.f"
    } else if (*lda < 1 || ! none && *lda < *n) {
#line 234 "MB04DY.f"
	*info = -4;
#line 235 "MB04DY.f"
    } else if (*ldqg < 1 || ! none && *ldqg < *n) {
#line 236 "MB04DY.f"
	*info = -6;
#line 237 "MB04DY.f"
    }

#line 239 "MB04DY.f"
    if (*info != 0) {
#line 240 "MB04DY.f"
	i__1 = -(*info);
#line 240 "MB04DY.f"
	xerbla_("MB04DY", &i__1, (ftnlen)6);
#line 241 "MB04DY.f"
	return 0;
#line 242 "MB04DY.f"
    }

/*     Quick return if possible. */

#line 246 "MB04DY.f"
    if (*n == 0 || none) {
#line 246 "MB04DY.f"
	return 0;
#line 246 "MB04DY.f"
    }

/*     Set some machine dependant constants. */

#line 251 "MB04DY.f"
    base = dlamch_("Base", (ftnlen)4);
#line 252 "MB04DY.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 253 "MB04DY.f"
    ufl = dlamch_("Safe minimum", (ftnlen)12);
#line 254 "MB04DY.f"
    ofl = 1. / ufl;
#line 255 "MB04DY.f"
    dlabad_(&ufl, &ofl);
#line 256 "MB04DY.f"
    sfmax = eps / base / ufl;
#line 257 "MB04DY.f"
    sfmin = 1. / sfmax;

#line 259 "MB04DY.f"
    if (norm) {

/*        Compute norms. */

#line 263 "MB04DY.f"
	anrm = dlange_("1-norm", n, n, &a[a_offset], lda, &dwork[1], (ftnlen)
		6);
#line 264 "MB04DY.f"
	gnrm = dlansy_("1-norm", "Upper", n, &qg[(qg_dim1 << 1) + 1], ldqg, &
		dwork[1], (ftnlen)6, (ftnlen)5);
#line 265 "MB04DY.f"
	qnrm = dlansy_("1-norm", "Lower", n, &qg[qg_offset], ldqg, &dwork[1], 
		(ftnlen)6, (ftnlen)5);
/* Computing MAX */
#line 266 "MB04DY.f"
	d__1 = max(1.,anrm), d__1 = max(d__1,gnrm);
#line 266 "MB04DY.f"
	y = max(d__1,qnrm);
#line 267 "MB04DY.f"
	tau = 1.;

/*        WHILE ( TAU < Y ) DO */
#line 270 "MB04DY.f"
L10:
#line 271 "MB04DY.f"
	if (tau < y && tau < sqrt(sfmax)) {
#line 272 "MB04DY.f"
	    tau *= base;
#line 273 "MB04DY.f"
	    goto L10;
#line 274 "MB04DY.f"
	}
/*        END WHILE 10 */
#line 276 "MB04DY.f"
	if (tau > 1.) {
#line 277 "MB04DY.f"
	    if ((d__1 = tau / base - y, abs(d__1)) < (d__2 = tau - y, abs(
		    d__2))) {
#line 277 "MB04DY.f"
		tau /= base;
#line 277 "MB04DY.f"
	    }
#line 279 "MB04DY.f"
	    dlascl_("General", &c__0, &c__0, &tau, &c_b19, n, n, &a[a_offset],
		     lda, &ierr, (ftnlen)7);
#line 280 "MB04DY.f"
	    dlascl_("Upper", &c__0, &c__0, &tau, &c_b19, n, n, &qg[(qg_dim1 <<
		     1) + 1], ldqg, &ierr, (ftnlen)5);
#line 282 "MB04DY.f"
	    dlascl_("Upper", &c__0, &c__0, &tau, &c_b19, n, n, &qg[(qg_dim1 <<
		     1) + 1], ldqg, &ierr, (ftnlen)5);
#line 284 "MB04DY.f"
	}

#line 286 "MB04DY.f"
	d__[1] = tau;

#line 288 "MB04DY.f"
    } else {
#line 289 "MB04DY.f"
	dgebal_("Scale", n, &a[a_offset], lda, &ilo, &ihi, &d__[1], &ierr, (
		ftnlen)5);

#line 291 "MB04DY.f"
	i__1 = *n;
#line 291 "MB04DY.f"
	for (j = 1; j <= i__1; ++j) {

#line 293 "MB04DY.f"
	    i__2 = *n;
#line 293 "MB04DY.f"
	    for (i__ = j; i__ <= i__2; ++i__) {
#line 294 "MB04DY.f"
		qg[i__ + j * qg_dim1] = qg[i__ + j * qg_dim1] * d__[j] * d__[
			i__];
#line 295 "MB04DY.f"
/* L20: */
#line 295 "MB04DY.f"
	    }

#line 297 "MB04DY.f"
/* L30: */
#line 297 "MB04DY.f"
	}

#line 299 "MB04DY.f"
	i__1 = *n + 1;
#line 299 "MB04DY.f"
	for (j = 2; j <= i__1; ++j) {

#line 301 "MB04DY.f"
	    i__2 = j - 1;
#line 301 "MB04DY.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 302 "MB04DY.f"
		qg[i__ + j * qg_dim1] = qg[i__ + j * qg_dim1] / d__[j - 1] / 
			d__[i__];
#line 303 "MB04DY.f"
/* L40: */
#line 303 "MB04DY.f"
	    }

#line 305 "MB04DY.f"
/* L50: */
#line 305 "MB04DY.f"
	}

#line 307 "MB04DY.f"
	gnrm = dlansy_("1-norm", "Upper", n, &qg[(qg_dim1 << 1) + 1], ldqg, &
		dwork[1], (ftnlen)6, (ftnlen)5);
#line 308 "MB04DY.f"
	qnrm = dlansy_("1-norm", "Lower", n, &qg[qg_offset], ldqg, &dwork[1], 
		(ftnlen)6, (ftnlen)5);
#line 309 "MB04DY.f"
	if (gnrm == 0.) {
#line 310 "MB04DY.f"
	    if (qnrm == 0.) {
#line 311 "MB04DY.f"
		rho = 1.;
#line 312 "MB04DY.f"
	    } else {
#line 313 "MB04DY.f"
		rho = sfmax;
#line 314 "MB04DY.f"
	    }
#line 315 "MB04DY.f"
	} else if (qnrm == 0.) {
#line 316 "MB04DY.f"
	    rho = sfmin;
#line 317 "MB04DY.f"
	} else {
#line 318 "MB04DY.f"
	    rho = sqrt(qnrm) / sqrt(gnrm);
#line 319 "MB04DY.f"
	}

#line 321 "MB04DY.f"
	dlascl_("Lower", &c__0, &c__0, &rho, &c_b19, n, n, &qg[qg_offset], 
		ldqg, &ierr, (ftnlen)5);
#line 322 "MB04DY.f"
	dlascl_("Upper", &c__0, &c__0, &c_b19, &rho, n, n, &qg[(qg_dim1 << 1) 
		+ 1], ldqg, &ierr, (ftnlen)5);
#line 324 "MB04DY.f"
	d__1 = sqrt(rho);
#line 324 "MB04DY.f"
	drscl_(n, &d__1, &d__[1], &c__1);
#line 325 "MB04DY.f"
    }

#line 327 "MB04DY.f"
    return 0;
/*     *** Last line of MB04DY *** */
} /* mb04dy_ */

