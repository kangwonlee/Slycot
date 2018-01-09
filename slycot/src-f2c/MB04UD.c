#line 1 "MB04UD.f"
/* MB04UD.f -- translated by f2c (version 20100827).
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

#line 1 "MB04UD.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static doublereal c_b11 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb04ud_(char *jobq, char *jobz, integer *m, integer *n, 
	doublereal *a, integer *lda, doublereal *e, integer *lde, doublereal *
	q, integer *ldq, doublereal *z__, integer *ldz, integer *ranke, 
	integer *istair, doublereal *tol, doublereal *dwork, integer *info, 
	ftnlen jobq_len, ftnlen jobz_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, l, lk, km1, nr1, mnk;
    static doublereal emx, tau;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal toler;
    static logical lzero;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical ljobqi, ljobzi, updatq;
    static doublereal emxnrm;
    static logical updatz;


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

/*     To compute orthogonal transformations Q and Z such that the */
/*     transformed pencil Q'(sE-A)Z has the E matrix in column echelon */
/*     form, where E and A are M-by-N matrices. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBQ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Q the unitary row permutations, as follows: */
/*             = 'N':  Do not form Q; */
/*             = 'I':  Q is initialized to the unit matrix and the */
/*                     unitary row permutation matrix Q is returned; */
/*             = 'U':  The given matrix Q is updated by the unitary */
/*                     row permutations used in the reduction. */

/*     JOBZ    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the unitary column transformations, as follows: */
/*             = 'N':  Do not form Z; */
/*             = 'I':  Z is initialized to the unit matrix and the */
/*                     unitary transformation matrix Z is returned; */
/*             = 'U':  The given matrix Z is updated by the unitary */
/*                     transformations used in the reduction. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows in the matrices A, E and the order of */
/*             the matrix Q.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns in the matrices A, E and the order */
/*             of the matrix Z.  N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the A matrix of the pencil sE-A. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the unitary transformed matrix Q' * A * Z. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the E matrix of the pencil sE-A, to be reduced to */
/*             column echelon form. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the unitary transformed matrix Q' * E * Z, which is in */
/*             column echelon form. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,M). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,*) */
/*             On entry, if JOBQ = 'U', then the leading M-by-M part of */
/*             this array must contain a given matrix Q (e.g. from a */
/*             previous call to another SLICOT routine), and on exit, the */
/*             leading M-by-M part of this array contains the product of */
/*             the input matrix Q and the row permutation matrix used to */
/*             transform the rows of matrix E. */
/*             On exit, if JOBQ = 'I', then the leading M-by-M part of */
/*             this array contains the matrix of accumulated unitary */
/*             row transformations performed. */
/*             If JOBQ = 'N', the array Q is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDQ = 1 and */
/*             declare this array to be Q(1,1) in the calling program). */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. If JOBQ = 'U' or */
/*             JOBQ = 'I', LDQ >= MAX(1,M); if JOBQ = 'N', LDQ >= 1. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,*) */
/*             On entry, if JOBZ = 'U', then the leading N-by-N part of */
/*             this array must contain a given matrix Z (e.g. from a */
/*             previous call to another SLICOT routine), and on exit, the */
/*             leading N-by-N part of this array contains the product of */
/*             the input matrix Z and the column transformation matrix */
/*             used to transform the columns of matrix E. */
/*             On exit, if JOBZ = 'I', then the leading N-by-N part of */
/*             this array contains the matrix of accumulated unitary */
/*             column transformations performed. */
/*             If JOBZ = 'N', the array Z is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDZ = 1 and */
/*             declare this array to be Z(1,1) in the calling program). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. If JOBZ = 'U' or */
/*             JOBZ = 'I', LDZ >= MAX(1,N); if JOBZ = 'N', LDZ >= 1. */

/*     RANKE   (output) INTEGER */
/*             The computed rank of the unitary transformed matrix E. */

/*     ISTAIR  (output) INTEGER array, dimension (M) */
/*             This array contains information on the column echelon form */
/*             of the unitary transformed matrix E. Specifically, */
/*             ISTAIR(i) = +j if the first non-zero element E(i,j) */
/*             is a corner point and -j otherwise, for i = 1,2,...,M. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             A tolerance below which matrix elements are considered */
/*             to be zero. If the user sets TOL to be less than (or */
/*             equal to) zero then the tolerance is taken as */
/*             EPS * MAX(ABS(E(I,J))), where EPS is the machine */
/*             precision (see LAPACK Library routine DLAMCH), */
/*             I = 1,2,...,M and J = 1,2,...,N. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension MAX(M,N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given an M-by-N matrix pencil sE-A with E not necessarily regular, */
/*     the routine computes a unitary transformed pencil Q'(sE-A)Z such */
/*     that the matrix Q' * E * Z is in column echelon form (trapezoidal */
/*     form).  Further details can be found in [1]. */

/*     [An M-by-N matrix E with rank(E) = r is said to be in column */
/*     echelon form if the following conditions are satisfied: */
/*     (a) the first (N - r) columns contain only zero elements; and */
/*     (b) if E(i(k),k) is the last nonzero element in column k for */
/*         k = N-r+1,...,N, i.e. E(i(k),k) <> 0 and E(j,k) = 0 for */
/*         j > i(k), then 1 <= i(N-r+1) < i(N-r+2) < ... < i(N) <= M.] */

/*     REFERENCES */

/*     [1] Beelen, Th. and Van Dooren, P. */
/*         An improved algorithm for the computation of Kronecker's */
/*         canonical form of a singular pencil. */
/*         Linear Algebra and Applications, 105, pp. 9-65, 1988. */

/*     NUMERICAL ASPECTS */

/*     It is shown in [1] that the algorithm is numerically backward */
/*     stable. The operations count is proportional to (MAX(M,N))**3. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998. */
/*     Based on Release 3.0 routine MB04SD modified by A. Varga, */
/*     German Aerospace Research Establishment, Oberpfaffenhofen, */
/*     Germany, Dec. 1997, to transform also the matrix A. */

/*     REVISIONS */

/*     A. Varga, DLR Oberpfaffenhofen, June 2005. */

/*     KEYWORDS */

/*     Echelon form, orthogonal transformation, staircase form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 214 "MB04UD.f"
    /* Parameter adjustments */
#line 214 "MB04UD.f"
    a_dim1 = *lda;
#line 214 "MB04UD.f"
    a_offset = 1 + a_dim1;
#line 214 "MB04UD.f"
    a -= a_offset;
#line 214 "MB04UD.f"
    e_dim1 = *lde;
#line 214 "MB04UD.f"
    e_offset = 1 + e_dim1;
#line 214 "MB04UD.f"
    e -= e_offset;
#line 214 "MB04UD.f"
    q_dim1 = *ldq;
#line 214 "MB04UD.f"
    q_offset = 1 + q_dim1;
#line 214 "MB04UD.f"
    q -= q_offset;
#line 214 "MB04UD.f"
    z_dim1 = *ldz;
#line 214 "MB04UD.f"
    z_offset = 1 + z_dim1;
#line 214 "MB04UD.f"
    z__ -= z_offset;
#line 214 "MB04UD.f"
    --istair;
#line 214 "MB04UD.f"
    --dwork;
#line 214 "MB04UD.f"

#line 214 "MB04UD.f"
    /* Function Body */
#line 214 "MB04UD.f"
    *info = 0;
#line 215 "MB04UD.f"
    ljobqi = lsame_(jobq, "I", (ftnlen)1, (ftnlen)1);
#line 216 "MB04UD.f"
    updatq = ljobqi || lsame_(jobq, "U", (ftnlen)1, (ftnlen)1);
#line 217 "MB04UD.f"
    ljobzi = lsame_(jobz, "I", (ftnlen)1, (ftnlen)1);
#line 218 "MB04UD.f"
    updatz = ljobzi || lsame_(jobz, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 222 "MB04UD.f"
    if (! updatq && ! lsame_(jobq, "N", (ftnlen)1, (ftnlen)1)) {
#line 223 "MB04UD.f"
	*info = -1;
#line 224 "MB04UD.f"
    } else if (! updatz && ! lsame_(jobz, "N", (ftnlen)1, (ftnlen)1)) {
#line 225 "MB04UD.f"
	*info = -2;
#line 226 "MB04UD.f"
    } else if (*m < 0) {
#line 227 "MB04UD.f"
	*info = -3;
#line 228 "MB04UD.f"
    } else if (*n < 0) {
#line 229 "MB04UD.f"
	*info = -4;
#line 230 "MB04UD.f"
    } else if (*lda < max(1,*m)) {
#line 231 "MB04UD.f"
	*info = -6;
#line 232 "MB04UD.f"
    } else if (*lde < max(1,*m)) {
#line 233 "MB04UD.f"
	*info = -8;
#line 234 "MB04UD.f"
    } else if (! updatq && *ldq < 1 || updatq && *ldq < max(1,*m)) {
#line 236 "MB04UD.f"
	*info = -10;
#line 237 "MB04UD.f"
    } else if (! updatz && *ldz < 1 || updatz && *ldz < max(1,*n)) {
#line 239 "MB04UD.f"
	*info = -12;
#line 240 "MB04UD.f"
    }

#line 242 "MB04UD.f"
    if (*info != 0) {

/*        Error return. */

#line 246 "MB04UD.f"
	i__1 = -(*info);
#line 246 "MB04UD.f"
	xerbla_("MB04UD", &i__1, (ftnlen)6);
#line 247 "MB04UD.f"
	return 0;
#line 248 "MB04UD.f"
    }

/*     Initialize Q and Z to the identity matrices, if needed. */

#line 252 "MB04UD.f"
    if (ljobqi) {
#line 252 "MB04UD.f"
	dlaset_("Full", m, m, &c_b10, &c_b11, &q[q_offset], ldq, (ftnlen)4);
#line 252 "MB04UD.f"
    }
#line 254 "MB04UD.f"
    if (ljobzi) {
#line 254 "MB04UD.f"
	dlaset_("Full", n, n, &c_b10, &c_b11, &z__[z_offset], ldz, (ftnlen)4);
#line 254 "MB04UD.f"
    }

/*     Quick return if possible. */

#line 259 "MB04UD.f"
    *ranke = min(*m,*n);

#line 261 "MB04UD.f"
    if (*ranke == 0) {
#line 261 "MB04UD.f"
	return 0;
#line 261 "MB04UD.f"
    }

#line 264 "MB04UD.f"
    toler = *tol;
#line 265 "MB04UD.f"
    if (toler <= 0.) {
#line 265 "MB04UD.f"
	toler = dlamch_("Epsilon", (ftnlen)7) * dlange_("M", m, n, &e[
		e_offset], lde, &dwork[1], (ftnlen)1);
#line 265 "MB04UD.f"
    }

#line 268 "MB04UD.f"
    k = *n;
#line 269 "MB04UD.f"
    lzero = FALSE_;

/*     WHILE ( ( K > 0 ) AND ( NOT a zero submatrix encountered ) ) DO */
#line 272 "MB04UD.f"
L20:
#line 272 "MB04UD.f"
    if (k > 0 && ! lzero) {

/*         Intermediate form of E */

/*                     <--k--><--n-k-> */
/*                l=1 |x....x|       | */
/*                    |      |       | */
/*                    |  Ek  |   X   | */
/*                    |      |       | */
/*            l=m-n+k |x....x|       | */
/*                    ---------------- */
/*                    |      |x ... x|  } */
/*                    |  O   |  x x x|  } */
/*                    |      |    x x|  } n-k */
/*                    |      |      x|  } */

/*        where submatrix Ek = E[1:m-n+k;1:k]. */

/*        Determine row LK in submatrix Ek with largest max-norm */
/*        (starting with row m-n+k). */

#line 293 "MB04UD.f"
	mnk = *m - *n + k;
#line 294 "MB04UD.f"
	emxnrm = 0.;
#line 295 "MB04UD.f"
	lk = mnk;

#line 297 "MB04UD.f"
	for (l = mnk; l >= 1; --l) {
#line 298 "MB04UD.f"
	    emx = (d__1 = e[l + idamax_(&k, &e[l + e_dim1], lde) * e_dim1], 
		    abs(d__1));
#line 299 "MB04UD.f"
	    if (emx > emxnrm) {
#line 300 "MB04UD.f"
		emxnrm = emx;
#line 301 "MB04UD.f"
		lk = l;
#line 302 "MB04UD.f"
	    }
#line 303 "MB04UD.f"
/* L40: */
#line 303 "MB04UD.f"
	}

#line 305 "MB04UD.f"
	if (emxnrm <= toler) {

/*           Set submatrix Ek to zero. */

#line 309 "MB04UD.f"
	    dlaset_("Full", &mnk, &k, &c_b10, &c_b10, &e[e_offset], lde, (
		    ftnlen)4);
#line 310 "MB04UD.f"
	    lzero = TRUE_;
#line 311 "MB04UD.f"
	    *ranke = *n - k;
#line 312 "MB04UD.f"
	} else {

/*           Submatrix Ek is not considered to be identically zero. */
/*           Check whether rows have to be interchanged. */

#line 317 "MB04UD.f"
	    if (lk != mnk) {

/*              Interchange rows lk and m-n+k in whole A- and E-matrix */
/*              and update the row transformation matrix Q, if needed. */
/*              (For Q, the number of elements involved is m.) */

#line 323 "MB04UD.f"
		dswap_(n, &e[lk + e_dim1], lde, &e[mnk + e_dim1], lde);
#line 324 "MB04UD.f"
		dswap_(n, &a[lk + a_dim1], lda, &a[mnk + a_dim1], lda);
#line 325 "MB04UD.f"
		if (updatq) {
#line 325 "MB04UD.f"
		    dswap_(m, &q[lk * q_dim1 + 1], &c__1, &q[mnk * q_dim1 + 1]
			    , &c__1);
#line 325 "MB04UD.f"
		}
#line 326 "MB04UD.f"
	    }

#line 328 "MB04UD.f"
	    km1 = k - 1;

/*           Determine a Householder transformation to annihilate */
/*           E(m-n+k,1:k-1) using E(m-n+k,k) as pivot. */
/*           Apply the transformation to the columns of A and Ek */
/*           (number of elements involved is m for A and m-n+k for Ek). */
/*           Update the column transformation matrix Z, if needed */
/*           (number of elements involved is n). */

#line 337 "MB04UD.f"
	    dlarfg_(&k, &e[mnk + k * e_dim1], &e[mnk + e_dim1], lde, &tau);
#line 338 "MB04UD.f"
	    emx = e[mnk + k * e_dim1];
#line 339 "MB04UD.f"
	    e[mnk + k * e_dim1] = 1.;
#line 340 "MB04UD.f"
	    i__1 = mnk - 1;
#line 340 "MB04UD.f"
	    dlarf_("Right", &i__1, &k, &e[mnk + e_dim1], lde, &tau, &e[
		    e_offset], lde, &dwork[1], (ftnlen)5);
#line 342 "MB04UD.f"
	    dlarf_("Right", m, &k, &e[mnk + e_dim1], lde, &tau, &a[a_offset], 
		    lda, &dwork[1], (ftnlen)5);
#line 344 "MB04UD.f"
	    if (updatz) {
#line 344 "MB04UD.f"
		dlarf_("Right", n, &k, &e[mnk + e_dim1], lde, &tau, &z__[
			z_offset], ldz, &dwork[1], (ftnlen)5);
#line 344 "MB04UD.f"
	    }
#line 346 "MB04UD.f"
	    e[mnk + k * e_dim1] = emx;
#line 347 "MB04UD.f"
	    dlaset_("Full", &c__1, &km1, &c_b10, &c_b10, &e[mnk + e_dim1], 
		    lde, (ftnlen)4);

#line 349 "MB04UD.f"
	    k = km1;
#line 350 "MB04UD.f"
	}
#line 351 "MB04UD.f"
	goto L20;
#line 352 "MB04UD.f"
    }
/*     END WHILE 20 */

/*     Initialise administration staircase form, i.e. */
/*     ISTAIR(i) =  j  if E(i,j) is a nonzero corner point */
/*               = -j  if E(i,j) is on the boundary but is no corner */
/*                     point. */
/*     Thus, */
/*     ISTAIR(m-k) =   n-k           for k=0,...,rank(E)-1 */
/*                 = -(n-rank(E)+1)  for k=rank(E),...,m-1. */

#line 363 "MB04UD.f"
    i__1 = *ranke - 1;
#line 363 "MB04UD.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 364 "MB04UD.f"
	istair[*m - i__] = *n - i__;
#line 365 "MB04UD.f"
/* L60: */
#line 365 "MB04UD.f"
    }

#line 367 "MB04UD.f"
    nr1 = -(*n - *ranke + 1);

#line 369 "MB04UD.f"
    i__1 = *m - *ranke;
#line 369 "MB04UD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 370 "MB04UD.f"
	istair[i__] = nr1;
#line 371 "MB04UD.f"
/* L80: */
#line 371 "MB04UD.f"
    }

#line 373 "MB04UD.f"
    return 0;
/* *** Last line of MB04UD *** */
} /* mb04ud_ */

