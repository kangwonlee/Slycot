#line 1 "SB08CD.f"
/* SB08CD.f -- translated by f2c (version 20100827).
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

#line 1 "SB08CD.f"
/* Table of constant values */

static doublereal c_b6 = 0.;
static doublereal c_b7 = 1.;
static integer c__1 = 1;

/* Subroutine */ int sb08cd_(char *dico, integer *n, integer *m, integer *p, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, integer *nq, 
	integer *nr, doublereal *br, integer *ldbr, doublereal *dr, integer *
	lddr, doublereal *tol, doublereal *dwork, integer *ldwork, integer *
	iwarn, integer *info, ftnlen dico_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, br_dim1, br_offset, c_dim1, 
	    c_offset, d_dim1, d_offset, dr_dim1, dr_offset, i__1, i__2, i__3, 
	    i__4, i__5, i__6;

    /* Local variables */
    static integer i__, kw, kbr;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ma02bd_(char *, integer *, integer *, doublereal *, integer *, 
	    ftnlen), ab07md_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen), sb08dd_(
	    char *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

/*     To construct, for a given system G = (A,B,C,D), an output */
/*     injection matrix H, an orthogonal transformation matrix Z, and a */
/*     gain matrix V, such that the systems */

/*          Q = (Z'*(A+H*C)*Z, Z'*(B+H*D), V*C*Z, V*D) */
/*     and */
/*          R = (Z'*(A+H*C)*Z, Z'*H, V*C*Z, V) */

/*     provide a stable left coprime factorization of G in the form */
/*                   -1 */
/*              G = R  * Q, */

/*     where G, Q and R are the corresponding transfer-function matrices */
/*     and the denominator R is co-inner, that is, R(s)*R'(-s) = I in */
/*     the continuous-time case, or R(z)*R'(1/z) = I in the discrete-time */
/*     case. The Z matrix is not explicitly computed. */

/*     Note: G must have no observable poles on the imaginary axis */
/*     for a continuous-time system, or on the unit circle for a */
/*     discrete-time system. If the given state-space representation */
/*     is not detectable, the undetectable part of the original */
/*     system is automatically deflated and the order of the systems */
/*     Q and R is accordingly reduced. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the state vector, i.e. the order of the */
/*             matrix A, and also the number of rows of the matrices B */
/*             and BR, and the number of columns of the matrix C. */
/*             N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of input vector, i.e. the number of columns */
/*             of the matrices B and D.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of output vector, i.e. the number of rows */
/*             of the matrices C, D and DR, and the number of columns */
/*             of the matrices BR and DR.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. The matrix A must not */
/*             have observable eigenvalues on the imaginary axis, if */
/*             DICO = 'C', or on the unit circle, if DICO = 'D'. */
/*             On exit, the leading NQ-by-NQ part of this array contains */
/*             the leading NQ-by-NQ part of the matrix Z'*(A+H*C)*Z, the */
/*             state dynamics matrix of the numerator factor Q, in a */
/*             real Schur form. The leading NR-by-NR part of this matrix */
/*             represents the state dynamics matrix of a minimal */
/*             realization of the denominator factor R. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,MAX(M,P)) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input/state matrix. */
/*             On exit, the leading NQ-by-M part of this array contains */
/*             the leading NQ-by-M part of the matrix Z'*(B+H*D), the */
/*             input/state matrix of the numerator factor Q. */
/*             The remaining part of this array is needed as workspace. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C. */
/*             On exit, the leading P-by-NQ part of this array contains */
/*             the leading P-by-NQ part of the matrix V*C*Z, the */
/*             state/output matrix of the numerator factor Q. */
/*             The first NR columns of this array represent the */
/*             state/output matrix of a minimal realization of the */
/*             denominator factor R. */
/*             The remaining part of this array is needed as workspace. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,M,P), if N > 0. */
/*             LDC >= 1,          if N = 0. */

/*     D       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDD,MAX(M,P)) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the input/output matrix. */
/*             On exit, the leading P-by-M part of this array contains */
/*             the matrix V*D representing the input/output matrix */
/*             of the numerator factor Q. */
/*             The remaining part of this array is needed as workspace. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,M,P). */

/*     NQ      (output) INTEGER */
/*             The order of the resulting factors Q and R. */
/*             Generally, NQ = N - NS, where NS is the number of */
/*             unobservable eigenvalues outside the stability region. */

/*     NR      (output) INTEGER */
/*             The order of the minimal realization of the factor R. */
/*             Generally, NR is the number of observable eigenvalues */
/*             of A outside the stability region (the number of modified */
/*             eigenvalues). */

/*     BR      (output) DOUBLE PRECISION array, dimension (LDBR,P) */
/*             The leading NQ-by-P part of this array contains the */
/*             leading NQ-by-P part of the output injection matrix */
/*             Z'*H, which reflects the eigenvalues of A lying outside */
/*             the stable region to values which are symmetric with */
/*             respect to the imaginary axis (if DICO = 'C') or the unit */
/*             circle (if DICO = 'D'). The first NR rows of this matrix */
/*             form the input/state matrix of a minimal realization of */
/*             the denominator factor R. */

/*     LDBR    INTEGER */
/*             The leading dimension of array BR.  LDBR >= MAX(1,N). */

/*     DR      (output) DOUBLE PRECISION array, dimension (LDDR,P) */
/*             The leading P-by-P part of this array contains the lower */
/*             triangular matrix V representing the input/output matrix */
/*             of the denominator factor R. */

/*     LDDR    INTEGER */
/*             The leading dimension of array DR.  LDDR >= MAX(1,P). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The absolute tolerance level below which the elements of */
/*             C are considered zero (used for observability tests). */
/*             If the user sets TOL <= 0, then an implicitly computed, */
/*             default tolerance, defined by  TOLDEF = N*EPS*NORM(C), */
/*             is used instead, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH) and NORM(C) denotes */
/*             the infinity-norm of C. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of working array DWORK. */
/*             LDWORK >= MAX( 1, P*N + MAX( N*(N+5),P*(P+2),4*P,4*M ) ). */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = K:  K violations of the numerical stability condition */
/*                   NORM(H) <= 10*NORM(A)/NORM(C) occured during the */
/*                   assignment of eigenvalues. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the reduction of A to a real Schur form failed; */
/*             = 2:  a failure was detected during the ordering of the */
/*                   real Schur form of A, or in the iterative process */
/*                   for reordering the eigenvalues of Z'*(A + H*C)*Z */
/*                   along the diagonal; */
/*             = 3:  if DICO = 'C' and the matrix A has an observable */
/*                   eigenvalue on the imaginary axis, or DICO = 'D' and */
/*                   A has an observable eigenvalue on the unit circle. */

/*     METHOD */

/*     The subroutine uses the right coprime factorization algorithm with */
/*     inner denominator of [1] applied to G'. */

/*     REFERENCES */

/*     [1] Varga A. */
/*         A Schur method for computing coprime factorizations with */
/*         inner denominators and applications in model reduction. */
/*         Proc. ACC'93, San Francisco, CA, pp. 2130-2131, 1993. */

/*     NUMERICAL ASPECTS */
/*                                            3 */
/*     The algorithm requires no more than 14N  floating point */
/*     operations. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine LCFID. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     May  2003, A. Varga, DLR Oberpfaffenhofen. */
/*     Nov  2003, A. Varga, DLR Oberpfaffenhofen. */

/*     KEYWORDS */

/*     Coprime factorization, eigenvalue, eigenvalue assignment, */
/*     feedback control, pole placement, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 267 "SB08CD.f"
    /* Parameter adjustments */
#line 267 "SB08CD.f"
    a_dim1 = *lda;
#line 267 "SB08CD.f"
    a_offset = 1 + a_dim1;
#line 267 "SB08CD.f"
    a -= a_offset;
#line 267 "SB08CD.f"
    b_dim1 = *ldb;
#line 267 "SB08CD.f"
    b_offset = 1 + b_dim1;
#line 267 "SB08CD.f"
    b -= b_offset;
#line 267 "SB08CD.f"
    c_dim1 = *ldc;
#line 267 "SB08CD.f"
    c_offset = 1 + c_dim1;
#line 267 "SB08CD.f"
    c__ -= c_offset;
#line 267 "SB08CD.f"
    d_dim1 = *ldd;
#line 267 "SB08CD.f"
    d_offset = 1 + d_dim1;
#line 267 "SB08CD.f"
    d__ -= d_offset;
#line 267 "SB08CD.f"
    br_dim1 = *ldbr;
#line 267 "SB08CD.f"
    br_offset = 1 + br_dim1;
#line 267 "SB08CD.f"
    br -= br_offset;
#line 267 "SB08CD.f"
    dr_dim1 = *lddr;
#line 267 "SB08CD.f"
    dr_offset = 1 + dr_dim1;
#line 267 "SB08CD.f"
    dr -= dr_offset;
#line 267 "SB08CD.f"
    --dwork;
#line 267 "SB08CD.f"

#line 267 "SB08CD.f"
    /* Function Body */
#line 267 "SB08CD.f"
    *iwarn = 0;
#line 268 "SB08CD.f"
    *info = 0;

/*     Check the scalar input parameters. */

#line 272 "SB08CD.f"
    if (! lsame_(dico, "C", (ftnlen)1, (ftnlen)1) && ! lsame_(dico, "D", (
	    ftnlen)1, (ftnlen)1)) {
#line 274 "SB08CD.f"
	*info = -1;
#line 275 "SB08CD.f"
    } else if (*n < 0) {
#line 276 "SB08CD.f"
	*info = -2;
#line 277 "SB08CD.f"
    } else if (*m < 0) {
#line 278 "SB08CD.f"
	*info = -3;
#line 279 "SB08CD.f"
    } else if (*p < 0) {
#line 280 "SB08CD.f"
	*info = -4;
#line 281 "SB08CD.f"
    } else if (*lda < max(1,*n)) {
#line 282 "SB08CD.f"
	*info = -6;
#line 283 "SB08CD.f"
    } else if (*ldb < max(1,*n)) {
#line 284 "SB08CD.f"
	*info = -8;
#line 285 "SB08CD.f"
    } else if (*ldc < 1 || *n > 0 && *ldc < max(*m,*p)) {
#line 287 "SB08CD.f"
	*info = -10;
#line 288 "SB08CD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 288 "SB08CD.f"
	i__1 = max(1,*m);
#line 288 "SB08CD.f"
	if (*ldd < max(i__1,*p)) {
#line 289 "SB08CD.f"
	    *info = -12;
#line 290 "SB08CD.f"
	} else if (*ldbr < max(1,*n)) {
#line 291 "SB08CD.f"
	    *info = -16;
#line 292 "SB08CD.f"
	} else if (*lddr < max(1,*p)) {
#line 293 "SB08CD.f"
	    *info = -18;
#line 294 "SB08CD.f"
	} else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 294 "SB08CD.f"
	    i__3 = *n * (*n + 5), i__4 = *p * (*p + 2), i__3 = max(i__3,i__4),
		     i__4 = *p << 2, i__3 = max(i__3,i__4), i__4 = *m << 2;
#line 294 "SB08CD.f"
	    i__1 = 1, i__2 = *p * *n + max(i__3,i__4);
#line 294 "SB08CD.f"
	    if (*ldwork < max(i__1,i__2)) {
#line 296 "SB08CD.f"
		*info = -21;
#line 297 "SB08CD.f"
	    }
#line 297 "SB08CD.f"
	}
#line 297 "SB08CD.f"
    }
#line 298 "SB08CD.f"
    if (*info != 0) {

/*        Error return. */

#line 302 "SB08CD.f"
	i__1 = -(*info);
#line 302 "SB08CD.f"
	xerbla_("SB08CD", &i__1, (ftnlen)6);
#line 303 "SB08CD.f"
	return 0;
#line 304 "SB08CD.f"
    }

/*     Quick return if possible. */

#line 308 "SB08CD.f"
    if (min(*n,*p) == 0) {
#line 309 "SB08CD.f"
	*nq = 0;
#line 310 "SB08CD.f"
	*nr = 0;
#line 311 "SB08CD.f"
	dwork[1] = 1.;
#line 312 "SB08CD.f"
	dlaset_("Full", p, p, &c_b6, &c_b7, &dr[dr_offset], lddr, (ftnlen)4);
#line 313 "SB08CD.f"
	return 0;
#line 314 "SB08CD.f"
    }

/*     Compute the dual system G' = (A',C',B',D'). */

#line 318 "SB08CD.f"
    ab07md_("D", n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[c_offset]
	    , ldc, &d__[d_offset], ldd, info, (ftnlen)1);

/*     Compute the right coprime factorization with inner */
/*     denominator of G'. */

/*     Workspace needed:      P*N; */
/*     Additional workspace:  need  MAX( N*(N+5), P*(P+2), 4*P, 4*M ); */
/*                            prefer larger. */

#line 328 "SB08CD.f"
    kbr = 1;
#line 329 "SB08CD.f"
    kw = kbr + *p * *n;
#line 330 "SB08CD.f"
    i__1 = *ldwork - kw + 1;
#line 330 "SB08CD.f"
    sb08dd_(dico, n, p, m, &a[a_offset], lda, &b[b_offset], ldb, &c__[
	    c_offset], ldc, &d__[d_offset], ldd, nq, nr, &dwork[kbr], p, &dr[
	    dr_offset], lddr, tol, &dwork[kw], &i__1, iwarn, info, (ftnlen)1);
#line 333 "SB08CD.f"
    if (*info == 0) {

/*        Determine the elements of the left coprime factorization from */
/*        those of the computed right coprime factorization and make the */
/*        state-matrix upper real Schur. */

/* Computing MAX */
#line 339 "SB08CD.f"
	i__2 = 0, i__3 = *nq - 1;
#line 339 "SB08CD.f"
	i__1 = max(i__2,i__3);
/* Computing MAX */
#line 339 "SB08CD.f"
	i__5 = 0, i__6 = *nq - 1;
#line 339 "SB08CD.f"
	i__4 = max(i__5,i__6);
#line 339 "SB08CD.f"
	tb01xd_("D", nq, p, m, &i__1, &i__4, &a[a_offset], lda, &b[b_offset], 
		ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, info, (ftnlen)
		1);

#line 342 "SB08CD.f"
	ma02ad_("Full", p, nq, &dwork[kbr], p, &br[br_offset], ldbr, (ftnlen)
		4);
#line 343 "SB08CD.f"
	ma02bd_("Left", nq, p, &br[br_offset], ldbr, (ftnlen)4);

#line 345 "SB08CD.f"
	i__1 = *p;
#line 345 "SB08CD.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 346 "SB08CD.f"
	    i__2 = i__ - 1;
#line 346 "SB08CD.f"
	    dswap_(&i__2, &dr[i__ + dr_dim1], lddr, &dr[i__ * dr_dim1 + 1], &
		    c__1);
#line 347 "SB08CD.f"
/* L10: */
#line 347 "SB08CD.f"
	}

#line 349 "SB08CD.f"
    }

#line 351 "SB08CD.f"
    dwork[1] = dwork[kw] + (doublereal) (kw - 1);

#line 353 "SB08CD.f"
    return 0;
/* *** Last line of SB08CD *** */
} /* sb08cd_ */

