#line 1 "MB04LD.f"
/* MB04LD.f -- translated by f2c (version 20100827).
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

#line 1 "MB04LD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b7 = 1.;
static doublereal c_b12 = 0.;

/* Subroutine */ int mb04ld_(char *uplo, integer *n, integer *m, integer *p, 
	doublereal *l, integer *ldl, doublereal *a, integer *lda, doublereal *
	b, integer *ldb, doublereal *c__, integer *ldc, doublereal *tau, 
	doublereal *dwork, ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, l_dim1, 
	    l_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, im;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *), dscal_(integer *, doublereal *, doublereal *, integer 
	    *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), daxpy_(integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *)
	    ;
    static logical luplo;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);


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

/*     To calculate an LQ factorization of the first block row and apply */
/*     the orthogonal transformations (from the right) also to the second */
/*     block row of a structured matrix, as follows */
/*                        _ */
/*        [ L   A ]     [ L   0 ] */
/*        [       ]*Q = [       ] */
/*        [ 0   B ]     [ C   D ] */
/*                 _ */
/*     where L and L are lower triangular. The matrix A can be full or */
/*     lower trapezoidal/triangular. The problem structure is exploited. */
/*     This computation is useful, for instance, in combined measurement */
/*     and time update of one iteration of the Kalman filter (square */
/*     root covariance filter). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPLO    CHARACTER*1 */
/*             Indicates if the matrix A is or not triangular as follows: */
/*             = 'L':  Matrix A is lower trapezoidal/triangular; */
/*             = 'F':  Matrix A is full. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER                 _ */
/*             The order of the matrices L and L.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrices A, B and D.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of rows of the matrices B, C and D.  P >= 0. */

/*     L       (input/output) DOUBLE PRECISION array, dimension (LDL,N) */
/*             On entry, the leading N-by-N lower triangular part of this */
/*             array must contain the lower triangular matrix L. */
/*             On exit, the leading N-by-N lower triangular part of this */
/*                                                        _ */
/*             array contains the lower triangular matrix L. */
/*             The strict upper triangular part of this array is not */
/*             referenced. */

/*     LDL     INTEGER */
/*             The leading dimension of array L.  LDL >= MAX(1,N). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,M) */
/*             On entry, if UPLO = 'F', the leading N-by-M part of this */
/*             array must contain the matrix A. If UPLO = 'L', the */
/*             leading N-by-MIN(N,M) part of this array must contain the */
/*             lower trapezoidal (lower triangular if N <= M) matrix A, */
/*             and the elements above the diagonal are not referenced. */
/*             On exit, the leading N-by-M part (lower trapezoidal or */
/*             triangular, if UPLO = 'L') of this array contains the */
/*             trailing components (the vectors v, see Method) of the */
/*             elementary reflectors used in the factorization. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the matrix B. */
/*             On exit, the leading P-by-M part of this array contains */
/*             the computed matrix D. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,P). */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array contains the */
/*             computed matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     TAU     (output) DOUBLE PRECISION array, dimension (N) */
/*             The scalar factors of the elementary reflectors used. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     METHOD */

/*     The routine uses N Householder transformations exploiting the zero */
/*     pattern of the block matrix.  A Householder matrix has the form */

/*                                     ( 1 ), */
/*        H  = I - tau *u *u',    u  = ( v ) */
/*         i          i  i  i      i   (  i) */

/*     where v  is an M-vector, if UPLO = 'F', or an min(i,M)-vector, if */
/*            i */
/*     UPLO = 'L'.  The components of v  are stored in the i-th row of A, */
/*                                     i */
/*     and tau  is stored in TAU(i). */
/*            i */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1997. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary reflector, LQ factorization, orthogonal transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 162 "MB04LD.f"
    /* Parameter adjustments */
#line 162 "MB04LD.f"
    l_dim1 = *ldl;
#line 162 "MB04LD.f"
    l_offset = 1 + l_dim1;
#line 162 "MB04LD.f"
    l -= l_offset;
#line 162 "MB04LD.f"
    a_dim1 = *lda;
#line 162 "MB04LD.f"
    a_offset = 1 + a_dim1;
#line 162 "MB04LD.f"
    a -= a_offset;
#line 162 "MB04LD.f"
    b_dim1 = *ldb;
#line 162 "MB04LD.f"
    b_offset = 1 + b_dim1;
#line 162 "MB04LD.f"
    b -= b_offset;
#line 162 "MB04LD.f"
    c_dim1 = *ldc;
#line 162 "MB04LD.f"
    c_offset = 1 + c_dim1;
#line 162 "MB04LD.f"
    c__ -= c_offset;
#line 162 "MB04LD.f"
    --tau;
#line 162 "MB04LD.f"
    --dwork;
#line 162 "MB04LD.f"

#line 162 "MB04LD.f"
    /* Function Body */
#line 162 "MB04LD.f"
    if (min(*m,*n) == 0) {
#line 162 "MB04LD.f"
	return 0;
#line 162 "MB04LD.f"
    }

#line 165 "MB04LD.f"
    luplo = lsame_(uplo, "L", (ftnlen)1, (ftnlen)1);
#line 166 "MB04LD.f"
    im = *m;

#line 168 "MB04LD.f"
    i__1 = *n;
#line 168 "MB04LD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {

/*        Annihilate the I-th row of A and apply the transformations to */
/*        the entire block matrix, exploiting its structure. */

#line 173 "MB04LD.f"
	if (luplo) {
#line 173 "MB04LD.f"
	    im = min(i__,*m);
#line 173 "MB04LD.f"
	}
#line 174 "MB04LD.f"
	i__2 = im + 1;
#line 174 "MB04LD.f"
	dlarfg_(&i__2, &l[i__ + i__ * l_dim1], &a[i__ + a_dim1], lda, &tau[
		i__]);
#line 175 "MB04LD.f"
	if (tau[i__] != 0.) {

/*           [    w   ]    [ L(I+1:N,I) A(I+1:N,1:IM) ]   [ 1 ] */
/*           [        ] := [                          ] * [   ] */
/*           [ C(:,I) ]    [      0        B(:,1:IM)  ]   [ v ] */

#line 181 "MB04LD.f"
	    if (i__ < *n) {
#line 182 "MB04LD.f"
		i__2 = *n - i__;
#line 182 "MB04LD.f"
		dcopy_(&i__2, &l[i__ + 1 + i__ * l_dim1], &c__1, &dwork[1], &
			c__1);
#line 183 "MB04LD.f"
		i__2 = *n - i__;
#line 183 "MB04LD.f"
		dgemv_("No transpose", &i__2, &im, &c_b7, &a[i__ + 1 + a_dim1]
			, lda, &a[i__ + a_dim1], lda, &c_b7, &dwork[1], &c__1,
			 (ftnlen)12);
#line 185 "MB04LD.f"
	    }
#line 186 "MB04LD.f"
	    dgemv_("No transpose", p, &im, &c_b7, &b[b_offset], ldb, &a[i__ + 
		    a_dim1], lda, &c_b12, &c__[i__ * c_dim1 + 1], &c__1, (
		    ftnlen)12);

/*           [ L(I+1:N,I) A(I+1:N,1:IM) ]    [ L(I+1:N,I) A(I+1:N,1:IM) ] */
/*           [                          ] := [                          ] */
/*           [   C(:,I)      D(:,1:IM)  ]    [      0        B(:,1:IM)  ] */

/*                                                 [    w   ] */
/*                                         - tau * [        ] * [ 1 , v'] */
/*                                                 [ C(:,I) ] */

#line 197 "MB04LD.f"
	    if (i__ < *n) {
#line 198 "MB04LD.f"
		i__2 = *n - i__;
#line 198 "MB04LD.f"
		d__1 = -tau[i__];
#line 198 "MB04LD.f"
		daxpy_(&i__2, &d__1, &dwork[1], &c__1, &l[i__ + 1 + i__ * 
			l_dim1], &c__1);
#line 199 "MB04LD.f"
		i__2 = *n - i__;
#line 199 "MB04LD.f"
		d__1 = -tau[i__];
#line 199 "MB04LD.f"
		dger_(&i__2, &im, &d__1, &dwork[1], &c__1, &a[i__ + a_dim1], 
			lda, &a[i__ + 1 + a_dim1], lda);
#line 201 "MB04LD.f"
	    }
#line 202 "MB04LD.f"
	    d__1 = -tau[i__];
#line 202 "MB04LD.f"
	    dscal_(p, &d__1, &c__[i__ * c_dim1 + 1], &c__1);
#line 203 "MB04LD.f"
	    dger_(p, &im, &c_b7, &c__[i__ * c_dim1 + 1], &c__1, &a[i__ + 
		    a_dim1], lda, &b[b_offset], ldb);
#line 204 "MB04LD.f"
	}
#line 205 "MB04LD.f"
/* L10: */
#line 205 "MB04LD.f"
    }

#line 207 "MB04LD.f"
    return 0;
/* *** Last line of MB04LD *** */
} /* mb04ld_ */

