#line 1 "AB07ND.f"
/* AB07ND.f -- translated by f2c (version 20100827).
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

#line 1 "AB07ND.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b15 = -1.;
static doublereal c_b16 = 0.;
static doublereal c_b31 = 1.;

/* Subroutine */ int ab07nd_(integer *n, integer *m, doublereal *a, integer *
	lda, doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *rcond, integer *iwork, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, bl, ierr;
    static logical blas3;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical block;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static integer chunk;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal dnorm;
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
	    integer *, doublereal *, integer *, integer *);
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

/*     To compute the inverse (Ai,Bi,Ci,Di) of a given system (A,B,C,D). */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs and outputs.  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state matrix A of the original system. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the state matrix Ai of the inverse system. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B of the original system. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the input matrix Bi of the inverse system. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the output matrix C of the original system. */
/*             On exit, the leading M-by-N part of this array contains */
/*             the output matrix Ci of the inverse system. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= MAX(1,M). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the feedthrough matrix D of the original system. */
/*             On exit, the leading M-by-M part of this array contains */
/*             the feedthrough matrix Di of the inverse system. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= MAX(1,M). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The estimated reciprocal condition number of the */
/*             feedthrough matrix D of the original system. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0 or M+1, DWORK(1) returns the optimal */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,4*M). */
/*             For good performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = i:  the matrix D is exactly singular; the (i,i) diagonal */
/*                   element is zero, i <= M; RCOND was set to zero; */
/*             = M+1:  the matrix D is numerically singular, i.e., RCOND */
/*                   is less than the relative machine precision, EPS */
/*                   (see LAPACK Library routine DLAMCH). The */
/*                   calculations have been completed, but the results */
/*                   could be very inaccurate. */

/*     METHOD */

/*     The matrices of the inverse system are computed with the formulas: */
/*                   -1              -1         -1           -1 */
/*       Ai = A - B*D  *C,  Bi = -B*D  ,  Ci = D  *C,  Di = D  . */

/*     NUMERICAL ASPECTS */

/*     The accuracy depends mainly on the condition number of the matrix */
/*     D to be inverted. The estimated reciprocal condition number is */
/*     returned in RCOND. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, March 2000. */
/*     D. Sima, University of Bucharest, April 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */
/*     Based on the routine SYSINV, A. Varga, 1992. */

/*     REVISIONS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, July 2000. */

/*     KEYWORDS */

/*     Inverse system, state-space model, state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 156 "AB07ND.f"
    /* Parameter adjustments */
#line 156 "AB07ND.f"
    a_dim1 = *lda;
#line 156 "AB07ND.f"
    a_offset = 1 + a_dim1;
#line 156 "AB07ND.f"
    a -= a_offset;
#line 156 "AB07ND.f"
    b_dim1 = *ldb;
#line 156 "AB07ND.f"
    b_offset = 1 + b_dim1;
#line 156 "AB07ND.f"
    b -= b_offset;
#line 156 "AB07ND.f"
    c_dim1 = *ldc;
#line 156 "AB07ND.f"
    c_offset = 1 + c_dim1;
#line 156 "AB07ND.f"
    c__ -= c_offset;
#line 156 "AB07ND.f"
    d_dim1 = *ldd;
#line 156 "AB07ND.f"
    d_offset = 1 + d_dim1;
#line 156 "AB07ND.f"
    d__ -= d_offset;
#line 156 "AB07ND.f"
    --iwork;
#line 156 "AB07ND.f"
    --dwork;
#line 156 "AB07ND.f"

#line 156 "AB07ND.f"
    /* Function Body */
#line 156 "AB07ND.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 160 "AB07ND.f"
    if (*n < 0) {
#line 161 "AB07ND.f"
	*info = -1;
#line 162 "AB07ND.f"
    } else if (*m < 0) {
#line 163 "AB07ND.f"
	*info = -2;
#line 164 "AB07ND.f"
    } else if (*lda < max(1,*n)) {
#line 165 "AB07ND.f"
	*info = -4;
#line 166 "AB07ND.f"
    } else if (*ldb < max(1,*n)) {
#line 167 "AB07ND.f"
	*info = -6;
#line 168 "AB07ND.f"
    } else if (*ldc < max(1,*m)) {
#line 169 "AB07ND.f"
	*info = -8;
#line 170 "AB07ND.f"
    } else if (*ldd < max(1,*m)) {
#line 171 "AB07ND.f"
	*info = -10;
#line 172 "AB07ND.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 172 "AB07ND.f"
	i__1 = 1, i__2 = *m << 2;
#line 172 "AB07ND.f"
	if (*ldwork < max(i__1,i__2)) {
#line 173 "AB07ND.f"
	    *info = -14;
#line 174 "AB07ND.f"
	}
#line 174 "AB07ND.f"
    }

#line 176 "AB07ND.f"
    if (*info != 0) {

/*        Error return. */

#line 180 "AB07ND.f"
	i__1 = -(*info);
#line 180 "AB07ND.f"
	xerbla_("AB07ND", &i__1, (ftnlen)6);
#line 181 "AB07ND.f"
	return 0;
#line 182 "AB07ND.f"
    }

/*     Quick return if possible. */

#line 186 "AB07ND.f"
    if (*m == 0) {
#line 187 "AB07ND.f"
	*rcond = 1.;
#line 188 "AB07ND.f"
	dwork[1] = 1.;
#line 189 "AB07ND.f"
	return 0;
#line 190 "AB07ND.f"
    }

/*     Factorize D. */

#line 194 "AB07ND.f"
    dgetrf_(m, m, &d__[d_offset], ldd, &iwork[1], info);
#line 195 "AB07ND.f"
    if (*info != 0) {
#line 196 "AB07ND.f"
	*rcond = 0.;
#line 197 "AB07ND.f"
	return 0;
#line 198 "AB07ND.f"
    }

/*     Compute the reciprocal condition number of the matrix D. */
/*     Workspace: need   4*M. */
/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*      minimal amount of workspace needed at that point in the code, */
/*      as well as the preferred amount for good performance. */
/*      NB refers to the optimal block size for the immediately */
/*      following subroutine, as returned by ILAENV.) */

#line 208 "AB07ND.f"
    dnorm = dlange_("1-norm", m, m, &d__[d_offset], ldd, &dwork[1], (ftnlen)6)
	    ;
#line 209 "AB07ND.f"
    dgecon_("1-norm", m, &d__[d_offset], ldd, &dnorm, rcond, &dwork[1], &
	    iwork[*m + 1], &ierr, (ftnlen)6);
#line 211 "AB07ND.f"
    if (*rcond < dlamch_("Epsilon", (ftnlen)7)) {
#line 211 "AB07ND.f"
	*info = *m + 1;
#line 211 "AB07ND.f"
    }
/*                   -1 */
/*     Compute Di = D  . */
/*     Workspace: need   M; */
/*                prefer M*NB. */

/* Computing MAX */
#line 218 "AB07ND.f"
    i__1 = *m << 2, i__2 = *m * ilaenv_(&c__1, "DGETRI", " ", m, &c_n1, &c_n1,
	     &c_n1, (ftnlen)6, (ftnlen)1);
#line 218 "AB07ND.f"
    maxwrk = max(i__1,i__2);
#line 219 "AB07ND.f"
    dgetri_(m, &d__[d_offset], ldd, &iwork[1], &dwork[1], ldwork, &ierr);
#line 220 "AB07ND.f"
    if (*n > 0) {
#line 221 "AB07ND.f"
	chunk = *ldwork / *m;
#line 222 "AB07ND.f"
	blas3 = chunk >= *n && *m > 1;
#line 223 "AB07ND.f"
	block = min(chunk,*m) > 1;
/*                          -1 */
/*        Compute  Bi = -B*D  . */

#line 227 "AB07ND.f"
	if (blas3) {

/*           Enough workspace for a fast BLAS 3 algorithm. */

#line 231 "AB07ND.f"
	    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[1], n, (ftnlen)4);
#line 232 "AB07ND.f"
	    dgemm_("NoTranspose", "NoTranspose", n, m, m, &c_b15, &dwork[1], 
		    n, &d__[d_offset], ldd, &c_b16, &b[b_offset], ldb, (
		    ftnlen)11, (ftnlen)11);

#line 235 "AB07ND.f"
	} else if (block) {

/*           Use as many rows of B as possible. */

#line 239 "AB07ND.f"
	    i__1 = *n;
#line 239 "AB07ND.f"
	    i__2 = chunk;
#line 239 "AB07ND.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
/* Computing MIN */
#line 240 "AB07ND.f"
		i__3 = *n - i__ + 1;
#line 240 "AB07ND.f"
		bl = min(i__3,chunk);
#line 241 "AB07ND.f"
		dlacpy_("Full", &bl, m, &b[i__ + b_dim1], ldb, &dwork[1], &bl,
			 (ftnlen)4);
#line 242 "AB07ND.f"
		dgemm_("NoTranspose", "NoTranspose", &bl, m, m, &c_b15, &
			dwork[1], &bl, &d__[d_offset], ldd, &c_b16, &b[i__ + 
			b_dim1], ldb, (ftnlen)11, (ftnlen)11);
#line 244 "AB07ND.f"
/* L10: */
#line 244 "AB07ND.f"
	    }

#line 246 "AB07ND.f"
	} else {

/*           Use a BLAS 2 algorithm. */

#line 250 "AB07ND.f"
	    i__2 = *n;
#line 250 "AB07ND.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 251 "AB07ND.f"
		dcopy_(m, &b[i__ + b_dim1], ldb, &dwork[1], &c__1);
#line 252 "AB07ND.f"
		dgemv_("Transpose", m, m, &c_b15, &d__[d_offset], ldd, &dwork[
			1], &c__1, &c_b16, &b[i__ + b_dim1], ldb, (ftnlen)9);
#line 254 "AB07ND.f"
/* L20: */
#line 254 "AB07ND.f"
	    }

#line 256 "AB07ND.f"
	}

/*        Compute  Ai = A + Bi*C. */

#line 260 "AB07ND.f"
	dgemm_("NoTranspose", "NoTranspose", n, n, m, &c_b31, &b[b_offset], 
		ldb, &c__[c_offset], ldc, &c_b31, &a[a_offset], lda, (ftnlen)
		11, (ftnlen)11);
/*                        -1 */
/*        Compute  C <-- D  *C. */

#line 265 "AB07ND.f"
	if (blas3) {

/*           Enough workspace for a fast BLAS 3 algorithm. */

#line 269 "AB07ND.f"
	    dlacpy_("Full", m, n, &c__[c_offset], ldc, &dwork[1], m, (ftnlen)
		    4);
#line 270 "AB07ND.f"
	    dgemm_("NoTranspose", "NoTranspose", m, n, m, &c_b31, &d__[
		    d_offset], ldd, &dwork[1], m, &c_b16, &c__[c_offset], ldc,
		     (ftnlen)11, (ftnlen)11);

#line 273 "AB07ND.f"
	} else if (block) {

/*           Use as many columns of C as possible. */

#line 277 "AB07ND.f"
	    i__2 = *n;
#line 277 "AB07ND.f"
	    i__1 = chunk;
#line 277 "AB07ND.f"
	    for (j = 1; i__1 < 0 ? j >= i__2 : j <= i__2; j += i__1) {
/* Computing MIN */
#line 278 "AB07ND.f"
		i__3 = *n - j + 1;
#line 278 "AB07ND.f"
		bl = min(i__3,chunk);
#line 279 "AB07ND.f"
		dlacpy_("Full", m, &bl, &c__[j * c_dim1 + 1], ldc, &dwork[1], 
			m, (ftnlen)4);
#line 280 "AB07ND.f"
		dgemm_("NoTranspose", "NoTranspose", m, &bl, m, &c_b31, &d__[
			d_offset], ldd, &dwork[1], m, &c_b16, &c__[j * c_dim1 
			+ 1], ldc, (ftnlen)11, (ftnlen)11);
#line 282 "AB07ND.f"
/* L30: */
#line 282 "AB07ND.f"
	    }

#line 284 "AB07ND.f"
	} else {

/*           Use a BLAS 2 algorithm. */

#line 288 "AB07ND.f"
	    i__1 = *n;
#line 288 "AB07ND.f"
	    for (j = 1; j <= i__1; ++j) {
#line 289 "AB07ND.f"
		dcopy_(m, &c__[j * c_dim1 + 1], &c__1, &dwork[1], &c__1);
#line 290 "AB07ND.f"
		dgemv_("NoTranspose", m, m, &c_b31, &d__[d_offset], ldd, &
			dwork[1], &c__1, &c_b16, &c__[j * c_dim1 + 1], &c__1, 
			(ftnlen)11);
#line 292 "AB07ND.f"
/* L40: */
#line 292 "AB07ND.f"
	    }

#line 294 "AB07ND.f"
	}
#line 295 "AB07ND.f"
    }

/*     Return optimal workspace in DWORK(1). */

/* Computing MAX */
#line 299 "AB07ND.f"
    i__1 = maxwrk, i__2 = *n * *m;
#line 299 "AB07ND.f"
    dwork[1] = (doublereal) max(i__1,i__2);
#line 300 "AB07ND.f"
    return 0;

/* *** Last line of AB07ND *** */
} /* ab07nd_ */

