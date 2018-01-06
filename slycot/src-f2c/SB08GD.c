#line 1 "SB08GD.f"
/* SB08GD.f -- translated by f2c (version 20100827).
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

#line 1 "SB08GD.f"
/* Table of constant values */

static doublereal c_b7 = -1.;
static doublereal c_b8 = 1.;

/* Subroutine */ int sb08gd_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *br, integer *
	ldbr, doublereal *dr, integer *lddr, integer *iwork, doublereal *
	dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, br_dim1, br_offset, c_dim1, 
	    c_offset, d_dim1, d_offset, dr_dim1, dr_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal rcond;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen), dgetrs_(char *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen);
    static doublereal drnorm;


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

/*     To construct the state-space representation for the system */
/*     G = (A,B,C,D) from the factors Q = (AQR,BQ,CQR,DQ) and */
/*     R = (AQR,BR,CQR,DR) of its left coprime factorization */
/*                   -1 */
/*              G = R  * Q, */

/*     where G, Q and R are the corresponding transfer-function matrices. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. Also the number of rows of the */
/*             matrices B and BR and the number of columns of the matrix */
/*             C. N represents the order of the systems Q and R.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of input vector, i.e. the number of columns */
/*             of the matrices B and D.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of output vector, i.e. the number of rows of */
/*             the matrices C, D and DR and the number of columns of the */
/*             matrices BR and DR.  P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix AQR of the systems */
/*             Q and R. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the state dynamics matrix of the system G. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input/state matrix BQ of the system Q. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the input/state matrix of the system G. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix CQR of the systems */
/*             Q and R. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the state/output matrix of the system G. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the input/output matrix DQ of the system Q. */
/*             On exit, the leading P-by-M part of this array contains */
/*             the input/output matrix of the system G. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     BR      (input) DOUBLE PRECISION array, dimension (LDBR,P) */
/*             The leading N-by-P part of this array must contain the */
/*             input/state matrix BR of the system R. */

/*     LDBR    INTEGER */
/*             The leading dimension of array BR.  LDBR >= MAX(1,N). */

/*     DR      (input/output) DOUBLE PRECISION array, dimension (LDDR,P) */
/*             On entry, the leading P-by-P part of this array must */
/*             contain the input/output matrix DR of the system R. */
/*             On exit, the leading P-by-P part of this array contains */
/*             the LU factorization of the matrix DR, as computed by */
/*             LAPACK Library routine DGETRF. */

/*     LDDR    INTEGER */
/*             The leading dimension of array DR.  LDDR >= MAX(1,P). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (P) */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(1,4*P)) */
/*             On exit, DWORK(1) contains an estimate of the reciprocal */
/*             condition number of the matrix DR. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the matrix DR is singular; */
/*             = 2:  the matrix DR is numerically singular (warning); */
/*                   the calculations continued. */

/*     METHOD */

/*     The subroutine computes the matrices of the state-space */
/*     representation G = (A,B,C,D) by using the formulas: */

/*                      -1              -1 */
/*     A = AQR - BR * DR  * CQR,  C = DR  * CQR, */
/*                      -1              -1 */
/*     B = BQ  - BR * DR  * DQ,   D = DR  * DQ. */

/*     REFERENCES */

/*     [1] Varga A. */
/*         Coprime factors model reduction method based on */
/*         square-root balancing-free techniques. */
/*         System Analysis, Modelling and Simulation, */
/*         vol. 11, pp. 303-311, 1993. */

/*     CONTRIBUTOR */

/*     C. Oara and A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine LCFI. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */

/*     KEYWORDS */

/*     Coprime factorization, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 178 "SB08GD.f"
    /* Parameter adjustments */
#line 178 "SB08GD.f"
    a_dim1 = *lda;
#line 178 "SB08GD.f"
    a_offset = 1 + a_dim1;
#line 178 "SB08GD.f"
    a -= a_offset;
#line 178 "SB08GD.f"
    b_dim1 = *ldb;
#line 178 "SB08GD.f"
    b_offset = 1 + b_dim1;
#line 178 "SB08GD.f"
    b -= b_offset;
#line 178 "SB08GD.f"
    c_dim1 = *ldc;
#line 178 "SB08GD.f"
    c_offset = 1 + c_dim1;
#line 178 "SB08GD.f"
    c__ -= c_offset;
#line 178 "SB08GD.f"
    d_dim1 = *ldd;
#line 178 "SB08GD.f"
    d_offset = 1 + d_dim1;
#line 178 "SB08GD.f"
    d__ -= d_offset;
#line 178 "SB08GD.f"
    br_dim1 = *ldbr;
#line 178 "SB08GD.f"
    br_offset = 1 + br_dim1;
#line 178 "SB08GD.f"
    br -= br_offset;
#line 178 "SB08GD.f"
    dr_dim1 = *lddr;
#line 178 "SB08GD.f"
    dr_offset = 1 + dr_dim1;
#line 178 "SB08GD.f"
    dr -= dr_offset;
#line 178 "SB08GD.f"
    --iwork;
#line 178 "SB08GD.f"
    --dwork;
#line 178 "SB08GD.f"

#line 178 "SB08GD.f"
    /* Function Body */
#line 178 "SB08GD.f"
    *info = 0;

/*     Check the scalar input parameters. */

#line 182 "SB08GD.f"
    if (*n < 0) {
#line 183 "SB08GD.f"
	*info = -1;
#line 184 "SB08GD.f"
    } else if (*m < 0) {
#line 185 "SB08GD.f"
	*info = -2;
#line 186 "SB08GD.f"
    } else if (*p < 0) {
#line 187 "SB08GD.f"
	*info = -3;
#line 188 "SB08GD.f"
    } else if (*lda < max(1,*n)) {
#line 189 "SB08GD.f"
	*info = -5;
#line 190 "SB08GD.f"
    } else if (*ldb < max(1,*n)) {
#line 191 "SB08GD.f"
	*info = -7;
#line 192 "SB08GD.f"
    } else if (*ldc < max(1,*p)) {
#line 193 "SB08GD.f"
	*info = -9;
#line 194 "SB08GD.f"
    } else if (*ldd < max(1,*p)) {
#line 195 "SB08GD.f"
	*info = -11;
#line 196 "SB08GD.f"
    } else if (*ldbr < max(1,*n)) {
#line 197 "SB08GD.f"
	*info = -13;
#line 198 "SB08GD.f"
    } else if (*lddr < max(1,*p)) {
#line 199 "SB08GD.f"
	*info = -15;
#line 200 "SB08GD.f"
    }
#line 201 "SB08GD.f"
    if (*info != 0) {

/*        Error return. */

#line 205 "SB08GD.f"
	i__1 = -(*info);
#line 205 "SB08GD.f"
	xerbla_("SB08GD", &i__1, (ftnlen)6);
#line 206 "SB08GD.f"
	return 0;
#line 207 "SB08GD.f"
    }

/*     Quick return if possible. */

#line 211 "SB08GD.f"
    if (*p == 0) {
#line 212 "SB08GD.f"
	dwork[1] = 1.;
#line 213 "SB08GD.f"
	return 0;
#line 214 "SB08GD.f"
    }

/*     Factor the matrix  DR.  First, compute the 1-norm. */

#line 218 "SB08GD.f"
    drnorm = dlange_("1-norm", p, p, &dr[dr_offset], lddr, &dwork[1], (ftnlen)
	    6);
#line 219 "SB08GD.f"
    dgetrf_(p, p, &dr[dr_offset], lddr, &iwork[1], info);
#line 220 "SB08GD.f"
    if (*info != 0) {
#line 221 "SB08GD.f"
	*info = 1;
#line 222 "SB08GD.f"
	dwork[1] = 0.;
#line 223 "SB08GD.f"
	return 0;
#line 224 "SB08GD.f"
    }
/*                   -1 */
/*     Compute C = DR  * CQR. */

#line 228 "SB08GD.f"
    dgetrs_("NoTranspose", p, n, &dr[dr_offset], lddr, &iwork[1], &c__[
	    c_offset], ldc, info, (ftnlen)11);
/*                              -1 */
/*     Compute A = AQR - BR * DR  * CQR. */

#line 232 "SB08GD.f"
    dgemm_("NoTranspose", "NoTranspose", n, n, p, &c_b7, &br[br_offset], ldbr,
	     &c__[c_offset], ldc, &c_b8, &a[a_offset], lda, (ftnlen)11, (
	    ftnlen)11);
/*                   -1 */
/*     Compute D = DR  * DQ. */

#line 237 "SB08GD.f"
    dgetrs_("NoTranspose", p, m, &dr[dr_offset], lddr, &iwork[1], &d__[
	    d_offset], ldd, info, (ftnlen)11);
/*                             -1 */
/*     Compute B = BQ - BR * DR  * DQ. */

#line 241 "SB08GD.f"
    dgemm_("NoTranspose", "NoTranspose", n, m, p, &c_b7, &br[br_offset], ldbr,
	     &d__[d_offset], ldd, &c_b8, &b[b_offset], ldb, (ftnlen)11, (
	    ftnlen)11);

/*     Estimate the reciprocal condition number of DR. */
/*     Workspace  4*P. */

#line 247 "SB08GD.f"
    dgecon_("1-norm", p, &dr[dr_offset], lddr, &drnorm, &rcond, &dwork[1], &
	    iwork[1], info, (ftnlen)6);
#line 249 "SB08GD.f"
    if (rcond <= dlamch_("Epsilon", (ftnlen)7)) {
#line 249 "SB08GD.f"
	*info = 2;
#line 249 "SB08GD.f"
    }

#line 252 "SB08GD.f"
    dwork[1] = rcond;

#line 254 "SB08GD.f"
    return 0;
/* *** Last line of SB08GD *** */
} /* sb08gd_ */

