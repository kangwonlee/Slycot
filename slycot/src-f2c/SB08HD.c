#line 1 "SB08HD.f"
/* SB08HD.f -- translated by f2c (version 20100827).
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

#line 1 "SB08HD.f"
/* Table of constant values */

static doublereal c_b8 = 1.;
static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b18 = -1.;

/* Subroutine */ int sb08hd_(integer *n, integer *m, integer *p, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *cr, integer *
	ldcr, doublereal *dr, integer *lddr, integer *iwork, doublereal *
	dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, cr_dim1, 
	    cr_offset, d_dim1, d_offset, dr_dim1, dr_offset, i__1;

    /* Local variables */
    extern /* Subroutine */ int ma02gd_(integer *, doublereal *, integer *, 
	    integer *, integer *, integer *, integer *), dgemm_(char *, char *
	    , integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dgetrf_(integer *, integer *, doublereal *, 
	    integer *, integer *, integer *), xerbla_(char *, integer *, 
	    ftnlen);
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
/*     G = (A,B,C,D) from the factors Q = (AQR,BQR,CQ,DQ) and */
/*     R = (AQR,BQR,CR,DR) of its right coprime factorization */
/*                       -1 */
/*              G = Q * R  , */

/*     where G, Q and R are the corresponding transfer-function matrices. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. Also the number of rows of the */
/*             matrix B and the number of columns of the matrices C and */
/*             CR. N represents the order of the systems Q and R. */
/*             N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of input vector. Also the number of columns */
/*             of the matrices B, D and DR and the number of rows of the */
/*             matrices CR and DR.  M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of output vector. Also the number of rows */
/*             of the matrices C and D.  P >= 0. */

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
/*             contain the input/state matrix BQR of the systems Q and R. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the input/state matrix of the system G. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix CQ of the system Q. */
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

/*     CR      (input) DOUBLE PRECISION array, dimension (LDCR,N) */
/*             The leading M-by-N part of this array must contain the */
/*             state/output matrix CR of the system R. */

/*     LDCR    INTEGER */
/*             The leading dimension of array CR.  LDCR >= MAX(1,M). */

/*     DR      (input/output) DOUBLE PRECISION array, dimension (LDDR,M) */
/*             On entry, the leading M-by-M part of this array must */
/*             contain the input/output matrix DR of the system R. */
/*             On exit, the leading M-by-M part of this array contains */
/*             the LU factorization of the matrix DR, as computed by */
/*             LAPACK Library routine DGETRF. */

/*     LDDR    INTEGER */
/*             The leading dimension of array DR.  LDDR >= MAX(1,M). */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(1,4*M)) */
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

/*                       -1                   -1 */
/*     A = AQR - BQR * DR  * CR,  B = BQR * DR  , */
/*                      -1                   -1 */
/*     C = CQ  - DQ * DR  * CR,   D = DQ * DR  . */

/*     REFERENCES */

/*     [1] Varga A. */
/*         Coprime factors model reduction method based on */
/*         square-root balancing-free techniques. */
/*         System Analysis, Modelling and Simulation, */
/*         vol. 11, pp. 303-311, 1993. */

/*     CONTRIBUTOR */

/*     C. Oara and A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine RCFI. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Nov. 1998, */
/*     full BLAS 3 version. */

/*     REVISIONS */

/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     Mar. 2000, V. Sima, Research Institute for Informatics, Bucharest. */

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

#line 181 "SB08HD.f"
    /* Parameter adjustments */
#line 181 "SB08HD.f"
    a_dim1 = *lda;
#line 181 "SB08HD.f"
    a_offset = 1 + a_dim1;
#line 181 "SB08HD.f"
    a -= a_offset;
#line 181 "SB08HD.f"
    b_dim1 = *ldb;
#line 181 "SB08HD.f"
    b_offset = 1 + b_dim1;
#line 181 "SB08HD.f"
    b -= b_offset;
#line 181 "SB08HD.f"
    c_dim1 = *ldc;
#line 181 "SB08HD.f"
    c_offset = 1 + c_dim1;
#line 181 "SB08HD.f"
    c__ -= c_offset;
#line 181 "SB08HD.f"
    d_dim1 = *ldd;
#line 181 "SB08HD.f"
    d_offset = 1 + d_dim1;
#line 181 "SB08HD.f"
    d__ -= d_offset;
#line 181 "SB08HD.f"
    cr_dim1 = *ldcr;
#line 181 "SB08HD.f"
    cr_offset = 1 + cr_dim1;
#line 181 "SB08HD.f"
    cr -= cr_offset;
#line 181 "SB08HD.f"
    dr_dim1 = *lddr;
#line 181 "SB08HD.f"
    dr_offset = 1 + dr_dim1;
#line 181 "SB08HD.f"
    dr -= dr_offset;
#line 181 "SB08HD.f"
    --iwork;
#line 181 "SB08HD.f"
    --dwork;
#line 181 "SB08HD.f"

#line 181 "SB08HD.f"
    /* Function Body */
#line 181 "SB08HD.f"
    *info = 0;

/*     Check the scalar input parameters. */

#line 185 "SB08HD.f"
    if (*n < 0) {
#line 186 "SB08HD.f"
	*info = -1;
#line 187 "SB08HD.f"
    } else if (*m < 0) {
#line 188 "SB08HD.f"
	*info = -2;
#line 189 "SB08HD.f"
    } else if (*p < 0) {
#line 190 "SB08HD.f"
	*info = -3;
#line 191 "SB08HD.f"
    } else if (*lda < max(1,*n)) {
#line 192 "SB08HD.f"
	*info = -5;
#line 193 "SB08HD.f"
    } else if (*ldb < max(1,*n)) {
#line 194 "SB08HD.f"
	*info = -7;
#line 195 "SB08HD.f"
    } else if (*ldc < max(1,*p)) {
#line 196 "SB08HD.f"
	*info = -9;
#line 197 "SB08HD.f"
    } else if (*ldd < max(1,*p)) {
#line 198 "SB08HD.f"
	*info = -11;
#line 199 "SB08HD.f"
    } else if (*ldcr < max(1,*m)) {
#line 200 "SB08HD.f"
	*info = -13;
#line 201 "SB08HD.f"
    } else if (*lddr < max(1,*m)) {
#line 202 "SB08HD.f"
	*info = -15;
#line 203 "SB08HD.f"
    }
#line 204 "SB08HD.f"
    if (*info != 0) {

/*        Error return. */

#line 208 "SB08HD.f"
	i__1 = -(*info);
#line 208 "SB08HD.f"
	xerbla_("SB08HD", &i__1, (ftnlen)6);
#line 209 "SB08HD.f"
	return 0;
#line 210 "SB08HD.f"
    }

/*     Quick return if possible. */

#line 214 "SB08HD.f"
    if (*m == 0) {
#line 215 "SB08HD.f"
	dwork[1] = 1.;
#line 216 "SB08HD.f"
	return 0;
#line 217 "SB08HD.f"
    }

/*     Factor the matrix  DR.  First, compute the 1-norm. */

#line 221 "SB08HD.f"
    drnorm = dlange_("1-norm", m, m, &dr[dr_offset], lddr, &dwork[1], (ftnlen)
	    6);
#line 222 "SB08HD.f"
    dgetrf_(m, m, &dr[dr_offset], lddr, &iwork[1], info);
#line 223 "SB08HD.f"
    if (*info != 0) {
#line 224 "SB08HD.f"
	*info = 1;
#line 225 "SB08HD.f"
	dwork[1] = 0.;
#line 226 "SB08HD.f"
	return 0;
#line 227 "SB08HD.f"
    }
/*                         -1 */
/*     Compute B = BQR * DR  , using the factorization P*DR = L*U. */

#line 231 "SB08HD.f"
    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", n, m, &c_b8, &dr[
	    dr_offset], lddr, &b[b_offset], ldb, (ftnlen)5, (ftnlen)5, (
	    ftnlen)11, (ftnlen)7);
#line 233 "SB08HD.f"
    dtrsm_("Right", "Lower", "NoTranspose", "Unit", n, m, &c_b8, &dr[
	    dr_offset], lddr, &b[b_offset], ldb, (ftnlen)5, (ftnlen)5, (
	    ftnlen)11, (ftnlen)4);
#line 235 "SB08HD.f"
    ma02gd_(n, &b[b_offset], ldb, &c__1, m, &iwork[1], &c_n1);
/*                               -1 */
/*     Compute A = AQR - BQR * DR  * CR. */

#line 239 "SB08HD.f"
    dgemm_("NoTranspose", "NoTranspose", n, n, m, &c_b18, &b[b_offset], ldb, &
	    cr[cr_offset], ldcr, &c_b8, &a[a_offset], lda, (ftnlen)11, (
	    ftnlen)11);
/*                        -1 */
/*     Compute D = DQ * DR  . */

#line 244 "SB08HD.f"
    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", p, m, &c_b8, &dr[
	    dr_offset], lddr, &d__[d_offset], ldd, (ftnlen)5, (ftnlen)5, (
	    ftnlen)11, (ftnlen)7);
#line 246 "SB08HD.f"
    dtrsm_("Right", "Lower", "NoTranspose", "Unit", p, m, &c_b8, &dr[
	    dr_offset], lddr, &d__[d_offset], ldd, (ftnlen)5, (ftnlen)5, (
	    ftnlen)11, (ftnlen)4);
#line 248 "SB08HD.f"
    ma02gd_(p, &d__[d_offset], ldd, &c__1, m, &iwork[1], &c_n1);
/*                             -1 */
/*     Compute C = CQ - DQ * DR  * CR. */

#line 252 "SB08HD.f"
    dgemm_("NoTranspose", "NoTranspose", p, n, m, &c_b18, &d__[d_offset], ldd,
	     &cr[cr_offset], ldcr, &c_b8, &c__[c_offset], ldc, (ftnlen)11, (
	    ftnlen)11);

/*     Estimate the reciprocal condition number of DR. */
/*     Workspace  4*M. */

#line 258 "SB08HD.f"
    dgecon_("1-norm", m, &dr[dr_offset], lddr, &drnorm, &rcond, &dwork[1], &
	    iwork[1], info, (ftnlen)6);
#line 260 "SB08HD.f"
    if (rcond <= dlamch_("Epsilon", (ftnlen)7)) {
#line 260 "SB08HD.f"
	*info = 2;
#line 260 "SB08HD.f"
    }

#line 263 "SB08HD.f"
    dwork[1] = rcond;

#line 265 "SB08HD.f"
    return 0;
/* *** Last line of SB08HD *** */
} /* sb08hd_ */
