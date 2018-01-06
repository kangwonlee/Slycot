#line 1 "NF01BW.f"
/* NF01BW.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BW.f"
/* Table of constant values */

static doublereal c_b4 = 1.;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int nf01bw_(integer *n, integer *ipar, integer *lipar, 
	doublereal *dpar, integer *ldpar, doublereal *j, integer *ldj, 
	doublereal *x, integer *incx, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1, i__2;

    /* Local variables */
    static doublereal c__;
    static integer m, bn, jl, ix, xl, st, bsm, bsn, ibsm, ibsn, nths;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dgemv_(char *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen), dcopy_(integer *, doublereal *, 
	    integer *, doublereal *, integer *), xerbla_(char *, integer *, 
	    ftnlen);


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

/*     To compute the matrix-vector product x <-- (J'*J + c*I)*x, for the */
/*     Jacobian J as received from SLICOT Library routine NF01BD: */

/*          /  dy(1)/dwb(1)  |  dy(1)/dtheta  \ */
/*     Jc = |       :        |       :        | . */
/*          \  dy(L)/dwb(L)  |  dy(L)/dtheta  / */

/*     This is a compressed representation of the actual structure */

/*         /   J_1    0    ..   0   |  L_1  \ */
/*         |    0    J_2   ..   0   |  L_2  | */
/*     J = |    :     :    ..   :   |   :   | . */
/*         |    :     :    ..   :   |   :   | */
/*         \    0     0    ..  J_L  |  L_L  / */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of the vector x. */
/*             N = BN*BSN + ST >= 0.  (See parameter description below.) */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix J, as follows: */
/*             IPAR(1) must contain ST, the number of parameters */
/*                     corresponding to the linear part.  ST >= 0. */
/*             IPAR(2) must contain BN, the number of blocks, BN = L, */
/*                     for the parameters corresponding to the nonlinear */
/*                     part.  BN >= 0. */
/*             IPAR(3) must contain BSM, the number of rows of the blocks */
/*                     J_k = dy(k)/dwb(k), k = 1:BN, if BN > 0, or the */
/*                     number of rows of the matrix J, if BN <= 1. */
/*             IPAR(4) must contain BSN, the number of columns of the */
/*                     blocks J_k, k = 1:BN.  BSN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 4. */

/*     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*             The real parameters needed for solving the problem. */
/*             The entry DPAR(1) must contain the real scalar c. */

/*     LDPAR   (input) INTEGER */
/*             The length of the array DPAR.  LDPAR >= 1. */

/*     J       (input) DOUBLE PRECISION array, dimension (LDJ, NC) */
/*             where NC = N if BN <= 1, and NC = BSN+ST, if BN > 1. */
/*             The leading NR-by-NC part of this array must contain */
/*             the (compressed) representation (Jc) of the Jacobian */
/*             matrix J, where NR = BSM if BN <= 1, and NR = BN*BSM, */
/*             if BN > 1. */

/*     LDJ     (input) INTEGER */
/*             The leading dimension of array J.  LDJ >= MAX(1,NR). */

/*     X       (input/output) DOUBLE PRECISION array, dimension */
/*             (1+(N-1)*INCX) */
/*             On entry, this incremented array must contain the */
/*             vector x. */
/*             On exit, this incremented array contains the value of the */
/*             matrix-vector product (J'*J + c*I)*x. */

/*     INCX    (input) INTEGER */
/*             The increment for the elements of X.  INCX >= 1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= NR. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The associativity of matrix multiplications is used; the result */
/*     is obtained as:  x_out = J'*( J*x ) + c*x. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Mar. 2001, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001, */
/*     Mar. 2002. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations, */
/*     Wiener system. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 146 "NF01BW.f"
    /* Parameter adjustments */
#line 146 "NF01BW.f"
    --ipar;
#line 146 "NF01BW.f"
    --dpar;
#line 146 "NF01BW.f"
    j_dim1 = *ldj;
#line 146 "NF01BW.f"
    j_offset = 1 + j_dim1;
#line 146 "NF01BW.f"
    j -= j_offset;
#line 146 "NF01BW.f"
    --x;
#line 146 "NF01BW.f"
    --dwork;
#line 146 "NF01BW.f"

#line 146 "NF01BW.f"
    /* Function Body */
#line 146 "NF01BW.f"
    *info = 0;

#line 148 "NF01BW.f"
    if (*n < 0) {
#line 149 "NF01BW.f"
	*info = -1;
#line 150 "NF01BW.f"
    } else if (*lipar < 4) {
#line 151 "NF01BW.f"
	*info = -3;
#line 152 "NF01BW.f"
    } else if (*ldpar < 1) {
#line 153 "NF01BW.f"
	*info = -5;
#line 154 "NF01BW.f"
    } else if (*incx < 1) {
#line 155 "NF01BW.f"
	*info = -9;
#line 156 "NF01BW.f"
    } else {
#line 157 "NF01BW.f"
	st = ipar[1];
#line 158 "NF01BW.f"
	bn = ipar[2];
#line 159 "NF01BW.f"
	bsm = ipar[3];
#line 160 "NF01BW.f"
	bsn = ipar[4];
#line 161 "NF01BW.f"
	nths = bn * bsn;
#line 162 "NF01BW.f"
	if (bn > 1) {
#line 163 "NF01BW.f"
	    m = bn * bsm;
#line 164 "NF01BW.f"
	} else {
#line 165 "NF01BW.f"
	    m = bsm;
#line 166 "NF01BW.f"
	}
/* Computing MIN */
#line 167 "NF01BW.f"
	i__1 = min(st,bn), i__1 = min(i__1,bsm);
#line 167 "NF01BW.f"
	if (min(i__1,bsn) < 0) {
#line 168 "NF01BW.f"
	    *info = -2;
#line 169 "NF01BW.f"
	} else if (*n != nths + st) {
#line 170 "NF01BW.f"
	    *info = -1;
#line 171 "NF01BW.f"
	} else if (*ldj < max(1,m)) {
#line 172 "NF01BW.f"
	    *info = -7;
#line 173 "NF01BW.f"
	} else if (*ldwork < m) {
#line 174 "NF01BW.f"
	    *info = -11;
#line 175 "NF01BW.f"
	}
#line 176 "NF01BW.f"
    }

/*     Return if there are illegal arguments. */

#line 180 "NF01BW.f"
    if (*info != 0) {
#line 181 "NF01BW.f"
	i__1 = -(*info);
#line 181 "NF01BW.f"
	xerbla_("NF01BW", &i__1, (ftnlen)6);
#line 182 "NF01BW.f"
	return 0;
#line 183 "NF01BW.f"
    }

/*     Quick return if possible. */

#line 187 "NF01BW.f"
    if (*n == 0) {
#line 187 "NF01BW.f"
	return 0;
#line 187 "NF01BW.f"
    }

#line 190 "NF01BW.f"
    c__ = dpar[1];

#line 192 "NF01BW.f"
    if (m == 0) {

/*        Special case, void Jacobian: x <-- c*x. */

#line 196 "NF01BW.f"
	dscal_(n, &c__, &x[1], incx);
#line 197 "NF01BW.f"
	return 0;
#line 198 "NF01BW.f"
    }

#line 200 "NF01BW.f"
    if (bn <= 1 || bsn == 0) {

/*        Special case, l <= 1 or BSN = 0: the Jacobian is represented */
/*        as a full matrix. Adapted code from NF01BX is included in-line. */

#line 205 "NF01BW.f"
	dgemv_("NoTranspose", &m, n, &c_b4, &j[j_offset], ldj, &x[1], incx, &
		c_b5, &dwork[1], &c__1, (ftnlen)11);
#line 207 "NF01BW.f"
	dgemv_("Transpose", &m, n, &c_b4, &j[j_offset], ldj, &dwork[1], &c__1,
		 &c__, &x[1], incx, (ftnlen)9);
#line 209 "NF01BW.f"
	return 0;
#line 210 "NF01BW.f"
    }

/*     General case: l > 1, BSN > 0, BSM > 0. */

#line 214 "NF01BW.f"
    jl = bsn + 1;
#line 215 "NF01BW.f"
    ix = bsn * *incx;
#line 216 "NF01BW.f"
    xl = bn * ix + 1;

#line 218 "NF01BW.f"
    if (st > 0) {
#line 219 "NF01BW.f"
	dgemv_("NoTranspose", &m, &st, &c_b4, &j[jl * j_dim1 + 1], ldj, &x[xl]
		, incx, &c_b5, &dwork[1], &c__1, (ftnlen)11);
#line 221 "NF01BW.f"
    } else {
#line 222 "NF01BW.f"
	dwork[1] = 0.;
#line 223 "NF01BW.f"
	dcopy_(&m, &dwork[1], &c__0, &dwork[1], &c__1);
#line 224 "NF01BW.f"
    }
#line 225 "NF01BW.f"
    ibsn = 1;

#line 227 "NF01BW.f"
    i__1 = m;
#line 227 "NF01BW.f"
    i__2 = bsm;
#line 227 "NF01BW.f"
    for (ibsm = 1; i__2 < 0 ? ibsm >= i__1 : ibsm <= i__1; ibsm += i__2) {
#line 228 "NF01BW.f"
	dgemv_("NoTranspose", &bsm, &bsn, &c_b4, &j[ibsm + j_dim1], ldj, &x[
		ibsn], incx, &c_b4, &dwork[ibsm], &c__1, (ftnlen)11);
#line 230 "NF01BW.f"
	dgemv_("Transpose", &bsm, &bsn, &c_b4, &j[ibsm + j_dim1], ldj, &dwork[
		ibsm], &c__1, &c__, &x[ibsn], incx, (ftnlen)9);
#line 232 "NF01BW.f"
	ibsn += ix;
#line 233 "NF01BW.f"
/* L10: */
#line 233 "NF01BW.f"
    }

#line 235 "NF01BW.f"
    if (st > 0) {
#line 235 "NF01BW.f"
	dgemv_("Transpose", &m, &st, &c_b4, &j[jl * j_dim1 + 1], ldj, &dwork[
		1], &c__1, &c__, &x[xl], incx, (ftnlen)9);
#line 235 "NF01BW.f"
    }

#line 239 "NF01BW.f"
    return 0;

/* *** Last line of NF01BW *** */
} /* nf01bw_ */

