#line 1 "SB10WD.f"
/* SB10WD.f -- translated by f2c (version 20100827).
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

#line 1 "SB10WD.f"
/* Table of constant values */

static doublereal c_b5 = 1.;
static doublereal c_b6 = 0.;
static doublereal c_b22 = -1.;

/* Subroutine */ int sb10wd_(integer *n, integer *m, integer *np, integer *
	ncon, integer *nmeas, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *d__, integer 
	*ldd, doublereal *f, integer *ldf, doublereal *h__, integer *ldh, 
	doublereal *tu, integer *ldtu, doublereal *ty, integer *ldty, 
	doublereal *ak, integer *ldak, doublereal *bk, integer *ldbk, 
	doublereal *ck, integer *ldck, doublereal *dk, integer *lddk, integer 
	*info)
{
    /* System generated locals */
    integer a_dim1, a_offset, ak_dim1, ak_offset, b_dim1, b_offset, bk_dim1, 
	    bk_offset, c_dim1, c_offset, ck_dim1, ck_offset, d_dim1, d_offset,
	     dk_dim1, dk_offset, f_dim1, f_offset, h_dim1, h_offset, tu_dim1, 
	    tu_offset, ty_dim1, ty_offset, i__1;

    /* Local variables */
    static integer m1, m2, np1, np2;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);


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

/*     To compute the matrices of the H2 optimal controller */

/*              | AK | BK | */
/*          K = |----|----|, */
/*              | CK | DK | */

/*     from the state feedback matrix F and output injection matrix H as */
/*     determined by the SLICOT Library routine SB10VD. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     NP      (input) INTEGER */
/*             The row size of the matrix C.  NP >= 0. */

/*     NCON    (input) INTEGER */
/*             The number of control inputs (M2).  M >= NCON >= 0. */
/*             NP-NMEAS >= NCON. */

/*     NMEAS   (input) INTEGER */
/*             The number of measurements (NP2).  NP >= NMEAS >= 0. */
/*             M-NCON >= NMEAS. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B. Only the submatrix */
/*             B2 = B(:,M-M2+1:M) is used. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading NP-by-N part of this array must contain the */
/*             system output matrix C. Only the submatrix */
/*             C2 = C(NP-NP2+1:NP,:) is used. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,NP). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading NP-by-M part of this array must contain the */
/*             system input/output matrix D. Only the submatrix */
/*             D22 = D(NP-NP2+1:NP,M-M2+1:M) is used. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D.  LDD >= max(1,NP). */

/*     F       (input) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading NCON-by-N part of this array must contain the */
/*             state feedback matrix F. */

/*     LDF     INTEGER */
/*             The leading dimension of the array F.  LDF >= max(1,NCON). */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,NMEAS) */
/*             The leading N-by-NMEAS part of this array must contain the */
/*             output injection matrix H. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

/*     TU      (input) DOUBLE PRECISION array, dimension (LDTU,M2) */
/*             The leading M2-by-M2 part of this array must contain the */
/*             control transformation matrix TU, as obtained by the */
/*             SLICOT Library routine SB10UD. */

/*     LDTU    INTEGER */
/*             The leading dimension of the array TU.  LDTU >= max(1,M2). */

/*     TY      (input) DOUBLE PRECISION array, dimension (LDTY,NP2) */
/*             The leading NP2-by-NP2 part of this array must contain the */
/*             measurement transformation matrix TY, as obtained by the */
/*             SLICOT Library routine SB10UD. */

/*     LDTY    INTEGER */
/*             The leading dimension of the array TY. */
/*             LDTY >= max(1,NP2). */

/*     AK      (output) DOUBLE PRECISION array, dimension (LDAK,N) */
/*             The leading N-by-N part of this array contains the */
/*             controller state matrix AK. */

/*     LDAK    INTEGER */
/*             The leading dimension of the array AK.  LDAK >= max(1,N). */

/*     BK      (output) DOUBLE PRECISION array, dimension (LDBK,NMEAS) */
/*             The leading N-by-NMEAS part of this array contains the */
/*             controller input matrix BK. */

/*     LDBK    INTEGER */
/*             The leading dimension of the array BK.  LDBK >= max(1,N). */

/*     CK      (output) DOUBLE PRECISION array, dimension (LDCK,N) */
/*             The leading NCON-by-N part of this array contains the */
/*             controller output matrix CK. */

/*     LDCK    INTEGER */
/*             The leading dimension of the array CK. */
/*             LDCK >= max(1,NCON). */

/*     DK      (output) DOUBLE PRECISION array, dimension (LDDK,NMEAS) */
/*             The leading NCON-by-NMEAS part of this array contains the */
/*             controller input/output matrix DK. */

/*     LDDK    INTEGER */
/*             The leading dimension of the array DK. */
/*             LDDK >= max(1,NCON). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine implements the formulas given in [1], [2]. */

/*     REFERENCES */

/*     [1] Zhou, K., Doyle, J.C., and Glover, K. */
/*         Robust and Optimal Control. */
/*         Prentice-Hall, Upper Saddle River, NJ, 1996. */

/*     [2] Balas, G.J., Doyle, J.C., Glover, K., Packard, A., and */
/*         Smith, R. */
/*         mu-Analysis and Synthesis Toolbox. */
/*         The MathWorks Inc., Natick, Mass., 1995. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of the result depends on the condition numbers of the */
/*     input and output transformations. */

/*     CONTRIBUTORS */

/*     P.Hr. Petkov, D.W. Gu and M.M. Konstantinov, October 1998. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 1999. */

/*     KEYWORDS */

/*     Algebraic Riccati equation, H2 optimal control, LQG, LQR, optimal */
/*     regulator, robust control. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode and Test input parameters. */

#line 217 "SB10WD.f"
    /* Parameter adjustments */
#line 217 "SB10WD.f"
    a_dim1 = *lda;
#line 217 "SB10WD.f"
    a_offset = 1 + a_dim1;
#line 217 "SB10WD.f"
    a -= a_offset;
#line 217 "SB10WD.f"
    b_dim1 = *ldb;
#line 217 "SB10WD.f"
    b_offset = 1 + b_dim1;
#line 217 "SB10WD.f"
    b -= b_offset;
#line 217 "SB10WD.f"
    c_dim1 = *ldc;
#line 217 "SB10WD.f"
    c_offset = 1 + c_dim1;
#line 217 "SB10WD.f"
    c__ -= c_offset;
#line 217 "SB10WD.f"
    d_dim1 = *ldd;
#line 217 "SB10WD.f"
    d_offset = 1 + d_dim1;
#line 217 "SB10WD.f"
    d__ -= d_offset;
#line 217 "SB10WD.f"
    f_dim1 = *ldf;
#line 217 "SB10WD.f"
    f_offset = 1 + f_dim1;
#line 217 "SB10WD.f"
    f -= f_offset;
#line 217 "SB10WD.f"
    h_dim1 = *ldh;
#line 217 "SB10WD.f"
    h_offset = 1 + h_dim1;
#line 217 "SB10WD.f"
    h__ -= h_offset;
#line 217 "SB10WD.f"
    tu_dim1 = *ldtu;
#line 217 "SB10WD.f"
    tu_offset = 1 + tu_dim1;
#line 217 "SB10WD.f"
    tu -= tu_offset;
#line 217 "SB10WD.f"
    ty_dim1 = *ldty;
#line 217 "SB10WD.f"
    ty_offset = 1 + ty_dim1;
#line 217 "SB10WD.f"
    ty -= ty_offset;
#line 217 "SB10WD.f"
    ak_dim1 = *ldak;
#line 217 "SB10WD.f"
    ak_offset = 1 + ak_dim1;
#line 217 "SB10WD.f"
    ak -= ak_offset;
#line 217 "SB10WD.f"
    bk_dim1 = *ldbk;
#line 217 "SB10WD.f"
    bk_offset = 1 + bk_dim1;
#line 217 "SB10WD.f"
    bk -= bk_offset;
#line 217 "SB10WD.f"
    ck_dim1 = *ldck;
#line 217 "SB10WD.f"
    ck_offset = 1 + ck_dim1;
#line 217 "SB10WD.f"
    ck -= ck_offset;
#line 217 "SB10WD.f"
    dk_dim1 = *lddk;
#line 217 "SB10WD.f"
    dk_offset = 1 + dk_dim1;
#line 217 "SB10WD.f"
    dk -= dk_offset;
#line 217 "SB10WD.f"

#line 217 "SB10WD.f"
    /* Function Body */
#line 217 "SB10WD.f"
    m1 = *m - *ncon;
#line 218 "SB10WD.f"
    m2 = *ncon;
#line 219 "SB10WD.f"
    np1 = *np - *nmeas;
#line 220 "SB10WD.f"
    np2 = *nmeas;

#line 222 "SB10WD.f"
    *info = 0;
#line 223 "SB10WD.f"
    if (*n < 0) {
#line 224 "SB10WD.f"
	*info = -1;
#line 225 "SB10WD.f"
    } else if (*m < 0) {
#line 226 "SB10WD.f"
	*info = -2;
#line 227 "SB10WD.f"
    } else if (*np < 0) {
#line 228 "SB10WD.f"
	*info = -3;
#line 229 "SB10WD.f"
    } else if (*ncon < 0 || m1 < 0 || m2 > np1) {
#line 230 "SB10WD.f"
	*info = -4;
#line 231 "SB10WD.f"
    } else if (*nmeas < 0 || np1 < 0 || np2 > m1) {
#line 232 "SB10WD.f"
	*info = -5;
#line 233 "SB10WD.f"
    } else if (*lda < max(1,*n)) {
#line 234 "SB10WD.f"
	*info = -7;
#line 235 "SB10WD.f"
    } else if (*ldb < max(1,*n)) {
#line 236 "SB10WD.f"
	*info = -9;
#line 237 "SB10WD.f"
    } else if (*ldc < max(1,*np)) {
#line 238 "SB10WD.f"
	*info = -11;
#line 239 "SB10WD.f"
    } else if (*ldd < max(1,*np)) {
#line 240 "SB10WD.f"
	*info = -13;
#line 241 "SB10WD.f"
    } else if (*ldf < max(1,m2)) {
#line 242 "SB10WD.f"
	*info = -15;
#line 243 "SB10WD.f"
    } else if (*ldh < max(1,*n)) {
#line 244 "SB10WD.f"
	*info = -17;
#line 245 "SB10WD.f"
    } else if (*ldtu < max(1,m2)) {
#line 246 "SB10WD.f"
	*info = -19;
#line 247 "SB10WD.f"
    } else if (*ldty < max(1,np2)) {
#line 248 "SB10WD.f"
	*info = -21;
#line 249 "SB10WD.f"
    } else if (*ldak < max(1,*n)) {
#line 250 "SB10WD.f"
	*info = -23;
#line 251 "SB10WD.f"
    } else if (*ldbk < max(1,*n)) {
#line 252 "SB10WD.f"
	*info = -25;
#line 253 "SB10WD.f"
    } else if (*ldck < max(1,m2)) {
#line 254 "SB10WD.f"
	*info = -27;
#line 255 "SB10WD.f"
    } else if (*lddk < max(1,m2)) {
#line 256 "SB10WD.f"
	*info = -29;
#line 257 "SB10WD.f"
    }
#line 258 "SB10WD.f"
    if (*info != 0) {
#line 259 "SB10WD.f"
	i__1 = -(*info);
#line 259 "SB10WD.f"
	xerbla_("SB10WD", &i__1, (ftnlen)6);
#line 260 "SB10WD.f"
	return 0;
#line 261 "SB10WD.f"
    }

/*     Quick return if possible. */

#line 265 "SB10WD.f"
    if (*n == 0 || *m == 0 || *np == 0 || m1 == 0 || m2 == 0 || np1 == 0 || 
	    np2 == 0) {
#line 265 "SB10WD.f"
	return 0;
#line 265 "SB10WD.f"
    }

/*     Compute the transpose of D22*F . BK is used as workspace. */

#line 270 "SB10WD.f"
    dgemm_("T", "T", n, &np2, &m2, &c_b5, &f[f_offset], ldf, &d__[np1 + 1 + (
	    m1 + 1) * d_dim1], ldd, &c_b6, &bk[bk_offset], ldbk, (ftnlen)1, (
	    ftnlen)1);

/*     Find AK = A + H*C2 + B2*F + H*D22*F . */

#line 275 "SB10WD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &ak[ak_offset], ldak, (ftnlen)4);
#line 276 "SB10WD.f"
    dgemm_("N", "N", n, n, &np2, &c_b5, &h__[h_offset], ldh, &c__[np1 + 1 + 
	    c_dim1], ldc, &c_b5, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);
#line 278 "SB10WD.f"
    dgemm_("N", "N", n, n, &m2, &c_b5, &b[(m1 + 1) * b_dim1 + 1], ldb, &f[
	    f_offset], ldf, &c_b5, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1)
	    ;
#line 280 "SB10WD.f"
    dgemm_("N", "T", n, n, &np2, &c_b5, &h__[h_offset], ldh, &bk[bk_offset], 
	    ldbk, &c_b5, &ak[ak_offset], ldak, (ftnlen)1, (ftnlen)1);

/*     Find BK = -H*Ty . */

#line 285 "SB10WD.f"
    dgemm_("N", "N", n, &np2, &np2, &c_b22, &h__[h_offset], ldh, &ty[
	    ty_offset], ldty, &c_b6, &bk[bk_offset], ldbk, (ftnlen)1, (ftnlen)
	    1);

/*     Find CK = Tu*F . */

#line 290 "SB10WD.f"
    dgemm_("N", "N", &m2, n, &m2, &c_b5, &tu[tu_offset], ldtu, &f[f_offset], 
	    ldf, &c_b6, &ck[ck_offset], ldck, (ftnlen)1, (ftnlen)1);

/*     Find DK . */

#line 295 "SB10WD.f"
    dlaset_("Full", &m2, &np2, &c_b6, &c_b6, &dk[dk_offset], lddk, (ftnlen)4);

#line 297 "SB10WD.f"
    return 0;
/* *** Last line of SB10WD *** */
} /* sb10wd_ */

