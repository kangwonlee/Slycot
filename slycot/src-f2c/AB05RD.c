#line 1 "AB05RD.f"
/* AB05RD.f -- translated by f2c (version 20100827).
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

#line 1 "AB05RD.f"
/* Table of constant values */

static doublereal c_b9 = 1.;
static doublereal c_b16 = 0.;

/* Subroutine */ int ab05rd_(char *fbtype, char *jobd, integer *n, integer *m,
	 integer *p, integer *mv, integer *pz, doublereal *alpha, doublereal *
	beta, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *f, integer *ldf, doublereal *k, integer *ldk, doublereal *
	g, integer *ldg, doublereal *h__, integer *ldh, doublereal *rcond, 
	doublereal *bc, integer *ldbc, doublereal *cc, integer *ldcc, 
	doublereal *dc, integer *lddc, integer *iwork, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen fbtype_len, ftnlen jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, bc_dim1, bc_offset, c_dim1, 
	    c_offset, cc_dim1, cc_offset, d_dim1, d_offset, dc_dim1, 
	    dc_offset, f_dim1, f_offset, g_dim1, g_offset, h_dim1, h_offset, 
	    k_dim1, k_offset, i__1, i__2;

    /* Local variables */
    static integer ldwp;
    extern /* Subroutine */ int ab05sd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen), dgemm_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen, ftnlen);
    static logical ljobd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical unitf, outpf;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To construct for a given state space system (A,B,C,D) the closed- */
/*     loop system (Ac,Bc,Cc,Dc) corresponding to the mixed output and */
/*     state feedback control law */

/*          u = alpha*F*y + beta*K*x + G*v */
/*          z = H*y. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FBTYPE  CHARACTER*1 */
/*             Specifies the type of the feedback law as follows: */
/*             = 'I':  Unitary output feedback (F = I); */
/*             = 'O':  General output feedback. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears */
/*             in the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The dimension of state vector x, i.e. the order of the */
/*             matrix A, the number of rows of B and the number of */
/*             columns of C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The dimension of input vector u, i.e. the number of */
/*             columns of matrices B and D, and the number of rows of F. */
/*             M >= 0. */

/*     P       (input) INTEGER */
/*             The dimension of output vector y, i.e. the number of rows */
/*             of matrices C and D, and the number of columns of F. */
/*             P >= 0 and P = M if FBTYPE = 'I'. */

/*     MV      (input) INTEGER */
/*             The dimension of the new input vector v, i.e. the number */
/*             of columns of matrix G.  MV >= 0. */

/*     PZ      (input) INTEGER. */
/*             The dimension of the new output vector z, i.e. the number */
/*             of rows of matrix H.  PZ >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The coefficient alpha in the output feedback law. */

/*     BETA    (input) DOUBLE PRECISION. */
/*             The coefficient beta in the state feedback law. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the system state transition matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the state matrix Ac of the closed-loop system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the system input matrix B. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the intermediary input matrix B1 (see METHOD). */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the system output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the intermediary output matrix C1+BETA*D1*K (see METHOD). */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,P) if N > 0. */
/*             LDC >= 1 if N = 0. */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, if JOBD = 'D', the leading P-by-M part of this */
/*             array must contain the system direct input/output */
/*             transmission matrix D. */
/*             On exit, the leading P-by-M part of this array contains */
/*             the intermediary direct input/output transmission matrix */
/*             D1 (see METHOD). */
/*             The array D is not referenced if JOBD = 'Z'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P) if JOBD = 'D'. */
/*             LDD >= 1 if JOBD = 'Z'. */

/*     F       (input) DOUBLE PRECISION array, dimension (LDF,P) */
/*             If FBTYPE = 'O', the leading M-by-P part of this array */
/*             must contain the output feedback matrix F. */
/*             If FBTYPE = 'I', then the feedback matrix is assumed to be */
/*             an M x M order identity matrix. */
/*             The array F is not referenced if FBTYPE = 'I'  or */
/*             ALPHA = 0. */

/*     LDF     INTEGER */
/*             The leading dimension of array F. */
/*             LDF >= MAX(1,M) if FBTYPE = 'O' and ALPHA <> 0. */
/*             LDF >= 1 if FBTYPE = 'I' or ALPHA = 0. */

/*     K       (input) DOUBLE PRECISION array, dimension (LDK,N) */
/*             The leading M-by-N part of this array must contain the */
/*             state feedback matrix K. */
/*             The array K is not referenced if BETA = 0. */

/*     LDK     INTEGER */
/*             The leading dimension of the array K. */
/*             LDK >= MAX(1,M) if BETA <> 0. */
/*             LDK >= 1 if BETA = 0. */

/*     G       (input) DOUBLE PRECISION array, dimension (LDG,MV) */
/*             The leading M-by-MV part of this array must contain the */
/*             system input scaling matrix G. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,M). */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,P) */
/*             The leading PZ-by-P part of this array must contain the */
/*             system output scaling matrix H. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= MAX(1,PZ). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The reciprocal condition number of the matrix */
/*             I - alpha*D*F. */

/*     BC      (output) DOUBLE PRECISION array, dimension (LDBC,MV) */
/*             The leading N-by-MV part of this array contains the input */
/*             matrix Bc of the closed-loop system. */

/*     LDBC    INTEGER */
/*             The leading dimension of array BC.  LDBC >= MAX(1,N). */

/*     CC      (output) DOUBLE PRECISION array, dimension (LDCC,N) */
/*             The leading PZ-by-N part of this array contains the */
/*             system output matrix Cc of the closed-loop system. */

/*     LDCC    INTEGER */
/*             The leading dimension of array CC. */
/*             LDCC >= MAX(1,PZ) if N > 0. */
/*             LDCC >= 1 if N = 0. */

/*     DC      (output) DOUBLE PRECISION array, dimension (LDDC,MV) */
/*             If JOBD = 'D', the leading PZ-by-MV part of this array */
/*             contains the direct input/output transmission matrix Dc */
/*             of the closed-loop system. */
/*             The array DC is not referenced if JOBD = 'Z'. */

/*     LDDC    INTEGER */
/*             The leading dimension of array DC. */
/*             LDDC >= MAX(1,PZ) if JOBD = 'D'. */
/*             LDDC >= 1 if JOBD = 'Z'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= MAX(1,2*P) if JOBD = 'D'. */
/*             LIWORK >= 1 if JOBD = 'Z'. */
/*             IWORK is not referenced if JOBD = 'Z'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= wspace, where */
/*                   wspace = MAX( 1, M, P*MV, P*P + 4*P ) if JOBD = 'D', */
/*                   wspace = MAX( 1, M ) if JOBD = 'Z'. */
/*             For best performance, LDWORK >= MAX( wspace, N*M, N*P ). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix I - alpha*D*F is numerically singular. */

/*     METHOD */

/*     The matrices of the closed-loop system have the expressions: */

/*     Ac = A1 + beta*B1*K,      Bc = B1*G, */
/*     Cc = H*(C1 + beta*D1*K),  Dc = H*D1*G, */

/*     where */

/*     A1 = A + alpha*B*F*E*C,   B1 = B + alpha*B*F*E*D, */
/*     C1 = E*C,                 D1 = E*D, */

/*     with E = (I - alpha*D*F)**-1. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of computations basically depends on the conditioning */
/*     of the matrix I - alpha*D*F. If RCOND is very small, it is likely */
/*     that the computed results are inaccurate. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Research Establishment, */
/*     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven, */
/*     Belgium, Nov. 1996. */

/*     REVISIONS */

/*     January 14, 1997, February 18, 1998. */
/*     V. Sima, Research Institute for Informatics, Bucharest, July 2003, */
/*     Jan. 2005. */

/*     KEYWORDS */

/*     Multivariable system, state-space model, state-space */
/*     representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External functions .. */
/*     .. External subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the input scalar arguments. */

#line 281 "AB05RD.f"
    /* Parameter adjustments */
#line 281 "AB05RD.f"
    a_dim1 = *lda;
#line 281 "AB05RD.f"
    a_offset = 1 + a_dim1;
#line 281 "AB05RD.f"
    a -= a_offset;
#line 281 "AB05RD.f"
    b_dim1 = *ldb;
#line 281 "AB05RD.f"
    b_offset = 1 + b_dim1;
#line 281 "AB05RD.f"
    b -= b_offset;
#line 281 "AB05RD.f"
    c_dim1 = *ldc;
#line 281 "AB05RD.f"
    c_offset = 1 + c_dim1;
#line 281 "AB05RD.f"
    c__ -= c_offset;
#line 281 "AB05RD.f"
    d_dim1 = *ldd;
#line 281 "AB05RD.f"
    d_offset = 1 + d_dim1;
#line 281 "AB05RD.f"
    d__ -= d_offset;
#line 281 "AB05RD.f"
    f_dim1 = *ldf;
#line 281 "AB05RD.f"
    f_offset = 1 + f_dim1;
#line 281 "AB05RD.f"
    f -= f_offset;
#line 281 "AB05RD.f"
    k_dim1 = *ldk;
#line 281 "AB05RD.f"
    k_offset = 1 + k_dim1;
#line 281 "AB05RD.f"
    k -= k_offset;
#line 281 "AB05RD.f"
    g_dim1 = *ldg;
#line 281 "AB05RD.f"
    g_offset = 1 + g_dim1;
#line 281 "AB05RD.f"
    g -= g_offset;
#line 281 "AB05RD.f"
    h_dim1 = *ldh;
#line 281 "AB05RD.f"
    h_offset = 1 + h_dim1;
#line 281 "AB05RD.f"
    h__ -= h_offset;
#line 281 "AB05RD.f"
    bc_dim1 = *ldbc;
#line 281 "AB05RD.f"
    bc_offset = 1 + bc_dim1;
#line 281 "AB05RD.f"
    bc -= bc_offset;
#line 281 "AB05RD.f"
    cc_dim1 = *ldcc;
#line 281 "AB05RD.f"
    cc_offset = 1 + cc_dim1;
#line 281 "AB05RD.f"
    cc -= cc_offset;
#line 281 "AB05RD.f"
    dc_dim1 = *lddc;
#line 281 "AB05RD.f"
    dc_offset = 1 + dc_dim1;
#line 281 "AB05RD.f"
    dc -= dc_offset;
#line 281 "AB05RD.f"
    --iwork;
#line 281 "AB05RD.f"
    --dwork;
#line 281 "AB05RD.f"

#line 281 "AB05RD.f"
    /* Function Body */
#line 281 "AB05RD.f"
    unitf = lsame_(fbtype, "I", (ftnlen)1, (ftnlen)1);
#line 282 "AB05RD.f"
    outpf = lsame_(fbtype, "O", (ftnlen)1, (ftnlen)1);
#line 283 "AB05RD.f"
    ljobd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);

#line 285 "AB05RD.f"
    *info = 0;

#line 287 "AB05RD.f"
    if (! unitf && ! outpf) {
#line 288 "AB05RD.f"
	*info = -1;
#line 289 "AB05RD.f"
    } else if (! ljobd && ! lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1)) {
#line 290 "AB05RD.f"
	*info = -2;
#line 291 "AB05RD.f"
    } else if (*n < 0) {
#line 292 "AB05RD.f"
	*info = -3;
#line 293 "AB05RD.f"
    } else if (*m < 0) {
#line 294 "AB05RD.f"
	*info = -4;
#line 295 "AB05RD.f"
    } else if (*p < 0 || unitf && *p != *m) {
#line 296 "AB05RD.f"
	*info = -5;
#line 297 "AB05RD.f"
    } else if (*mv < 0) {
#line 298 "AB05RD.f"
	*info = -6;
#line 299 "AB05RD.f"
    } else if (*pz < 0) {
#line 300 "AB05RD.f"
	*info = -7;
#line 301 "AB05RD.f"
    } else if (*lda < max(1,*n)) {
#line 302 "AB05RD.f"
	*info = -11;
#line 303 "AB05RD.f"
    } else if (*ldb < max(1,*n)) {
#line 304 "AB05RD.f"
	*info = -13;
#line 305 "AB05RD.f"
    } else if (*n > 0 && *ldc < max(1,*p) || *n == 0 && *ldc < 1) {
#line 307 "AB05RD.f"
	*info = -15;
#line 308 "AB05RD.f"
    } else if (ljobd && *ldd < max(1,*p) || ! ljobd && *ldd < 1) {
#line 310 "AB05RD.f"
	*info = -17;
#line 311 "AB05RD.f"
    } else if (outpf && *alpha != 0. && *ldf < max(1,*m) || (unitf || *alpha 
	    == 0.) && *ldf < 1) {
#line 313 "AB05RD.f"
	*info = -19;
#line 314 "AB05RD.f"
    } else if (*beta != 0. && *ldk < max(1,*m) || *beta == 0. && *ldk < 1) {
#line 316 "AB05RD.f"
	*info = -21;
#line 317 "AB05RD.f"
    } else if (*ldg < max(1,*m)) {
#line 318 "AB05RD.f"
	*info = -23;
#line 319 "AB05RD.f"
    } else if (*ldh < max(1,*pz)) {
#line 320 "AB05RD.f"
	*info = -25;
#line 321 "AB05RD.f"
    } else if (*ldbc < max(1,*n)) {
#line 322 "AB05RD.f"
	*info = -28;
#line 323 "AB05RD.f"
    } else if (*n > 0 && *ldcc < max(1,*pz) || *n == 0 && *ldcc < 1) {
#line 325 "AB05RD.f"
	*info = -30;
#line 326 "AB05RD.f"
    } else if (ljobd && *lddc < max(1,*pz) || ! ljobd && *lddc < 1) {
#line 328 "AB05RD.f"
	*info = -32;
#line 329 "AB05RD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 329 "AB05RD.f"
	i__1 = max(1,*m), i__2 = *p * *mv, i__1 = max(i__1,i__2), i__2 = *p * 
		*p + (*p << 2);
#line 329 "AB05RD.f"
	if (ljobd && *ldwork < max(i__1,i__2) || ! ljobd && *ldwork < max(1,*
		m)) {
#line 331 "AB05RD.f"
	    *info = -35;
#line 332 "AB05RD.f"
	}
#line 332 "AB05RD.f"
    }

#line 334 "AB05RD.f"
    if (*info != 0) {

/*        Error return. */

#line 338 "AB05RD.f"
	i__1 = -(*info);
#line 338 "AB05RD.f"
	xerbla_("AB05RD", &i__1, (ftnlen)6);
#line 339 "AB05RD.f"
	return 0;
#line 340 "AB05RD.f"
    }

/*     Quick return if possible. */

/* Computing MAX */
#line 344 "AB05RD.f"
    i__1 = *n, i__2 = min(*m,*p), i__1 = max(i__1,i__2), i__2 = min(*mv,*pz);
#line 344 "AB05RD.f"
    if (max(i__1,i__2) == 0) {
#line 345 "AB05RD.f"
	*rcond = 1.;
#line 346 "AB05RD.f"
	return 0;
#line 347 "AB05RD.f"
    }

/*     Apply the partial output feedback u = alpha*F*y + v1 */

#line 351 "AB05RD.f"
    ab05sd_(fbtype, jobd, n, m, p, alpha, &a[a_offset], lda, &b[b_offset], 
	    ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &f[f_offset], ldf, 
	    rcond, &iwork[1], &dwork[1], ldwork, info, (ftnlen)1, (ftnlen)1);
#line 354 "AB05RD.f"
    if (*info != 0) {
#line 354 "AB05RD.f"
	return 0;
#line 354 "AB05RD.f"
    }

/*     Apply the partial state feedback v1 = beta*K*x + v2. */

/*     Compute Ac = A1 + beta*B1*K and C1 <- C1 + beta*D1*K. */

#line 360 "AB05RD.f"
    if (*beta != 0. && *n > 0) {
#line 361 "AB05RD.f"
	dgemm_("N", "N", n, n, m, beta, &b[b_offset], ldb, &k[k_offset], ldk, 
		&c_b9, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
#line 363 "AB05RD.f"
	if (ljobd) {
#line 363 "AB05RD.f"
	    dgemm_("N", "N", p, n, m, beta, &d__[d_offset], ldd, &k[k_offset],
		     ldk, &c_b9, &c__[c_offset], ldc, (ftnlen)1, (ftnlen)1);
#line 363 "AB05RD.f"
	}
#line 366 "AB05RD.f"
    }

/*     Apply the input and output conversions v2 = G*v, z = H*y. */

/*     Compute Bc = B1*G. */

#line 372 "AB05RD.f"
    dgemm_("N", "N", n, mv, m, &c_b9, &b[b_offset], ldb, &g[g_offset], ldg, &
	    c_b16, &bc[bc_offset], ldbc, (ftnlen)1, (ftnlen)1);

/*     Compute Cc = H*C1. */

#line 377 "AB05RD.f"
    if (*n > 0) {
#line 377 "AB05RD.f"
	dgemm_("N", "N", pz, n, p, &c_b9, &h__[h_offset], ldh, &c__[c_offset],
		 ldc, &c_b16, &cc[cc_offset], ldcc, (ftnlen)1, (ftnlen)1);
#line 377 "AB05RD.f"
    }

/*     Compute Dc = H*D1*G. */

#line 383 "AB05RD.f"
    if (ljobd) {
#line 384 "AB05RD.f"
	ldwp = max(1,*p);
#line 385 "AB05RD.f"
	dgemm_("N", "N", p, mv, m, &c_b9, &d__[d_offset], ldd, &g[g_offset], 
		ldg, &c_b16, &dwork[1], &ldwp, (ftnlen)1, (ftnlen)1);
#line 387 "AB05RD.f"
	dgemm_("N", "N", pz, mv, p, &c_b9, &h__[h_offset], ldh, &dwork[1], &
		ldwp, &c_b16, &dc[dc_offset], lddc, (ftnlen)1, (ftnlen)1);
#line 389 "AB05RD.f"
    }

#line 391 "AB05RD.f"
    return 0;
/* *** Last line of AB05RD *** */
} /* ab05rd_ */

