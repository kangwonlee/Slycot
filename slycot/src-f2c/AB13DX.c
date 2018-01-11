#line 1 "AB13DX.f"
/* AB13DX.f -- translated by f2c (version 20100827).
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

#line 1 "AB13DX.f"
/* Table of constant values */

static doublecomplex c_b1 = {1.,0.};
static doublereal c_b17 = -1.;
static doublereal c_b18 = 1.;
static doublereal c_b24 = 0.;

doublereal ab13dx_(char *dico, char *jobe, char *jobd, integer *n, integer *m,
	 integer *p, doublereal *omega, doublereal *a, integer *lda, 
	doublereal *e, integer *lde, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, integer *iwork, 
	doublereal *dwork, integer *ldwork, doublecomplex *cwork, integer *
	lcwork, integer *info, ftnlen dico_len, ftnlen jobe_len, ftnlen 
	jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, i__1, i__2, i__3;
    doublereal ret_val, d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static integer i__, j, id, is, icb, icc, icd;
    static doublereal upd;
    static integer icwk, ierr, iwrk;
    extern /* Subroutine */ int mb02rd_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    integer *, ftnlen), mb02sd_(integer *, doublereal *, integer *, 
	    integer *, integer *), dgemm_(char *, char *, integer *, integer *
	    , integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr, specl, fulle;
    extern /* Subroutine */ int mb02rz_(char *, integer *, integer *, 
	    doublecomplex *, integer *, integer *, doublecomplex *, integer *,
	     integer *, ftnlen);
    static doublereal bnorm, cnorm;
    static logical withd;
    static integer minpm;
    extern /* Subroutine */ int mb02sz_(integer *, doublecomplex *, integer *,
	     integer *, integer *), zgemm_(char *, char *, integer *, integer 
	    *, integer *, doublecomplex *, doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, doublecomplex *, 
	    integer *, ftnlen, ftnlen);
    static logical nodyn;
    extern /* Subroutine */ int zlacp2_(char *, integer *, integer *, 
	    doublereal *, integer *, doublecomplex *, integer *, ftnlen);
    static doublereal lambdi;
    extern doublereal dlange_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    static doublereal lambdr;
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), xerbla_(char *, integer *, ftnlen);
    static integer mincwr;
    extern /* Subroutine */ int zgesvd_(char *, char *, integer *, integer *, 
	    doublecomplex *, integer *, doublereal *, doublecomplex *, 
	    integer *, doublecomplex *, integer *, doublecomplex *, integer *,
	     doublereal *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;


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

/*     To compute the maximum singular value of a given continuous-time */
/*     or discrete-time transfer-function matrix, either standard or in */
/*     the descriptor form, */

/*                                     -1 */
/*        G(lambda) = C*( lambda*E - A ) *B + D , */

/*     for a given complex value lambda, where lambda = j*omega, in the */
/*     continuous-time case, and lambda = exp(j*omega), in the */
/*     discrete-time case. The matrices A, E, B, C, and D are real */
/*     matrices of appropriate dimensions. Matrix A must be in an upper */
/*     Hessenberg form, and if JOBE ='G', the matrix E must be upper */
/*     triangular. The matrices B and C must correspond to the system */
/*     in (generalized) Hessenberg form. */

/*     FUNCTION VALUE */

/*     AB13DX   DOUBLE PRECISION */
/*              The maximum singular value of G(lambda). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system, as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOBE    CHARACTER*1 */
/*             Specifies whether E is an upper triangular or an identity */
/*             matrix, as follows: */
/*             = 'G':  E is a general upper triangular matrix; */
/*             = 'I':  E is the identity matrix. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The row size of the matrix C.  P >= 0. */

/*     OMEGA   (input) DOUBLE PRECISION */
/*             The frequency value for which the calculations should be */
/*             done. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N upper Hessenberg part of this */
/*             array must contain the state dynamics matrix A in upper */
/*             Hessenberg form. The elements below the subdiagonal are */
/*             not referenced. */
/*             On exit, if M > 0, P > 0, OMEGA = 0, DICO = 'C', B <> 0, */
/*             and C <> 0, the leading N-by-N upper Hessenberg part of */
/*             this array contains the factors L and U from the LU */
/*             factorization of A (A = P*L*U); the unit diagonal elements */
/*             of L are not stored, L is lower bidiagonal, and P is */
/*             stored in IWORK (see SLICOT Library routine MB02SD). */
/*             Otherwise, this array is unchanged on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             If JOBE = 'G', the leading N-by-N upper triangular part of */
/*             this array must contain the upper triangular descriptor */
/*             matrix E of the system. The elements of the strict lower */
/*             triangular part of this array are not referenced. */
/*             If JOBE = 'I', then E is assumed to be the identity */
/*             matrix and is not referenced. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E. */
/*             LDE >= MAX(1,N), if JOBE = 'G'; */
/*             LDE >= 1,        if JOBE = 'I'. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the system input matrix B. */
/*             On exit, if M > 0, P > 0, OMEGA = 0, DICO = 'C', B <> 0, */
/*             C <> 0, and INFO = 0 or N+1, the leading N-by-M part of */
/*             this array contains the solution of the system A*X = B. */
/*             Otherwise, this array is unchanged on exit. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             system output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, if JOBD = 'D', the leading P-by-M part of this */
/*             array must contain the direct transmission matrix D. */
/*             On exit, if (N = 0, or B = 0, or C = 0) and JOBD = 'D', */
/*             or (OMEGA = 0, DICO = 'C', JOBD = 'D', and INFO = 0 or */
/*             N+1), the contents of this array is destroyed. */
/*             Otherwise, this array is unchanged on exit. */
/*             This array is not referenced if JOBD = 'Z'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P), if JOBD = 'D'; */
/*             LDD >= 1,        if JOBD = 'Z'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), where */
/*             LIWORK = N, if N > 0, M > 0, P > 0, B <> 0, and C <> 0; */
/*             LIWORK = 0, otherwise. */
/*             This array contains the pivot indices in the LU */
/*             factorization of the matrix lambda*E - A; for 1 <= i <= N, */
/*             row i of the matrix was interchanged with row IWORK(i). */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK, and DWORK(2), ..., DWORK(MIN(P,M)) contain the */
/*             singular values of G(lambda), except for the first one, */
/*             which is returned in the function value AB13DX. */
/*             If (N = 0, or B = 0, or C = 0) and JOBD = 'Z', the last */
/*             MIN(P,M)-1 zero singular values of G(lambda) are not */
/*             stored in DWORK(2), ..., DWORK(MIN(P,M)). */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= MAX(1, LDW1 + LDW2 ), */
/*             LDW1 = P*M, if N > 0, B <> 0, C <> 0, OMEGA = 0, */
/*                            DICO = 'C', and JOBD = 'Z'; */
/*             LDW1 = 0,   otherwise; */
/*             LDW2 = MIN(P,M) + MAX(3*MIN(P,M) + MAX(P,M), 5*MIN(P,M)), */
/*                         if (N = 0, or B = 0, or C = 0) and JOBD = 'D', */
/*                         or (N > 0, B <> 0, C <> 0, OMEGA = 0, and */
/*                             DICO = 'C'); */
/*             LDW2 = 0,   if (N = 0, or B = 0, or C = 0) and JOBD = 'Z', */
/*                         or MIN(P,M) = 0; */
/*             LDW2 = 6*MIN(P,M), otherwise. */
/*             For good performance, LDWORK must generally be larger. */

/*     CWORK   COMPLEX*16 array, dimension (LCWORK) */
/*             On exit, if INFO = 0, CWORK(1) contains the optimal */
/*             LCWORK. */

/*     LCWORK  INTEGER */
/*             The dimension of the array CWORK. */
/*             LCWORK >= 1, if N = 0, or B = 0, or C = 0, or (OMEGA = 0 */
/*                             and DICO = 'C') or MIN(P,M) = 0; */
/*             LCWORK >= MAX(1, (N+M)*(N+P) + 2*MIN(P,M) + MAX(P,M)), */
/*                          otherwise. */
/*             For good performance, LCWORK must generally be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, U(i,i) is exactly zero; the LU */
/*                   factorization of the matrix lambda*E - A has been */
/*                   completed, but the factor U is exactly singular, */
/*                   i.e., the matrix lambda*E - A is exactly singular; */
/*             = N+1:  the SVD algorithm for computing singular values */
/*                   did not converge. */

/*     METHOD */

/*     The routine implements standard linear algebra calculations, */
/*     taking problem structure into account. LAPACK Library routines */
/*     DGESVD and ZGESVD are used for finding the singular values. */

/*     CONTRIBUTORS */

/*     D. Sima, University of Bucharest, May 2001. */
/*     V. Sima, Research Institute for Informatics, Bucharest, May 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Sep. 2005. */

/*     KEYWORDS */

/*     H-infinity optimal control, robust control, system norm. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */

/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar parameters. */

#line 262 "AB13DX.f"
    /* Parameter adjustments */
#line 262 "AB13DX.f"
    a_dim1 = *lda;
#line 262 "AB13DX.f"
    a_offset = 1 + a_dim1;
#line 262 "AB13DX.f"
    a -= a_offset;
#line 262 "AB13DX.f"
    e_dim1 = *lde;
#line 262 "AB13DX.f"
    e_offset = 1 + e_dim1;
#line 262 "AB13DX.f"
    e -= e_offset;
#line 262 "AB13DX.f"
    b_dim1 = *ldb;
#line 262 "AB13DX.f"
    b_offset = 1 + b_dim1;
#line 262 "AB13DX.f"
    b -= b_offset;
#line 262 "AB13DX.f"
    c_dim1 = *ldc;
#line 262 "AB13DX.f"
    c_offset = 1 + c_dim1;
#line 262 "AB13DX.f"
    c__ -= c_offset;
#line 262 "AB13DX.f"
    d_dim1 = *ldd;
#line 262 "AB13DX.f"
    d_offset = 1 + d_dim1;
#line 262 "AB13DX.f"
    d__ -= d_offset;
#line 262 "AB13DX.f"
    --iwork;
#line 262 "AB13DX.f"
    --dwork;
#line 262 "AB13DX.f"
    --cwork;
#line 262 "AB13DX.f"

#line 262 "AB13DX.f"
    /* Function Body */
#line 262 "AB13DX.f"
    *info = 0;
#line 263 "AB13DX.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 264 "AB13DX.f"
    fulle = lsame_(jobe, "G", (ftnlen)1, (ftnlen)1);
#line 265 "AB13DX.f"
    withd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);

#line 267 "AB13DX.f"
    if (! (discr || lsame_(dico, "C", (ftnlen)1, (ftnlen)1))) {
#line 268 "AB13DX.f"
	*info = -1;
#line 269 "AB13DX.f"
    } else if (! (fulle || lsame_(jobe, "I", (ftnlen)1, (ftnlen)1))) {
#line 270 "AB13DX.f"
	*info = -2;
#line 271 "AB13DX.f"
    } else if (! (withd || lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1))) {
#line 272 "AB13DX.f"
	*info = -3;
#line 273 "AB13DX.f"
    } else if (*n < 0) {
#line 274 "AB13DX.f"
	*info = -4;
#line 275 "AB13DX.f"
    } else if (*m < 0) {
#line 276 "AB13DX.f"
	*info = -5;
#line 277 "AB13DX.f"
    } else if (*p < 0) {
#line 278 "AB13DX.f"
	*info = -6;
#line 279 "AB13DX.f"
    } else if (*lda < max(1,*n)) {
#line 280 "AB13DX.f"
	*info = -9;
#line 281 "AB13DX.f"
    } else if (*lde < 1 || fulle && *lde < *n) {
#line 282 "AB13DX.f"
	*info = -11;
#line 283 "AB13DX.f"
    } else if (*ldb < max(1,*n)) {
#line 284 "AB13DX.f"
	*info = -13;
#line 285 "AB13DX.f"
    } else if (*ldc < max(1,*p)) {
#line 286 "AB13DX.f"
	*info = -15;
#line 287 "AB13DX.f"
    } else if (*ldd < 1 || withd && *ldd < *p) {
#line 288 "AB13DX.f"
	*info = -17;
#line 289 "AB13DX.f"
    } else {
#line 290 "AB13DX.f"
	bnorm = dlange_("1-norm", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)
		6);
#line 291 "AB13DX.f"
	cnorm = dlange_("1-norm", p, n, &c__[c_offset], ldc, &dwork[1], (
		ftnlen)6);
#line 292 "AB13DX.f"
	nodyn = *n == 0 || min(bnorm,cnorm) == 0.;
#line 293 "AB13DX.f"
	specl = ! nodyn && *omega == 0. && ! discr;
#line 294 "AB13DX.f"
	minpm = min(*p,*m);

/*        Compute workspace. */

#line 298 "AB13DX.f"
	if (minpm == 0) {
#line 299 "AB13DX.f"
	    minwrk = 0;
#line 300 "AB13DX.f"
	} else if (specl || nodyn && withd) {
/* Computing MAX */
#line 301 "AB13DX.f"
	    i__1 = minpm * 3 + max(*p,*m), i__2 = minpm * 5;
#line 301 "AB13DX.f"
	    minwrk = minpm + max(i__1,i__2);
#line 302 "AB13DX.f"
	    if (specl && ! withd) {
#line 302 "AB13DX.f"
		minwrk += *p * *m;
#line 302 "AB13DX.f"
	    }
#line 304 "AB13DX.f"
	} else if (nodyn && ! withd) {
#line 305 "AB13DX.f"
	    minwrk = 0;
#line 306 "AB13DX.f"
	} else {
#line 307 "AB13DX.f"
	    minwrk = minpm * 6;
#line 308 "AB13DX.f"
	}
#line 309 "AB13DX.f"
	minwrk = max(1,minwrk);

#line 311 "AB13DX.f"
	if (*ldwork < minwrk) {
#line 312 "AB13DX.f"
	    *info = -20;
#line 313 "AB13DX.f"
	} else {
#line 314 "AB13DX.f"
	    if (nodyn || *omega == 0. && ! discr || minpm == 0) {
#line 316 "AB13DX.f"
		mincwr = 1;
#line 317 "AB13DX.f"
	    } else {
/* Computing MAX */
#line 318 "AB13DX.f"
		i__1 = 1, i__2 = (*n + *m) * (*n + *p) + (minpm << 1) + max(*
			p,*m);
#line 318 "AB13DX.f"
		mincwr = max(i__1,i__2);
#line 320 "AB13DX.f"
	    }
#line 321 "AB13DX.f"
	    if (*lcwork < mincwr) {
#line 321 "AB13DX.f"
		*info = -22;
#line 321 "AB13DX.f"
	    }
#line 323 "AB13DX.f"
	}
#line 324 "AB13DX.f"
    }

#line 326 "AB13DX.f"
    if (*info != 0) {
#line 327 "AB13DX.f"
	i__1 = -(*info);
#line 327 "AB13DX.f"
	xerbla_("AB13DX", &i__1, (ftnlen)6);
#line 328 "AB13DX.f"
	return ret_val;
#line 329 "AB13DX.f"
    }

/*     Quick return if possible. */

#line 333 "AB13DX.f"
    if (minpm == 0) {
#line 334 "AB13DX.f"
	ret_val = 0.;

#line 336 "AB13DX.f"
	dwork[1] = 1.;
#line 337 "AB13DX.f"
	cwork[1].r = 1., cwork[1].i = 0.;
#line 338 "AB13DX.f"
	return ret_val;
#line 339 "AB13DX.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance.) */

#line 345 "AB13DX.f"
    is = 1;
#line 346 "AB13DX.f"
    iwrk = is + minpm;

#line 348 "AB13DX.f"
    if (nodyn) {

/*        No dynamics: Determine the maximum singular value of G = D . */

#line 352 "AB13DX.f"
	if (withd) {

/*           Workspace: need   MIN(P,M) + MAX(3*MIN(P,M) + MAX(P,M), */
/*                                            5*MIN(P,M)); */
/*                      prefer larger. */

#line 358 "AB13DX.f"
	    i__1 = *ldwork - iwrk + 1;
#line 358 "AB13DX.f"
	    dgesvd_("No Vectors", "No Vectors", p, m, &d__[d_offset], ldd, &
		    dwork[is], &dwork[1], p, &dwork[1], m, &dwork[iwrk], &
		    i__1, &ierr, (ftnlen)10, (ftnlen)10);
#line 361 "AB13DX.f"
	    if (ierr > 0) {
#line 362 "AB13DX.f"
		*info = *n + 1;
#line 363 "AB13DX.f"
		return ret_val;
#line 364 "AB13DX.f"
	    }
#line 365 "AB13DX.f"
	    ret_val = dwork[is];
#line 366 "AB13DX.f"
	    maxwrk = (integer) dwork[iwrk] + iwrk - 1;
#line 367 "AB13DX.f"
	} else {
#line 368 "AB13DX.f"
	    ret_val = 0.;
#line 369 "AB13DX.f"
	    maxwrk = 1;
#line 370 "AB13DX.f"
	}

#line 372 "AB13DX.f"
	dwork[1] = (doublereal) maxwrk;
#line 373 "AB13DX.f"
	cwork[1].r = 1., cwork[1].i = 0.;
#line 374 "AB13DX.f"
	return ret_val;
#line 375 "AB13DX.f"
    }

/*     Determine the maximum singular value of */
/*        G(lambda) = C*inv(lambda*E - A)*B + D. */
/*     The (generalized) Hessenberg form of the system is used. */

#line 381 "AB13DX.f"
    if (specl) {

/*        Special continuous-time case: */
/*        Determine the maximum singular value of the real matrix G(0). */
/*        Workspace: need   MIN(P,M) + MAX(3*MIN(P,M) + MAX(P,M), */
/*                                         5*MIN(P,M)); */
/*                   prefer larger. */

#line 389 "AB13DX.f"
	mb02sd_(n, &a[a_offset], lda, &iwork[1], &ierr);
#line 390 "AB13DX.f"
	if (ierr > 0) {
#line 391 "AB13DX.f"
	    *info = ierr;
#line 392 "AB13DX.f"
	    dwork[1] = 1.;
#line 393 "AB13DX.f"
	    cwork[1].r = 1., cwork[1].i = 0.;
#line 394 "AB13DX.f"
	    return ret_val;
#line 395 "AB13DX.f"
	}
#line 396 "AB13DX.f"
	mb02rd_("No Transpose", n, m, &a[a_offset], lda, &iwork[1], &b[
		b_offset], ldb, &ierr, (ftnlen)12);
#line 398 "AB13DX.f"
	if (withd) {
#line 399 "AB13DX.f"
	    dgemm_("No Transpose", "No Transpose", p, m, n, &c_b17, &c__[
		    c_offset], ldc, &b[b_offset], ldb, &c_b18, &d__[d_offset],
		     ldd, (ftnlen)12, (ftnlen)12);
#line 401 "AB13DX.f"
	    i__1 = *ldwork - iwrk + 1;
#line 401 "AB13DX.f"
	    dgesvd_("No Vectors", "No Vectors", p, m, &d__[d_offset], ldd, &
		    dwork[is], &dwork[1], p, &dwork[1], m, &dwork[iwrk], &
		    i__1, &ierr, (ftnlen)10, (ftnlen)10);
#line 404 "AB13DX.f"
	} else {

/*           Additional workspace: need   P*M. */

#line 408 "AB13DX.f"
	    id = iwrk;
#line 409 "AB13DX.f"
	    iwrk = id + *p * *m;
#line 410 "AB13DX.f"
	    dgemm_("No Transpose", "No Transpose", p, m, n, &c_b17, &c__[
		    c_offset], ldc, &b[b_offset], ldb, &c_b24, &dwork[id], p, 
		    (ftnlen)12, (ftnlen)12);
#line 412 "AB13DX.f"
	    i__1 = *ldwork - iwrk + 1;
#line 412 "AB13DX.f"
	    dgesvd_("No Vectors", "No Vectors", p, m, &dwork[id], p, &dwork[
		    is], &dwork[1], p, &dwork[1], m, &dwork[iwrk], &i__1, &
		    ierr, (ftnlen)10, (ftnlen)10);
#line 415 "AB13DX.f"
	}
#line 416 "AB13DX.f"
	if (ierr > 0) {
#line 417 "AB13DX.f"
	    *info = *n + 1;
#line 418 "AB13DX.f"
	    return ret_val;
#line 419 "AB13DX.f"
	}

#line 421 "AB13DX.f"
	ret_val = dwork[is];
#line 422 "AB13DX.f"
	dwork[1] = (doublereal) ((integer) dwork[iwrk] + iwrk - 1);
#line 423 "AB13DX.f"
	cwork[1].r = 1., cwork[1].i = 0.;
#line 424 "AB13DX.f"
	return ret_val;
#line 425 "AB13DX.f"
    }

/*     General case: Determine the maximum singular value of G(lambda). */
/*     Complex workspace:  need   N*N + N*M + P*N + P*M. */

#line 430 "AB13DX.f"
    icb = *n * *n + 1;
#line 431 "AB13DX.f"
    icc = icb + *n * *m;
#line 432 "AB13DX.f"
    icd = icc + *p * *n;
#line 433 "AB13DX.f"
    icwk = icd + *p * *m;

#line 435 "AB13DX.f"
    if (withd) {
#line 436 "AB13DX.f"
	upd = 1.;
#line 437 "AB13DX.f"
    } else {
#line 438 "AB13DX.f"
	upd = 0.;
#line 439 "AB13DX.f"
    }

#line 441 "AB13DX.f"
    if (discr) {
#line 442 "AB13DX.f"
	lambdr = cos(*omega);
#line 443 "AB13DX.f"
	lambdi = sin(*omega);

/*        Build lambda*E - A . */

#line 447 "AB13DX.f"
	if (fulle) {

#line 449 "AB13DX.f"
	    i__1 = *n;
#line 449 "AB13DX.f"
	    for (j = 1; j <= i__1; ++j) {

#line 451 "AB13DX.f"
		i__2 = j;
#line 451 "AB13DX.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 452 "AB13DX.f"
		    i__3 = i__ + (j - 1) * *n;
#line 452 "AB13DX.f"
		    d__1 = lambdr * e[i__ + j * e_dim1] - a[i__ + j * a_dim1];
#line 452 "AB13DX.f"
		    d__2 = lambdi * e[i__ + j * e_dim1];
#line 452 "AB13DX.f"
		    z__1.r = d__1, z__1.i = d__2;
#line 452 "AB13DX.f"
		    cwork[i__3].r = z__1.r, cwork[i__3].i = z__1.i;
#line 455 "AB13DX.f"
/* L10: */
#line 455 "AB13DX.f"
		}

#line 457 "AB13DX.f"
		if (j < *n) {
#line 457 "AB13DX.f"
		    i__2 = j + 1 + (j - 1) * *n;
#line 457 "AB13DX.f"
		    d__1 = -a[j + 1 + j * a_dim1];
#line 457 "AB13DX.f"
		    z__1.r = d__1, z__1.i = 0.;
#line 457 "AB13DX.f"
		    cwork[i__2].r = z__1.r, cwork[i__2].i = z__1.i;
#line 457 "AB13DX.f"
		}
#line 459 "AB13DX.f"
/* L20: */
#line 459 "AB13DX.f"
	    }

#line 461 "AB13DX.f"
	} else {

#line 463 "AB13DX.f"
	    i__1 = *n;
#line 463 "AB13DX.f"
	    for (j = 1; j <= i__1; ++j) {

/* Computing MIN */
#line 465 "AB13DX.f"
		i__3 = j + 1;
#line 465 "AB13DX.f"
		i__2 = min(i__3,*n);
#line 465 "AB13DX.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 466 "AB13DX.f"
		    i__3 = i__ + (j - 1) * *n;
#line 466 "AB13DX.f"
		    d__1 = -a[i__ + j * a_dim1];
#line 466 "AB13DX.f"
		    cwork[i__3].r = d__1, cwork[i__3].i = 0.;
#line 467 "AB13DX.f"
/* L30: */
#line 467 "AB13DX.f"
		}

#line 469 "AB13DX.f"
		i__2 = j + (j - 1) * *n;
#line 469 "AB13DX.f"
		d__1 = lambdr - a[j + j * a_dim1];
#line 469 "AB13DX.f"
		z__1.r = d__1, z__1.i = lambdi;
#line 469 "AB13DX.f"
		cwork[i__2].r = z__1.r, cwork[i__2].i = z__1.i;
#line 470 "AB13DX.f"
/* L40: */
#line 470 "AB13DX.f"
	    }

#line 472 "AB13DX.f"
	}

#line 474 "AB13DX.f"
    } else {

/*        Build j*omega*E - A. */

#line 478 "AB13DX.f"
	if (fulle) {

#line 480 "AB13DX.f"
	    i__1 = *n;
#line 480 "AB13DX.f"
	    for (j = 1; j <= i__1; ++j) {

#line 482 "AB13DX.f"
		i__2 = j;
#line 482 "AB13DX.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 483 "AB13DX.f"
		    i__3 = i__ + (j - 1) * *n;
#line 483 "AB13DX.f"
		    d__1 = -a[i__ + j * a_dim1];
#line 483 "AB13DX.f"
		    d__2 = *omega * e[i__ + j * e_dim1];
#line 483 "AB13DX.f"
		    z__1.r = d__1, z__1.i = d__2;
#line 483 "AB13DX.f"
		    cwork[i__3].r = z__1.r, cwork[i__3].i = z__1.i;
#line 485 "AB13DX.f"
/* L50: */
#line 485 "AB13DX.f"
		}

#line 487 "AB13DX.f"
		if (j < *n) {
#line 487 "AB13DX.f"
		    i__2 = j + 1 + (j - 1) * *n;
#line 487 "AB13DX.f"
		    d__1 = -a[j + 1 + j * a_dim1];
#line 487 "AB13DX.f"
		    z__1.r = d__1, z__1.i = 0.;
#line 487 "AB13DX.f"
		    cwork[i__2].r = z__1.r, cwork[i__2].i = z__1.i;
#line 487 "AB13DX.f"
		}
#line 489 "AB13DX.f"
/* L60: */
#line 489 "AB13DX.f"
	    }

#line 491 "AB13DX.f"
	} else {

#line 493 "AB13DX.f"
	    i__1 = *n;
#line 493 "AB13DX.f"
	    for (j = 1; j <= i__1; ++j) {

/* Computing MIN */
#line 495 "AB13DX.f"
		i__3 = j + 1;
#line 495 "AB13DX.f"
		i__2 = min(i__3,*n);
#line 495 "AB13DX.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 496 "AB13DX.f"
		    i__3 = i__ + (j - 1) * *n;
#line 496 "AB13DX.f"
		    d__1 = -a[i__ + j * a_dim1];
#line 496 "AB13DX.f"
		    cwork[i__3].r = d__1, cwork[i__3].i = 0.;
#line 497 "AB13DX.f"
/* L70: */
#line 497 "AB13DX.f"
		}

#line 499 "AB13DX.f"
		i__2 = j + (j - 1) * *n;
#line 499 "AB13DX.f"
		d__1 = -a[j + j * a_dim1];
#line 499 "AB13DX.f"
		z__1.r = d__1, z__1.i = *omega;
#line 499 "AB13DX.f"
		cwork[i__2].r = z__1.r, cwork[i__2].i = z__1.i;
#line 500 "AB13DX.f"
/* L80: */
#line 500 "AB13DX.f"
	    }

#line 502 "AB13DX.f"
	}

#line 504 "AB13DX.f"
    }

/*     Build G(lambda) . */

#line 508 "AB13DX.f"
    zlacp2_("Full", n, m, &b[b_offset], ldb, &cwork[icb], n, (ftnlen)4);
#line 509 "AB13DX.f"
    zlacp2_("Full", p, n, &c__[c_offset], ldc, &cwork[icc], p, (ftnlen)4);
#line 510 "AB13DX.f"
    if (withd) {
#line 510 "AB13DX.f"
	zlacp2_("Full", p, m, &d__[d_offset], ldd, &cwork[icd], p, (ftnlen)4);
#line 510 "AB13DX.f"
    }

#line 513 "AB13DX.f"
    mb02sz_(n, &cwork[1], n, &iwork[1], &ierr);
#line 514 "AB13DX.f"
    if (ierr > 0) {
#line 515 "AB13DX.f"
	*info = ierr;
#line 516 "AB13DX.f"
	dwork[1] = 1.;
#line 517 "AB13DX.f"
	i__1 = icwk - 1;
#line 517 "AB13DX.f"
	cwork[1].r = (doublereal) i__1, cwork[1].i = 0.;
#line 518 "AB13DX.f"
	return ret_val;
#line 519 "AB13DX.f"
    }
#line 520 "AB13DX.f"
    mb02rz_("No Transpose", n, m, &cwork[1], n, &iwork[1], &cwork[icb], n, &
	    ierr, (ftnlen)12);
#line 522 "AB13DX.f"
    z__1.r = upd, z__1.i = 0.;
#line 522 "AB13DX.f"
    zgemm_("No Transpose", "No Transpose", p, m, n, &c_b1, &cwork[icc], p, &
	    cwork[icb], n, &z__1, &cwork[icd], p, (ftnlen)12, (ftnlen)12);

/*     Additional workspace, complex: need   2*MIN(P,M) + MAX(P,M); */
/*                                    prefer larger; */
/*                           real:    need   5*MIN(P,M). */

#line 530 "AB13DX.f"
    i__1 = *lcwork - icwk + 1;
#line 530 "AB13DX.f"
    zgesvd_("No Vectors", "No Vectors", p, m, &cwork[icd], p, &dwork[is], &
	    cwork[1], p, &cwork[1], m, &cwork[icwk], &i__1, &dwork[iwrk], &
	    ierr, (ftnlen)10, (ftnlen)10);
#line 533 "AB13DX.f"
    if (ierr > 0) {
#line 534 "AB13DX.f"
	*info = *n + 1;
#line 535 "AB13DX.f"
	return ret_val;
#line 536 "AB13DX.f"
    }
#line 537 "AB13DX.f"
    ret_val = dwork[is];

#line 539 "AB13DX.f"
    dwork[1] = (doublereal) (minpm * 6);
#line 540 "AB13DX.f"
    i__2 = icwk;
#line 540 "AB13DX.f"
    i__1 = (integer) cwork[i__2].r + icwk - 1;
#line 540 "AB13DX.f"
    cwork[1].r = (doublereal) i__1, cwork[1].i = 0.;

#line 542 "AB13DX.f"
    return ret_val;
/* *** Last line of AB13DX *** */
} /* ab13dx_ */

