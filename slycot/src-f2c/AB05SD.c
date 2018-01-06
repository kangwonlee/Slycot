#line 1 "AB05SD.f"
/* AB05SD.f -- translated by f2c (version 20100827).
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

#line 1 "AB05SD.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b11 = 1.;
static doublereal c_b14 = 0.;
static integer c__1 = 1;

/* Subroutine */ int ab05sd_(char *fbtype, char *jobd, integer *n, integer *m,
	 integer *p, doublereal *alpha, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *f, integer *ldf, 
	doublereal *rcond, integer *iwork, doublereal *dwork, integer *ldwork,
	 integer *info, ftnlen fbtype_len, ftnlen jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, f_dim1, f_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, iw, ldwn, ldwp;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static logical ljobd;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal enorm;
    static logical unitf;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static doublereal dummy[1];
    static logical outpf;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgecon_(char *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dgetrf_(integer *, integer *, 
	    doublereal *, integer *, integer *, integer *), dlacpy_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, ftnlen), xerbla_(char *, integer *, ftnlen), dgetrs_(
	    char *, integer *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, ftnlen);


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
/*     loop system (Ac,Bc,Cc,Dc) corresponding to the output feedback */
/*     control law */

/*          u = alpha*F*y + v. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     FBTYPE  CHARACTER*1 */
/*             Specifies the type of the feedback law as follows: */
/*             = 'I':  Unitary output feedback (F = I); */
/*             = 'O':  General output feedback. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of state variables, i.e. the order of the */
/*             matrix A, the number of rows of B and the number of */
/*             columns of C.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of input variables, i.e. the number of columns */
/*             of matrices B and D, and the number of rows of F.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of output variables, i.e. the number of rows of */
/*             matrices C and D, and the number of columns of F.  P >= 0 */
/*             and P = M if FBTYPE = 'I'. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The coefficient alpha in the output feedback law. */

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
/*             the input matrix Bc of the closed-loop system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the system output matrix C. */
/*             On exit, the leading P-by-N part of this array contains */
/*             the output matrix Cc of the closed-loop system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,P) if N > 0. */
/*             LDC >= 1 if N = 0. */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the system direct input/output transmission */
/*             matrix D. */
/*             On exit, if JOBD = 'D', the leading P-by-M part of this */
/*             array contains the direct input/output transmission */
/*             matrix Dc of the closed-loop system. */
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
/*             The array F is not referenced if FBTYPE = 'I' or */
/*             ALPHA = 0. */

/*     LDF     INTEGER */
/*             The leading dimension of array F. */
/*             LDF >= MAX(1,M) if FBTYPE = 'O' and ALPHA <> 0. */
/*             LDF >= 1 if FBTYPE = 'I' or ALPHA = 0. */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The reciprocal condition number of the matrix */
/*             I - alpha*D*F. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= MAX(1,2*P) if JOBD = 'D'. */
/*             LIWORK >= 1 if JOBD = 'Z'. */
/*             IWORK is not referenced if JOBD = 'Z'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= wspace, where */
/*                       wspace = MAX( 1, M, P*P + 4*P ) if JOBD = 'D', */
/*                       wspace = MAX( 1, M ) if JOBD = 'Z'. */
/*             For best performance, LDWORK >= MAX( wspace, N*M, N*P ). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the matrix I - alpha*D*F is numerically singular. */

/*     METHOD */

/*     The matrices of the closed-loop system have the expressions: */

/*     Ac = A + alpha*B*F*E*C,  Bc = B + alpha*B*F*E*D, */
/*     Cc = E*C,                Dc = E*D, */

/*     where E = (I - alpha*D*F)**-1. */

/*     NUMERICAL ASPECTS */

/*     The accuracy of computations basically depends on the conditioning */
/*     of the matrix I - alpha*D*F.  If RCOND is very small, it is likely */
/*     that the computed results are inaccurate. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Research Establishment, */
/*     Oberpfaffenhofen, Germany, and V. Sima, Katholieke Univ. Leuven, */
/*     Belgium, Nov. 1996. */

/*     REVISIONS */

/*     January 14, 1997. */
/*     V. Sima, Research Institute for Informatics, Bucharest, July 2003. */

/*     KEYWORDS */

/*     Multivariable system, state-space model, state-space */
/*     representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External functions .. */
/*     .. External subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the input scalar arguments. */

#line 213 "AB05SD.f"
    /* Parameter adjustments */
#line 213 "AB05SD.f"
    a_dim1 = *lda;
#line 213 "AB05SD.f"
    a_offset = 1 + a_dim1;
#line 213 "AB05SD.f"
    a -= a_offset;
#line 213 "AB05SD.f"
    b_dim1 = *ldb;
#line 213 "AB05SD.f"
    b_offset = 1 + b_dim1;
#line 213 "AB05SD.f"
    b -= b_offset;
#line 213 "AB05SD.f"
    c_dim1 = *ldc;
#line 213 "AB05SD.f"
    c_offset = 1 + c_dim1;
#line 213 "AB05SD.f"
    c__ -= c_offset;
#line 213 "AB05SD.f"
    d_dim1 = *ldd;
#line 213 "AB05SD.f"
    d_offset = 1 + d_dim1;
#line 213 "AB05SD.f"
    d__ -= d_offset;
#line 213 "AB05SD.f"
    f_dim1 = *ldf;
#line 213 "AB05SD.f"
    f_offset = 1 + f_dim1;
#line 213 "AB05SD.f"
    f -= f_offset;
#line 213 "AB05SD.f"
    --iwork;
#line 213 "AB05SD.f"
    --dwork;
#line 213 "AB05SD.f"

#line 213 "AB05SD.f"
    /* Function Body */
#line 213 "AB05SD.f"
    unitf = lsame_(fbtype, "I", (ftnlen)1, (ftnlen)1);
#line 214 "AB05SD.f"
    outpf = lsame_(fbtype, "O", (ftnlen)1, (ftnlen)1);
#line 215 "AB05SD.f"
    ljobd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);
#line 216 "AB05SD.f"
    ldwn = max(1,*n);
#line 217 "AB05SD.f"
    ldwp = max(1,*p);

#line 219 "AB05SD.f"
    *info = 0;

#line 221 "AB05SD.f"
    if (! unitf && ! outpf) {
#line 222 "AB05SD.f"
	*info = -1;
#line 223 "AB05SD.f"
    } else if (! ljobd && ! lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1)) {
#line 224 "AB05SD.f"
	*info = -2;
#line 225 "AB05SD.f"
    } else if (*n < 0) {
#line 226 "AB05SD.f"
	*info = -3;
#line 227 "AB05SD.f"
    } else if (*m < 0) {
#line 228 "AB05SD.f"
	*info = -4;
#line 229 "AB05SD.f"
    } else if (*p < 0 || unitf && *p != *m) {
#line 230 "AB05SD.f"
	*info = -5;
#line 231 "AB05SD.f"
    } else if (*lda < ldwn) {
#line 232 "AB05SD.f"
	*info = -7;
#line 233 "AB05SD.f"
    } else if (*ldb < ldwn) {
#line 234 "AB05SD.f"
	*info = -9;
#line 235 "AB05SD.f"
    } else if (*n > 0 && *ldc < ldwp || *n == 0 && *ldc < 1) {
#line 237 "AB05SD.f"
	*info = -11;
#line 238 "AB05SD.f"
    } else if (ljobd && *ldd < ldwp || ! ljobd && *ldd < 1) {
#line 240 "AB05SD.f"
	*info = -13;
#line 241 "AB05SD.f"
    } else if (outpf && *alpha != 0. && *ldf < max(1,*m) || (unitf || *alpha 
	    == 0.) && *ldf < 1) {
#line 243 "AB05SD.f"
	*info = -16;
#line 244 "AB05SD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 244 "AB05SD.f"
	i__1 = max(1,*m), i__2 = *p * *p + (*p << 2);
#line 244 "AB05SD.f"
	if (ljobd && *ldwork < max(i__1,i__2) || ! ljobd && *ldwork < max(1,*
		m)) {
#line 246 "AB05SD.f"
	    *info = -20;
#line 247 "AB05SD.f"
	}
#line 247 "AB05SD.f"
    }

#line 249 "AB05SD.f"
    if (*info != 0) {

/*        Error return. */

#line 253 "AB05SD.f"
	i__1 = -(*info);
#line 253 "AB05SD.f"
	xerbla_("AB05SD", &i__1, (ftnlen)6);
#line 254 "AB05SD.f"
	return 0;
#line 255 "AB05SD.f"
    }

/*     Quick return if possible. */

#line 259 "AB05SD.f"
    *rcond = 1.;
/* Computing MAX */
#line 260 "AB05SD.f"
    i__1 = *n, i__2 = min(*m,*p);
#line 260 "AB05SD.f"
    if (max(i__1,i__2) == 0 || *alpha == 0.) {
#line 260 "AB05SD.f"
	return 0;
#line 260 "AB05SD.f"
    }

#line 263 "AB05SD.f"
    if (ljobd) {
#line 264 "AB05SD.f"
	iw = *p * *p + 1;

/*        Compute I - alpha*D*F. */

#line 268 "AB05SD.f"
	if (unitf) {
#line 269 "AB05SD.f"
	    dlacpy_("F", p, p, &d__[d_offset], ldd, &dwork[1], &ldwp, (ftnlen)
		    1);
#line 270 "AB05SD.f"
	    if (*alpha != -1.) {
#line 270 "AB05SD.f"
		d__1 = -(*alpha);
#line 270 "AB05SD.f"
		dlascl_("G", &c__0, &c__0, &c_b11, &d__1, p, p, &dwork[1], &
			ldwp, info, (ftnlen)1);
#line 270 "AB05SD.f"
	    }
#line 273 "AB05SD.f"
	} else {
#line 274 "AB05SD.f"
	    d__1 = -(*alpha);
#line 274 "AB05SD.f"
	    dgemm_("N", "N", p, p, m, &d__1, &d__[d_offset], ldd, &f[f_offset]
		    , ldf, &c_b14, &dwork[1], &ldwp, (ftnlen)1, (ftnlen)1);
#line 276 "AB05SD.f"
	}

#line 278 "AB05SD.f"
	dummy[0] = 1.;
#line 279 "AB05SD.f"
	i__1 = *p + 1;
#line 279 "AB05SD.f"
	daxpy_(p, &c_b11, dummy, &c__0, &dwork[1], &i__1);

/*        Compute Cc = E*C, Dc = E*D, where E = (I - alpha*D*F)**-1. */

#line 283 "AB05SD.f"
	enorm = dlange_("1", p, p, &dwork[1], &ldwp, &dwork[iw], (ftnlen)1);
#line 284 "AB05SD.f"
	dgetrf_(p, p, &dwork[1], &ldwp, &iwork[1], info);
#line 285 "AB05SD.f"
	if (*info > 0) {

/*           Error return. */

#line 289 "AB05SD.f"
	    *rcond = 0.;
#line 290 "AB05SD.f"
	    *info = 1;
#line 291 "AB05SD.f"
	    return 0;
#line 292 "AB05SD.f"
	}
#line 293 "AB05SD.f"
	dgecon_("1", p, &dwork[1], &ldwp, &enorm, rcond, &dwork[iw], &iwork[*
		p + 1], info, (ftnlen)1);
#line 295 "AB05SD.f"
	if (*rcond <= dlamch_("E", (ftnlen)1)) {

/*           Error return. */

#line 299 "AB05SD.f"
	    *info = 1;
#line 300 "AB05SD.f"
	    return 0;
#line 301 "AB05SD.f"
	}

#line 303 "AB05SD.f"
	if (*n > 0) {
#line 303 "AB05SD.f"
	    dgetrs_("N", p, n, &dwork[1], &ldwp, &iwork[1], &c__[c_offset], 
		    ldc, info, (ftnlen)1);
#line 303 "AB05SD.f"
	}
#line 305 "AB05SD.f"
	dgetrs_("N", p, m, &dwork[1], &ldwp, &iwork[1], &d__[d_offset], ldd, 
		info, (ftnlen)1);
#line 306 "AB05SD.f"
    }

#line 308 "AB05SD.f"
    if (*n == 0) {
#line 308 "AB05SD.f"
	return 0;
#line 308 "AB05SD.f"
    }

/*     Compute Ac = A + alpha*B*F*Cc and Bc = B + alpha*B*F*Dc. */

#line 313 "AB05SD.f"
    if (unitf) {
#line 314 "AB05SD.f"
	dgemm_("N", "N", n, n, m, alpha, &b[b_offset], ldb, &c__[c_offset], 
		ldc, &c_b11, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
#line 316 "AB05SD.f"
	if (ljobd) {

#line 318 "AB05SD.f"
	    if (*ldwork < *n * *m) {

/*              Not enough working space for using DGEMM. */

#line 322 "AB05SD.f"
		i__1 = *n;
#line 322 "AB05SD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 323 "AB05SD.f"
		    dcopy_(p, &b[i__ + b_dim1], ldb, &dwork[1], &c__1);
#line 324 "AB05SD.f"
		    dgemv_("T", p, p, alpha, &d__[d_offset], ldd, &dwork[1], &
			    c__1, &c_b11, &b[i__ + b_dim1], ldb, (ftnlen)1);
#line 326 "AB05SD.f"
/* L10: */
#line 326 "AB05SD.f"
		}

#line 328 "AB05SD.f"
	    } else {
#line 329 "AB05SD.f"
		dlacpy_("F", n, m, &b[b_offset], ldb, &dwork[1], &ldwn, (
			ftnlen)1);
#line 330 "AB05SD.f"
		dgemm_("N", "N", n, p, m, alpha, &dwork[1], &ldwn, &d__[
			d_offset], ldd, &c_b11, &b[b_offset], ldb, (ftnlen)1, 
			(ftnlen)1);
#line 332 "AB05SD.f"
	    }
#line 333 "AB05SD.f"
	}
#line 334 "AB05SD.f"
    } else {

#line 336 "AB05SD.f"
	if (*ldwork < *n * *p) {

/*           Not enough working space for using DGEMM. */

#line 340 "AB05SD.f"
	    i__1 = *n;
#line 340 "AB05SD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 341 "AB05SD.f"
		dgemv_("N", m, p, alpha, &f[f_offset], ldf, &c__[i__ * c_dim1 
			+ 1], &c__1, &c_b14, &dwork[1], &c__1, (ftnlen)1);
#line 343 "AB05SD.f"
		dgemv_("N", n, m, &c_b11, &b[b_offset], ldb, &dwork[1], &c__1,
			 &c_b11, &a[i__ * a_dim1 + 1], &c__1, (ftnlen)1);
#line 345 "AB05SD.f"
/* L20: */
#line 345 "AB05SD.f"
	    }

#line 347 "AB05SD.f"
	    if (ljobd) {

#line 349 "AB05SD.f"
		i__1 = *n;
#line 349 "AB05SD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 350 "AB05SD.f"
		    dgemv_("T", m, p, alpha, &f[f_offset], ldf, &b[i__ + 
			    b_dim1], ldb, &c_b14, &dwork[1], &c__1, (ftnlen)1)
			    ;
#line 352 "AB05SD.f"
		    dgemv_("T", p, m, &c_b11, &d__[d_offset], ldd, &dwork[1], 
			    &c__1, &c_b11, &b[i__ + b_dim1], ldb, (ftnlen)1);
#line 354 "AB05SD.f"
/* L30: */
#line 354 "AB05SD.f"
		}

#line 356 "AB05SD.f"
	    }
#line 357 "AB05SD.f"
	} else {

#line 359 "AB05SD.f"
	    dgemm_("N", "N", n, p, m, alpha, &b[b_offset], ldb, &f[f_offset], 
		    ldf, &c_b14, &dwork[1], &ldwn, (ftnlen)1, (ftnlen)1);
#line 361 "AB05SD.f"
	    dgemm_("N", "N", n, n, p, &c_b11, &dwork[1], &ldwn, &c__[c_offset]
		    , ldc, &c_b11, &a[a_offset], lda, (ftnlen)1, (ftnlen)1);
#line 363 "AB05SD.f"
	    if (ljobd) {
#line 363 "AB05SD.f"
		dgemm_("N", "N", n, m, p, &c_b11, &dwork[1], &ldwn, &d__[
			d_offset], ldd, &c_b11, &b[b_offset], ldb, (ftnlen)1, 
			(ftnlen)1);
#line 363 "AB05SD.f"
	    }
#line 366 "AB05SD.f"
	}
#line 367 "AB05SD.f"
    }

#line 369 "AB05SD.f"
    return 0;
/* *** Last line of AB05SD *** */
} /* ab05sd_ */

