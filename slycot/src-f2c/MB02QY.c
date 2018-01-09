#line 1 "MB02QY.f"
/* MB02QY.f -- translated by f2c (version 20100827).
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

#line 1 "MB02QY.f"
/* Table of constant values */

static integer c__0 = 0;
static doublereal c_b15 = 0.;
static doublereal c_b28 = 1.;
static integer c__1 = 1;

/* Subroutine */ int mb02qy_(integer *m, integer *n, integer *nrhs, integer *
	rank, doublereal *a, integer *lda, integer *jpvt, doublereal *b, 
	integer *ldb, doublereal *tau, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, mn;
    static doublereal anrm, bnrm;
    static integer iascl, ibscl;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dtrsm_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen), dlaset_(char *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static doublereal bignum;
    extern doublereal dlantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    static doublereal maxwrk, smlnum;
    extern /* Subroutine */ int dormrz_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dtzrzf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);


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

/*     To determine the minimum-norm solution to a real linear least */
/*     squares problem: */

/*         minimize || A * X - B ||, */

/*     using the rank-revealing QR factorization of a real general */
/*     M-by-N matrix  A,  computed by SLICOT Library routine  MB03OD. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrices A and B.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     NRHS    (input) INTEGER */
/*             The number of columns of the matrix B.  NRHS >= 0. */

/*     RANK    (input) INTEGER */
/*             The effective rank of  A,  as returned by SLICOT Library */
/*             routine  MB03OD.  min(M,N) >= RANK >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDA, N ) */
/*             On entry, the leading min(M,N)-by-N upper trapezoidal */
/*             part of this array contains the triangular factor  R,  as */
/*             returned by SLICOT Library routine  MB03OD.  The strict */
/*             lower trapezoidal part of  A  is not referenced. */
/*             On exit, if  RANK < N,  the leading  RANK-by-RANK  upper */
/*             triangular part of this array contains the upper */
/*             triangular matrix  R  of the complete orthogonal */
/*             factorization of  A,  and the submatrix  (1:RANK,RANK+1:N) */
/*             of this array, with the array  TAU,  represent the */
/*             orthogonal matrix  Z  (of the complete orthogonal */
/*             factorization of  A),  as a product of  RANK  elementary */
/*             reflectors. */
/*             On exit, if  RANK = N,  this array is unchanged. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     JPVT    (input) INTEGER array, dimension ( N ) */
/*             The recorded permutations performed by SLICOT Library */
/*             routine  MB03OD;  if  JPVT(i) = k,  then the i-th column */
/*             of  A*P  was the k-th column of the original matrix  A. */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDB, NRHS ) */
/*             On entry, if  NRHS > 0,  the leading M-by-NRHS part of */
/*             this array must contain the matrix  B  (corresponding to */
/*             the transformed matrix  A,  returned by SLICOT Library */
/*             routine  MB03OD). */
/*             On exit, if  NRHS > 0,  the leading N-by-NRHS part of this */
/*             array contains the solution matrix X. */
/*             If  M >= N  and  RANK = N,  the residual sum-of-squares */
/*             for the solution in the i-th column is given by the sum */
/*             of squares of elements  N+1:M  in that column. */
/*             If  NRHS = 0,  the array  B  is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= max(1,M,N),  if  NRHS > 0. */
/*             LDB >= 1,           if  NRHS = 0. */

/*     TAU     (output) DOUBLE PRECISION array, dimension ( min(M,N) ) */
/*             The scalar factors of the elementary reflectors. */
/*             If  RANK = N,  the array  TAU  is not referenced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension ( LDWORK ) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( 1, N, NRHS ). */
/*             For good performance,  LDWORK  should sometimes be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine uses a QR factorization with column pivoting: */

/*        A * P = Q * R = Q * [ R11 R12 ], */
/*                            [  0  R22 ] */

/*     where  R11  is an upper triangular submatrix of estimated rank */
/*     RANK,  the effective rank of  A.  The submatrix  R22  can be */
/*     considered as negligible. */

/*     If  RANK < N,  then  R12  is annihilated by orthogonal */
/*     transformations from the right, arriving at the complete */
/*     orthogonal factorization: */

/*        A * P = Q * [ T11 0 ] * Z. */
/*                    [  0  0 ] */

/*     The minimum-norm solution is then */

/*        X = P * Z' [ inv(T11)*Q1'*B ], */
/*                   [        0       ] */

/*     where Q1 consists of the first  RANK  columns of Q. */

/*     The input data for  MB02QY  are the transformed matrices  Q' * A */
/*     (returned by SLICOT Library routine  MB03OD)  and  Q' * B. */
/*     Matrix  Q  is not needed. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Least squares solutions; QR decomposition. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 182 "MB02QY.f"
    /* Parameter adjustments */
#line 182 "MB02QY.f"
    a_dim1 = *lda;
#line 182 "MB02QY.f"
    a_offset = 1 + a_dim1;
#line 182 "MB02QY.f"
    a -= a_offset;
#line 182 "MB02QY.f"
    --jpvt;
#line 182 "MB02QY.f"
    b_dim1 = *ldb;
#line 182 "MB02QY.f"
    b_offset = 1 + b_dim1;
#line 182 "MB02QY.f"
    b -= b_offset;
#line 182 "MB02QY.f"
    --tau;
#line 182 "MB02QY.f"
    --dwork;
#line 182 "MB02QY.f"

#line 182 "MB02QY.f"
    /* Function Body */
#line 182 "MB02QY.f"
    mn = min(*m,*n);

/*     Test the input scalar arguments. */

#line 186 "MB02QY.f"
    *info = 0;
#line 187 "MB02QY.f"
    if (*m < 0) {
#line 188 "MB02QY.f"
	*info = -1;
#line 189 "MB02QY.f"
    } else if (*n < 0) {
#line 190 "MB02QY.f"
	*info = -2;
#line 191 "MB02QY.f"
    } else if (*nrhs < 0) {
#line 192 "MB02QY.f"
	*info = -3;
#line 193 "MB02QY.f"
    } else if (*rank < 0 || *rank > mn) {
#line 194 "MB02QY.f"
	*info = -4;
#line 195 "MB02QY.f"
    } else if (*lda < max(1,*m)) {
#line 196 "MB02QY.f"
	*info = -6;
#line 197 "MB02QY.f"
    } else if (*ldb < 1 || *nrhs > 0 && *ldb < max(*m,*n)) {
#line 199 "MB02QY.f"
	*info = -9;
#line 200 "MB02QY.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 200 "MB02QY.f"
	i__1 = max(1,*n);
#line 200 "MB02QY.f"
	if (*ldwork < max(i__1,*nrhs)) {
#line 201 "MB02QY.f"
	    *info = -12;
#line 202 "MB02QY.f"
	}
#line 202 "MB02QY.f"
    }

#line 204 "MB02QY.f"
    if (*info != 0) {
#line 205 "MB02QY.f"
	i__1 = -(*info);
#line 205 "MB02QY.f"
	xerbla_("MB02QY", &i__1, (ftnlen)6);
#line 206 "MB02QY.f"
	return 0;
#line 207 "MB02QY.f"
    }

/*     Quick return if possible. */

#line 211 "MB02QY.f"
    if (min(mn,*nrhs) == 0) {
#line 212 "MB02QY.f"
	dwork[1] = 1.;
#line 213 "MB02QY.f"
	return 0;
#line 214 "MB02QY.f"
    }

/*     Logically partition R = [ R11 R12 ], */
/*                             [  0  R22 ] */

/*     where R11 = R(1:RANK,1:RANK).  If  RANK = N,  let  T11 = R11. */

#line 221 "MB02QY.f"
    maxwrk = (doublereal) (*n);
#line 222 "MB02QY.f"
    if (*rank < *n) {

/*        Get machine parameters. */

#line 226 "MB02QY.f"
	smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
		ftnlen)9);
#line 227 "MB02QY.f"
	bignum = 1. / smlnum;
#line 228 "MB02QY.f"
	dlabad_(&smlnum, &bignum);

/*        Scale A, B if max entries outside range [SMLNUM,BIGNUM]. */

#line 232 "MB02QY.f"
	anrm = dlantr_("MaxNorm", "Upper", "Non-unit", rank, n, &a[a_offset], 
		lda, &dwork[1], (ftnlen)7, (ftnlen)5, (ftnlen)8);
#line 234 "MB02QY.f"
	iascl = 0;
#line 235 "MB02QY.f"
	if (anrm > 0. && anrm < smlnum) {

/*           Scale matrix norm up to SMLNUM. */

#line 239 "MB02QY.f"
	    dlascl_("Upper", &c__0, &c__0, &anrm, &smlnum, rank, n, &a[
		    a_offset], lda, info, (ftnlen)5);
#line 241 "MB02QY.f"
	    iascl = 1;
#line 242 "MB02QY.f"
	} else if (anrm > bignum) {

/*           Scale matrix norm down to BIGNUM. */

#line 246 "MB02QY.f"
	    dlascl_("Upper", &c__0, &c__0, &anrm, &bignum, rank, n, &a[
		    a_offset], lda, info, (ftnlen)5);
#line 248 "MB02QY.f"
	    iascl = 2;
#line 249 "MB02QY.f"
	} else if (anrm == 0.) {

/*           Matrix all zero. Return zero solution. */

#line 253 "MB02QY.f"
	    dlaset_("Full", n, nrhs, &c_b15, &c_b15, &b[b_offset], ldb, (
		    ftnlen)4);
#line 254 "MB02QY.f"
	    dwork[1] = 1.;
#line 255 "MB02QY.f"
	    return 0;
#line 256 "MB02QY.f"
	}

#line 258 "MB02QY.f"
	bnrm = dlange_("MaxNorm", m, nrhs, &b[b_offset], ldb, &dwork[1], (
		ftnlen)7);
#line 259 "MB02QY.f"
	ibscl = 0;
#line 260 "MB02QY.f"
	if (bnrm > 0. && bnrm < smlnum) {

/*           Scale matrix norm up to SMLNUM. */

#line 264 "MB02QY.f"
	    dlascl_("General", &c__0, &c__0, &bnrm, &smlnum, m, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)7);
#line 266 "MB02QY.f"
	    ibscl = 1;
#line 267 "MB02QY.f"
	} else if (bnrm > bignum) {

/*           Scale matrix norm down to BIGNUM. */

#line 271 "MB02QY.f"
	    dlascl_("General", &c__0, &c__0, &bnrm, &bignum, m, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)7);
#line 273 "MB02QY.f"
	    ibscl = 2;
#line 274 "MB02QY.f"
	}

/*        [R11,R12] = [ T11, 0 ] * Z. */
/*        Details of Householder rotations are stored in TAU. */
/*        Workspace need RANK, prefer RANK*NB. */

#line 280 "MB02QY.f"
	dtzrzf_(rank, n, &a[a_offset], lda, &tau[1], &dwork[1], ldwork, info);
#line 281 "MB02QY.f"
	maxwrk = max(maxwrk,dwork[1]);
#line 282 "MB02QY.f"
    }

/*     B(1:RANK,1:NRHS) := inv(T11) * B(1:RANK,1:NRHS). */

#line 286 "MB02QY.f"
    dtrsm_("Left", "Upper", "No transpose", "Non-unit", rank, nrhs, &c_b28, &
	    a[a_offset], lda, &b[b_offset], ldb, (ftnlen)4, (ftnlen)5, (
	    ftnlen)12, (ftnlen)8);

#line 289 "MB02QY.f"
    if (*rank < *n) {

#line 291 "MB02QY.f"
	i__1 = *n - *rank;
#line 291 "MB02QY.f"
	dlaset_("Full", &i__1, nrhs, &c_b15, &c_b15, &b[*rank + 1 + b_dim1], 
		ldb, (ftnlen)4);

/*        B(1:N,1:NRHS) := Z' * B(1:N,1:NRHS). */
/*        Workspace need NRHS, prefer NRHS*NB. */

#line 297 "MB02QY.f"
	i__1 = *n - *rank;
#line 297 "MB02QY.f"
	dormrz_("Left", "Transpose", n, nrhs, rank, &i__1, &a[a_offset], lda, 
		&tau[1], &b[b_offset], ldb, &dwork[1], ldwork, info, (ftnlen)
		4, (ftnlen)9);
#line 299 "MB02QY.f"
	maxwrk = max(maxwrk,dwork[1]);

/*        Undo scaling. */

#line 303 "MB02QY.f"
	if (iascl == 1) {
#line 304 "MB02QY.f"
	    dlascl_("General", &c__0, &c__0, &anrm, &smlnum, n, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)7);
#line 306 "MB02QY.f"
	    dlascl_("Upper", &c__0, &c__0, &smlnum, &anrm, rank, rank, &a[
		    a_offset], lda, info, (ftnlen)5);
#line 308 "MB02QY.f"
	} else if (iascl == 2) {
#line 309 "MB02QY.f"
	    dlascl_("General", &c__0, &c__0, &anrm, &bignum, n, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)7);
#line 311 "MB02QY.f"
	    dlascl_("Upper", &c__0, &c__0, &bignum, &anrm, rank, rank, &a[
		    a_offset], lda, info, (ftnlen)5);
#line 313 "MB02QY.f"
	}
#line 314 "MB02QY.f"
	if (ibscl == 1) {
#line 315 "MB02QY.f"
	    dlascl_("General", &c__0, &c__0, &smlnum, &bnrm, n, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)7);
#line 317 "MB02QY.f"
	} else if (ibscl == 2) {
#line 318 "MB02QY.f"
	    dlascl_("General", &c__0, &c__0, &bignum, &bnrm, n, nrhs, &b[
		    b_offset], ldb, info, (ftnlen)7);
#line 320 "MB02QY.f"
	}
#line 321 "MB02QY.f"
    }

/*     B(1:N,1:NRHS) := P * B(1:N,1:NRHS). */
/*     Workspace N. */

#line 326 "MB02QY.f"
    i__1 = *nrhs;
#line 326 "MB02QY.f"
    for (j = 1; j <= i__1; ++j) {

#line 328 "MB02QY.f"
	i__2 = *n;
#line 328 "MB02QY.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 329 "MB02QY.f"
	    dwork[jpvt[i__]] = b[i__ + j * b_dim1];
#line 330 "MB02QY.f"
/* L10: */
#line 330 "MB02QY.f"
	}

#line 332 "MB02QY.f"
	dcopy_(n, &dwork[1], &c__1, &b[j * b_dim1 + 1], &c__1);
#line 333 "MB02QY.f"
/* L20: */
#line 333 "MB02QY.f"
    }

#line 335 "MB02QY.f"
    dwork[1] = maxwrk;
#line 336 "MB02QY.f"
    return 0;

/* *** Last line of MB02QY *** */
} /* mb02qy_ */

