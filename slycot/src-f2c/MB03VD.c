#line 1 "MB03VD.f"
/* MB03VD.f -- translated by f2c (version 20100827).
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

#line 1 "MB03VD.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;

/* Subroutine */ int mb03vd_(integer *n, integer *p, integer *ilo, integer *
	ihi, doublereal *a, integer *lda1, integer *lda2, doublereal *tau, 
	integer *ldtau, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_dim2, a_offset, tau_dim1, tau_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, i1, i2, nh;
    extern /* Subroutine */ int mb04py_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen), dcopy_(integer *, doublereal *, integer *, doublereal *,
	     integer *);
    static doublereal dummy[1];
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), xerbla_(char *, integer *, ftnlen);


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

/*     To reduce a product of p real general matrices A = A_1*A_2*...*A_p */
/*     to upper Hessenberg form, H = H_1*H_2*...*H_p, where H_1 is */
/*     upper Hessenberg, and H_2, ..., H_p are upper triangular, by using */
/*     orthogonal similarity transformations on A, */

/*             Q_1' * A_1 * Q_2 = H_1, */
/*             Q_2' * A_2 * Q_3 = H_2, */
/*                    ... */
/*             Q_p' * A_p * Q_1 = H_p. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the square matrices A_1, A_2, ..., A_p. */
/*             N >= 0. */

/*     P       (input) INTEGER */
/*             The number of matrices in the product A_1*A_2*...*A_p. */
/*             P >= 1. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that all matrices A_j, j = 2, ..., p, are */
/*             already upper triangular in rows and columns 1:ILO-1 and */
/*             IHI+1:N, and A_1 is upper Hessenberg in rows and columns */
/*             1:ILO-1 and IHI+1:N, with A_1(ILO,ILO-1) = 0 (unless */
/*             ILO = 1), and A_1(IHI+1,IHI) = 0 (unless IHI = N). */
/*             If this is not the case, ILO and IHI should be set to 1 */
/*             and N, respectively. */
/*             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDA1,LDA2,P) */
/*             On entry, the leading N-by-N-by-P part of this array must */
/*             contain the matrices of factors to be reduced; */
/*             specifically, A(*,*,j) must contain A_j, j = 1, ..., p. */
/*             On exit, the leading N-by-N upper triangle and the first */
/*             subdiagonal of A(*,*,1) contain the upper Hessenberg */
/*             matrix H_1, and the elements below the first subdiagonal, */
/*             with the first column of the array TAU represent the */
/*             orthogonal matrix Q_1 as a product of elementary */
/*             reflectors. See FURTHER COMMENTS. */
/*             For j > 1, the leading N-by-N upper triangle of A(*,*,j) */
/*             contains the upper triangular matrix H_j, and the elements */
/*             below the diagonal, with the j-th column of the array TAU */
/*             represent the orthogonal matrix Q_j as a product of */
/*             elementary reflectors. See FURTHER COMMENTS. */

/*     LDA1    INTEGER */
/*             The first leading dimension of the array A. */
/*             LDA1 >= max(1,N). */

/*     LDA2    INTEGER */
/*             The second leading dimension of the array A. */
/*             LDA2 >= max(1,N). */

/*     TAU     (output) DOUBLE PRECISION array, dimension (LDTAU,P) */
/*             The leading N-1 elements in the j-th column contain the */
/*             scalar factors of the elementary reflectors used to form */
/*             the matrix Q_j, j = 1, ..., P. See FURTHER COMMENTS. */

/*     LDTAU   INTEGER */
/*             The leading dimension of the array TAU. */
/*             LDTAU >= max(1,N-1). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The algorithm consists in ihi-ilo major steps. In each such */
/*     step i, ilo <= i <= ihi-1, the subdiagonal elements in the i-th */
/*     column of A_j are annihilated using a Householder transformation */
/*     from the left, which is also applied to A_(j-1) from the right, */
/*     for j = p:-1:2. Then, the elements below the subdiagonal of the */
/*     i-th column of A_1 are annihilated, and the Householder */
/*     transformation is also applied to A_p from the right. */
/*     See FURTHER COMMENTS. */

/*     REFERENCES */

/*     [1] Bojanczyk, A.W., Golub, G. and Van Dooren, P. */
/*         The periodic Schur decomposition: algorithms and applications. */
/*         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42, */
/*         1992. */

/*     [2] Sreedhar, J. and Van Dooren, P. */
/*         Periodic Schur form and some matrix equations. */
/*         Proc. of the Symposium on the Mathematical Theory of Networks */
/*         and Systems (MTNS'93), Regensburg, Germany (U. Helmke, */
/*         R. Mennicken and J. Saurer, Eds.), Vol. 1, pp. 339-362, 1994. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable. */

/*     FURTHER COMMENTS */

/*     Each matrix Q_j is represented as a product of (ihi-ilo) */
/*     elementary reflectors, */

/*        Q_j = H_j(ilo) H_j(ilo+1) . . . H_j(ihi-1). */

/*     Each H_j(i), i = ilo, ..., ihi-1, has the form */

/*        H_j(i) = I - tau_j * v_j * v_j', */

/*     where tau_j is a real scalar, and v_j is a real vector with */
/*     v_j(1:i) = 0, v_j(i+1) = 1 and v_j(ihi+1:n) = 0; v_j(i+2:ihi) */
/*     is stored on exit in A_j(i+2:ihi,i), and tau_j in TAU(i,j). */

/*     The contents of A_1 are illustrated by the following example */
/*     for n = 7, ilo = 2, and ihi = 6: */

/*     on entry                         on exit */

/*     ( a   a   a   a   a   a   a )    ( a   h   h   h   h   h   a ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   h   h   h   h   h   a ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   h   h   h   h   h   h ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   v2  h   h   h   h   h ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  h   h   h   h ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  v4  h   h   h ) */
/*     ( 0   0   0   0   0   0   a )    ( 0   0   0   0   0   0   a ) */

/*     where a denotes an element of the original matrix A_1, h denotes */
/*     a modified element of the upper Hessenberg matrix H_1, and vi */
/*     denotes an element of the vector defining H_1(i). */

/*     The contents of A_j, j > 1, are illustrated by the following */
/*     example for n = 7, ilo = 2, and ihi = 6: */

/*     on entry                         on exit */

/*     ( a   a   a   a   a   a   a )    ( a   h   h   h   h   h   a ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   h   h   h   h   h   h ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   v2  h   h   h   h   h ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  h   h   h   h ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  v4  h   h   h ) */
/*     ( 0   a   a   a   a   a   a )    ( 0   v2  v3  v4  v5  h   h ) */
/*     ( 0   0   0   0   0   0   a )    ( 0   0   0   0   0   0   a ) */

/*     where a denotes an element of the original matrix A_j, h denotes */
/*     a modified element of the upper triangular matrix H_j, and vi */
/*     denotes an element of the vector defining H_j(i). (The element */
/*     (1,2) in A_p is also unchanged for this example.) */

/*     Note that for P = 1, the LAPACK Library routine DGEHRD could be */
/*     more efficient on some computer architectures than this routine */
/*     (a BLAS 2 version). */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, and A. Varga, */
/*     German Aerospace Center, DLR Oberpfaffenhofen, February 1999. */
/*     Partly based on the routine PSHESS by A. Varga */
/*     (DLR Oberpfaffenhofen), November 26, 1995. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Hessenberg form, orthogonal transformation, periodic systems, */
/*     similarity transformation, triangular form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 228 "MB03VD.f"
    /* Parameter adjustments */
#line 228 "MB03VD.f"
    a_dim1 = *lda1;
#line 228 "MB03VD.f"
    a_dim2 = *lda2;
#line 228 "MB03VD.f"
    a_offset = 1 + a_dim1 * (1 + a_dim2);
#line 228 "MB03VD.f"
    a -= a_offset;
#line 228 "MB03VD.f"
    tau_dim1 = *ldtau;
#line 228 "MB03VD.f"
    tau_offset = 1 + tau_dim1;
#line 228 "MB03VD.f"
    tau -= tau_offset;
#line 228 "MB03VD.f"
    --dwork;
#line 228 "MB03VD.f"

#line 228 "MB03VD.f"
    /* Function Body */
#line 228 "MB03VD.f"
    *info = 0;
#line 229 "MB03VD.f"
    if (*n < 0) {
#line 230 "MB03VD.f"
	*info = -1;
#line 231 "MB03VD.f"
    } else if (*p < 1) {
#line 232 "MB03VD.f"
	*info = -2;
#line 233 "MB03VD.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 234 "MB03VD.f"
	*info = -3;
#line 235 "MB03VD.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 236 "MB03VD.f"
	*info = -4;
#line 237 "MB03VD.f"
    } else if (*lda1 < max(1,*n)) {
#line 238 "MB03VD.f"
	*info = -6;
#line 239 "MB03VD.f"
    } else if (*lda2 < max(1,*n)) {
#line 240 "MB03VD.f"
	*info = -7;
#line 241 "MB03VD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 241 "MB03VD.f"
	i__1 = 1, i__2 = *n - 1;
#line 241 "MB03VD.f"
	if (*ldtau < max(i__1,i__2)) {
#line 242 "MB03VD.f"
	    *info = -9;
#line 243 "MB03VD.f"
	}
#line 243 "MB03VD.f"
    }
#line 244 "MB03VD.f"
    if (*info != 0) {
#line 245 "MB03VD.f"
	i__1 = -(*info);
#line 245 "MB03VD.f"
	xerbla_("MB03VD", &i__1, (ftnlen)6);
#line 246 "MB03VD.f"
	return 0;
#line 247 "MB03VD.f"
    }

/*     Quick return if possible. */

#line 251 "MB03VD.f"
    nh = *ihi - *ilo + 1;
#line 252 "MB03VD.f"
    if (nh <= 1) {
#line 252 "MB03VD.f"
	return 0;
#line 252 "MB03VD.f"
    }

#line 255 "MB03VD.f"
    dummy[0] = 0.;

#line 257 "MB03VD.f"
    i__1 = *ihi - 1;
#line 257 "MB03VD.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 258 "MB03VD.f"
	i1 = i__ + 1;
/* Computing MIN */
#line 259 "MB03VD.f"
	i__2 = i__ + 2;
#line 259 "MB03VD.f"
	i2 = min(i__2,*n);

#line 261 "MB03VD.f"
	for (j = *p; j >= 2; --j) {

/*           Set the elements 1:ILO-1 and IHI:N-1 of TAU(*,J) to zero. */

#line 265 "MB03VD.f"
	    i__2 = *ilo - 1;
#line 265 "MB03VD.f"
	    dcopy_(&i__2, dummy, &c__0, &tau[j * tau_dim1 + 1], &c__1);
#line 266 "MB03VD.f"
	    if (*ihi < *n) {
#line 266 "MB03VD.f"
		i__2 = *n - *ihi;
#line 266 "MB03VD.f"
		dcopy_(&i__2, dummy, &c__0, &tau[*ihi + j * tau_dim1], &c__1);
#line 266 "MB03VD.f"
	    }

/*           Compute elementary reflector H_j(i) to annihilate */
/*           A_j(i+1:ihi,i). */

#line 272 "MB03VD.f"
	    i__2 = *ihi - i__ + 1;
#line 272 "MB03VD.f"
	    dlarfg_(&i__2, &a[i__ + (i__ + j * a_dim2) * a_dim1], &a[i1 + (
		    i__ + j * a_dim2) * a_dim1], &c__1, &tau[i__ + j * 
		    tau_dim1]);

/*           Apply H_j(i) to A_(j-1)(1:ihi,i:ihi) from the right. */

#line 277 "MB03VD.f"
	    i__2 = *ihi - i__ + 1;
#line 277 "MB03VD.f"
	    mb04py_("Right", ihi, &i__2, &a[i1 + (i__ + j * a_dim2) * a_dim1],
		     &tau[i__ + j * tau_dim1], &a[(i__ + (j - 1) * a_dim2) * 
		    a_dim1 + 1], lda1, &dwork[1], (ftnlen)5);

/*           Apply H_j(i) to A_j(i:ihi,i+1:n) from the left. */

#line 282 "MB03VD.f"
	    i__2 = *ihi - i__ + 1;
#line 282 "MB03VD.f"
	    i__3 = *n - i__;
#line 282 "MB03VD.f"
	    mb04py_("Left", &i__2, &i__3, &a[i1 + (i__ + j * a_dim2) * a_dim1]
		    , &tau[i__ + j * tau_dim1], &a[i__ + (i1 + j * a_dim2) * 
		    a_dim1], lda1, &dwork[1], (ftnlen)4);
#line 284 "MB03VD.f"
/* L10: */
#line 284 "MB03VD.f"
	}

/*        Compute elementary reflector H_1(i) to annihilate */
/*        A_1(i+2:ihi,i). */

#line 289 "MB03VD.f"
	i__2 = *ihi - i__;
#line 289 "MB03VD.f"
	dlarfg_(&i__2, &a[i1 + (i__ + a_dim2) * a_dim1], &a[i2 + (i__ + 
		a_dim2) * a_dim1], &c__1, &tau[i__ + tau_dim1]);

/*        Apply H_1(i) to A_p(1:ihi,i+1:ihi) from the right. */

#line 294 "MB03VD.f"
	i__2 = *ihi - i__;
#line 294 "MB03VD.f"
	mb04py_("Right", ihi, &i__2, &a[i2 + (i__ + a_dim2) * a_dim1], &tau[
		i__ + tau_dim1], &a[(i1 + *p * a_dim2) * a_dim1 + 1], lda1, &
		dwork[1], (ftnlen)5);

/*        Apply H_1(i) to A_1(i+1:ihi,i+1:n) from the left. */

#line 299 "MB03VD.f"
	i__2 = *ihi - i__;
#line 299 "MB03VD.f"
	i__3 = *n - i__;
#line 299 "MB03VD.f"
	mb04py_("Left", &i__2, &i__3, &a[i2 + (i__ + a_dim2) * a_dim1], &tau[
		i__ + tau_dim1], &a[i1 + (i1 + a_dim2) * a_dim1], lda1, &
		dwork[1], (ftnlen)4);
#line 301 "MB03VD.f"
/* L20: */
#line 301 "MB03VD.f"
    }

#line 303 "MB03VD.f"
    return 0;

/* *** Last line of MB03VD *** */
} /* mb03vd_ */

