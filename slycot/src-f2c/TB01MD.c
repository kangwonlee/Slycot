#line 1 "TB01MD.f"
/* TB01MD.f -- translated by f2c (version 20100827).
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

#line 1 "TB01MD.f"
/* Table of constant values */

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;
static integer c__1 = 1;

/* Subroutine */ int tb01md_(char *jobu, char *uplo, integer *n, integer *m, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	u, integer *ldu, doublereal *dwork, integer *info, ftnlen jobu_len, 
	ftnlen uplo_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, u_dim1, u_offset, i__1, i__2, 
	    i__3;

    /* Local variables */
    static integer j, m1, n1, ii, nj;
    static doublereal dz;
    static integer par1, par2, par3, par4, par5, par6;
    static logical ljoba, ljobi;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical luplo;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dlatzm_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen);


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

/*     To reduce the pair (B,A) to upper or lower controller Hessenberg */
/*     form using (and optionally accumulating) unitary state-space */
/*     transformations. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBU    CHARACTER*1 */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix U the unitary state-space transformations for */
/*             reducing the system, as follows: */
/*             = 'N':  Do not form U; */
/*             = 'I':  U is initialized to the unit matrix and the */
/*                     unitary transformation matrix U is returned; */
/*             = 'U':  The given matrix U is updated by the unitary */
/*                     transformations used in the reduction. */

/*     UPLO    CHARACTER*1 */
/*             Indicates whether the user wishes the pair (B,A) to be */
/*             reduced to upper or lower controller Hessenberg form as */
/*             follows: */
/*             = 'U':  Upper controller Hessenberg form; */
/*             = 'L':  Lower controller Hessenberg form. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The actual state dimension, i.e. the order of the */
/*             matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The actual input dimension, i.e. the number of columns of */
/*             the matrix B.  M >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state transition matrix A to be transformed. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the transformed state transition matrix U' * A * U. */
/*             The annihilated elements are set to zero. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input matrix B to be transformed. */
/*             On exit, the leading N-by-M part of this array contains */
/*             the transformed input matrix U' * B. */
/*             The annihilated elements are set to zero. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             On entry, if JOBU = 'U', then the leading N-by-N part of */
/*             this array must contain a given matrix U (e.g. from a */
/*             previous call to another SLICOT routine), and on exit, the */
/*             leading N-by-N part of this array contains the product of */
/*             the input matrix U and the state-space transformation */
/*             matrix which reduces the given pair to controller */
/*             Hessenberg form. */
/*             On exit, if JOBU = 'I', then the leading N-by-N part of */
/*             this array contains the matrix of accumulated unitary */
/*             similarity transformations which reduces the given pair */
/*             to controller Hessenberg form. */
/*             If JOBU = 'N', the array U is not referenced and can be */
/*             supplied as a dummy array (i.e. set parameter LDU = 1 and */
/*             declare this array to be U(1,1) in the calling program). */

/*     LDU     INTEGER */
/*             The leading dimension of array U. If JOBU = 'U' or */
/*             JOBU = 'I', LDU >= MAX(1,N); if JOBU = 'N', LDU >= 1. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (MAX(N,M-1)) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine computes a unitary state-space transformation U, which */
/*     reduces the pair (B,A) to one of the following controller */
/*     Hessenberg forms: */

/*                    |*  . . .  *|*  . . . . . .  *| */
/*                    |   .      .|.               .| */
/*                    |     .    .|.               .| */
/*                    |       .  .|.               .| */
/*       [U'B|U'AU] = |          *|.               .| N */
/*                    |           |*               .| */
/*                    |           |   .            .| */
/*                    |           |     .          .| */
/*                    |           |       .        .| */
/*                    |           |         * . .  *| */
/*                         M               N */

/*     if UPLO = 'U', or */

/*                    |*  . . *         |           | */
/*                    |.        .       |           | */
/*                    |.          .     |           | */
/*                    |.            .   |           | */
/*       [U'AU|U'B] = |.               *|           | N */
/*                    |.               .|*          | */
/*                    |.               .|.  .       | */
/*                    |.               .|.    .     | */
/*                    |.               .|.      .   | */
/*                    |*  . . . . . .  *|*  . . .  *| */
/*                            N               M */
/*     if UPLO = 'L'. */

/*     IF M >= N, then the matrix U'B is trapezoidal and U'AU is full. */

/*     REFERENCES */

/*     [1] Van Dooren, P. and Verhaegen, M.H.G. */
/*         On the use of unitary state-space transformations. */
/*         In : Contemporary Mathematics on Linear Algebra and its Role */
/*         in Systems Theory, 47, AMS, Providence, 1985. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O((N + M) x N**2) operations and is */
/*     backward stable (see [1]). */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TB01AD by M. Vanbegin, and */
/*     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium. */

/*     REVISIONS */

/*     February 1997. */

/*     KEYWORDS */

/*     Controllability, controller Hessenberg form, orthogonal */
/*     transformation, unitary transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 196 "TB01MD.f"
    /* Parameter adjustments */
#line 196 "TB01MD.f"
    a_dim1 = *lda;
#line 196 "TB01MD.f"
    a_offset = 1 + a_dim1;
#line 196 "TB01MD.f"
    a -= a_offset;
#line 196 "TB01MD.f"
    b_dim1 = *ldb;
#line 196 "TB01MD.f"
    b_offset = 1 + b_dim1;
#line 196 "TB01MD.f"
    b -= b_offset;
#line 196 "TB01MD.f"
    u_dim1 = *ldu;
#line 196 "TB01MD.f"
    u_offset = 1 + u_dim1;
#line 196 "TB01MD.f"
    u -= u_offset;
#line 196 "TB01MD.f"
    --dwork;
#line 196 "TB01MD.f"

#line 196 "TB01MD.f"
    /* Function Body */
#line 196 "TB01MD.f"
    *info = 0;
#line 197 "TB01MD.f"
    luplo = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);
#line 198 "TB01MD.f"
    ljobi = lsame_(jobu, "I", (ftnlen)1, (ftnlen)1);
#line 199 "TB01MD.f"
    ljoba = ljobi || lsame_(jobu, "U", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 203 "TB01MD.f"
    if (! ljoba && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
#line 204 "TB01MD.f"
	*info = -1;
#line 205 "TB01MD.f"
    } else if (! luplo && ! lsame_(uplo, "L", (ftnlen)1, (ftnlen)1)) {
#line 206 "TB01MD.f"
	*info = -2;
#line 207 "TB01MD.f"
    } else if (*n < 0) {
#line 208 "TB01MD.f"
	*info = -3;
#line 209 "TB01MD.f"
    } else if (*m < 0) {
#line 210 "TB01MD.f"
	*info = -4;
#line 211 "TB01MD.f"
    } else if (*lda < max(1,*n)) {
#line 212 "TB01MD.f"
	*info = -6;
#line 213 "TB01MD.f"
    } else if (*ldb < max(1,*n)) {
#line 214 "TB01MD.f"
	*info = -8;
#line 215 "TB01MD.f"
    } else if (! ljoba && *ldu < 1 || ljoba && *ldu < max(1,*n)) {
#line 217 "TB01MD.f"
	*info = -10;
#line 218 "TB01MD.f"
    }

#line 220 "TB01MD.f"
    if (*info != 0) {

/*        Error return */

#line 224 "TB01MD.f"
	i__1 = -(*info);
#line 224 "TB01MD.f"
	xerbla_("TB01MD", &i__1, (ftnlen)6);
#line 225 "TB01MD.f"
	return 0;
#line 226 "TB01MD.f"
    }

/*     Quick return if possible. */

#line 230 "TB01MD.f"
    if (*n == 0 || *m == 0) {
#line 230 "TB01MD.f"
	return 0;
#line 230 "TB01MD.f"
    }

#line 233 "TB01MD.f"
    m1 = *m + 1;
#line 234 "TB01MD.f"
    n1 = *n - 1;

#line 236 "TB01MD.f"
    if (ljobi) {

/*        Initialize U to the identity matrix. */

#line 240 "TB01MD.f"
	dlaset_("Full", n, n, &c_b9, &c_b10, &u[u_offset], ldu, (ftnlen)4);
#line 241 "TB01MD.f"
    }

/*     Perform transformations involving both B and A. */

#line 245 "TB01MD.f"
    i__1 = min(*m,n1);
#line 245 "TB01MD.f"
    for (j = 1; j <= i__1; ++j) {
#line 246 "TB01MD.f"
	nj = *n - j;
#line 247 "TB01MD.f"
	if (luplo) {
#line 248 "TB01MD.f"
	    par1 = j;
#line 249 "TB01MD.f"
	    par2 = j;
#line 250 "TB01MD.f"
	    par3 = j + 1;
#line 251 "TB01MD.f"
	    par4 = *m;
#line 252 "TB01MD.f"
	    par5 = *n;
#line 253 "TB01MD.f"
	} else {
#line 254 "TB01MD.f"
	    par1 = *m - j + 1;
#line 255 "TB01MD.f"
	    par2 = nj + 1;
#line 256 "TB01MD.f"
	    par3 = 1;
#line 257 "TB01MD.f"
	    par4 = *m - j;
#line 258 "TB01MD.f"
	    par5 = nj;
#line 259 "TB01MD.f"
	}

#line 261 "TB01MD.f"
	i__2 = nj + 1;
#line 261 "TB01MD.f"
	dlarfg_(&i__2, &b[par2 + par1 * b_dim1], &b[par3 + par1 * b_dim1], &
		c__1, &dz);

/*        Update A. */

#line 265 "TB01MD.f"
	i__2 = nj + 1;
#line 265 "TB01MD.f"
	dlatzm_("Left", &i__2, n, &b[par3 + par1 * b_dim1], &c__1, &dz, &a[
		par2 + a_dim1], &a[par3 + a_dim1], lda, &dwork[1], (ftnlen)4);
#line 267 "TB01MD.f"
	i__2 = nj + 1;
#line 267 "TB01MD.f"
	dlatzm_("Right", n, &i__2, &b[par3 + par1 * b_dim1], &c__1, &dz, &a[
		par2 * a_dim1 + 1], &a[par3 * a_dim1 + 1], lda, &dwork[1], (
		ftnlen)5);

#line 270 "TB01MD.f"
	if (ljoba) {

/*           Update U. */

#line 274 "TB01MD.f"
	    i__2 = nj + 1;
#line 274 "TB01MD.f"
	    dlatzm_("Right", n, &i__2, &b[par3 + par1 * b_dim1], &c__1, &dz, &
		    u[par2 * u_dim1 + 1], &u[par3 * u_dim1 + 1], ldu, &dwork[
		    1], (ftnlen)5);
#line 276 "TB01MD.f"
	}

#line 278 "TB01MD.f"
	if (j != *m) {

/*           Update B */

#line 282 "TB01MD.f"
	    i__2 = nj + 1;
#line 282 "TB01MD.f"
	    i__3 = par4 - par3 + 1;
#line 282 "TB01MD.f"
	    dlatzm_("Left", &i__2, &i__3, &b[par3 + par1 * b_dim1], &c__1, &
		    dz, &b[par2 + par3 * b_dim1], &b[par3 + par3 * b_dim1], 
		    ldb, &dwork[1], (ftnlen)4);
#line 284 "TB01MD.f"
	}

#line 286 "TB01MD.f"
	i__2 = par5;
#line 286 "TB01MD.f"
	for (ii = par3; ii <= i__2; ++ii) {
#line 287 "TB01MD.f"
	    b[ii + par1 * b_dim1] = 0.;
#line 288 "TB01MD.f"
/* L10: */
#line 288 "TB01MD.f"
	}

#line 290 "TB01MD.f"
/* L20: */
#line 290 "TB01MD.f"
    }

#line 292 "TB01MD.f"
    i__1 = n1;
#line 292 "TB01MD.f"
    for (j = m1; j <= i__1; ++j) {

/*        Perform next transformations only involving A. */

#line 296 "TB01MD.f"
	nj = *n - j;
#line 297 "TB01MD.f"
	if (luplo) {
#line 298 "TB01MD.f"
	    par1 = j - *m;
#line 299 "TB01MD.f"
	    par2 = j;
#line 300 "TB01MD.f"
	    par3 = j + 1;
#line 301 "TB01MD.f"
	    par4 = *n;
#line 302 "TB01MD.f"
	    par5 = j - *m + 1;
#line 303 "TB01MD.f"
	    par6 = *n;
#line 304 "TB01MD.f"
	} else {
#line 305 "TB01MD.f"
	    par1 = *n + m1 - j;
#line 306 "TB01MD.f"
	    par2 = nj + 1;
#line 307 "TB01MD.f"
	    par3 = 1;
#line 308 "TB01MD.f"
	    par4 = nj;
#line 309 "TB01MD.f"
	    par5 = 1;
#line 310 "TB01MD.f"
	    par6 = *n + *m - j;
#line 311 "TB01MD.f"
	}

#line 313 "TB01MD.f"
	i__2 = nj + 1;
#line 313 "TB01MD.f"
	dlarfg_(&i__2, &a[par2 + par1 * a_dim1], &a[par3 + par1 * a_dim1], &
		c__1, &dz);

/*        Update A. */

#line 317 "TB01MD.f"
	i__2 = nj + 1;
#line 317 "TB01MD.f"
	i__3 = par6 - par5 + 1;
#line 317 "TB01MD.f"
	dlatzm_("Left", &i__2, &i__3, &a[par3 + par1 * a_dim1], &c__1, &dz, &
		a[par2 + par5 * a_dim1], &a[par3 + par5 * a_dim1], lda, &
		dwork[1], (ftnlen)4);
#line 319 "TB01MD.f"
	i__2 = nj + 1;
#line 319 "TB01MD.f"
	dlatzm_("Right", n, &i__2, &a[par3 + par1 * a_dim1], &c__1, &dz, &a[
		par2 * a_dim1 + 1], &a[par3 * a_dim1 + 1], lda, &dwork[1], (
		ftnlen)5);

#line 322 "TB01MD.f"
	if (ljoba) {

/*           Update U. */

#line 326 "TB01MD.f"
	    i__2 = nj + 1;
#line 326 "TB01MD.f"
	    dlatzm_("Right", n, &i__2, &a[par3 + par1 * a_dim1], &c__1, &dz, &
		    u[par2 * u_dim1 + 1], &u[par3 * u_dim1 + 1], ldu, &dwork[
		    1], (ftnlen)5);
#line 328 "TB01MD.f"
	}

#line 330 "TB01MD.f"
	i__2 = par4;
#line 330 "TB01MD.f"
	for (ii = par3; ii <= i__2; ++ii) {
#line 331 "TB01MD.f"
	    a[ii + par1 * a_dim1] = 0.;
#line 332 "TB01MD.f"
/* L30: */
#line 332 "TB01MD.f"
	}

#line 334 "TB01MD.f"
/* L40: */
#line 334 "TB01MD.f"
    }

#line 336 "TB01MD.f"
    return 0;
/* *** Last line of TB01MD *** */
} /* tb01md_ */

