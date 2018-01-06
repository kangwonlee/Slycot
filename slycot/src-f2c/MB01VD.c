#line 1 "MB01VD.f"
/* MB01VD.f -- translated by f2c (version 20100827).
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

#line 1 "MB01VD.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int mb01vd_(char *trana, char *tranb, integer *ma, integer *
	na, integer *mb, integer *nb, doublereal *alpha, doublereal *beta, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, integer *mc, integer *nc, integer *info, ftnlen 
	trana_len, ftnlen tranb_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, j, k, l, ic, jc, lc, nz;
    static doublereal aij, dum[1];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static logical transa, transb, sparse;


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

/*     To perform the following matrix operation */

/*        C = alpha*kron( op(A), op(B) ) + beta*C, */

/*     where alpha and beta are real scalars, op(M) is either matrix M or */
/*     its transpose, M', and kron( X, Y ) denotes the Kronecker product */
/*     of the matrices X and Y. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TRANA   CHARACTER*1 */
/*             Specifies the form of op(A) to be used as follows: */
/*             = 'N':  op(A) = A; */
/*             = 'T':  op(A) = A'; */
/*             = 'C':  op(A) = A'. */

/*     TRANB   CHARACTER*1 */
/*             Specifies the form of op(B) to be used as follows: */
/*             = 'N':  op(B) = B; */
/*             = 'T':  op(B) = B'; */
/*             = 'C':  op(B) = B'. */

/*     Input/Output Parameters */

/*     MA      (input) INTEGER */
/*             The number of rows of the matrix op(A).  MA >= 0. */

/*     NA      (input) INTEGER */
/*             The number of columns of the matrix op(A).  NA >= 0. */

/*     MB      (input) INTEGER */
/*             The number of rows of the matrix op(B).  MB >= 0. */

/*     NB      (input) INTEGER */
/*             The number of columns of the matrix op(B).  NB >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             The scalar alpha. When alpha is zero then A and B need not */
/*             be set before entry. */

/*     BETA    (input) DOUBLE PRECISION */
/*             The scalar beta. When beta is zero then C need not be */
/*             set before entry. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,ka), */
/*             where ka is NA when TRANA = 'N', and is MA otherwise. */
/*             If TRANA = 'N', the leading MA-by-NA part of this array */
/*             must contain the matrix A; otherwise, the leading NA-by-MA */
/*             part of this array must contain the matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= max(1,MA), if TRANA = 'N'; */
/*             LDA >= max(1,NA), if TRANA = 'T' or 'C'. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,kb) */
/*             where kb is NB when TRANB = 'N', and is MB otherwise. */
/*             If TRANB = 'N', the leading MB-by-NB part of this array */
/*             must contain the matrix B; otherwise, the leading NB-by-MB */
/*             part of this array must contain the matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= max(1,MB), if TRANB = 'N'; */
/*             LDB >= max(1,NB), if TRANB = 'T' or 'C'. */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,NC) */
/*             On entry, if beta is nonzero, the leading MC-by-NC part of */
/*             this array must contain the given matric C, where */
/*             MC = MA*MB and NC = NA*NB. */
/*             On exit, the leading MC-by-NC part of this array contains */
/*             the computed matrix expression */
/*             C = alpha*kron( op(A), op(B) ) + beta*C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= max(1,MC). */

/*     MC      (output) INTEGER */
/*             The number of rows of the matrix C.  MC = MA*MB. */

/*     NC      (output) INTEGER */
/*             The number of columns of the matrix C.  NC = NA*NB. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The Kronecker product of the matrices op(A) and op(B) is computed */
/*     column by column. */

/*     FURTHER COMMENTS */

/*     The multiplications by zero elements in A are avoided, if the */
/*     matrix A is considered to be sparse, i.e., if */
/*     (number of zeros in A)/(MA*NA) >= SPARST = 0.8. The code makes */
/*     NB+1 passes through the matrix A, and MA*NA passes through the */
/*     matrix B. If LDA and/or LDB are very large, and op(A) = A' and/or */
/*     op(B) = B', it could be more efficient to transpose A and/or B */
/*     before calling this routine, and use the 'N' values for TRANA */
/*     and/or TRANB. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, February 2000. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 176 "MB01VD.f"
    /* Parameter adjustments */
#line 176 "MB01VD.f"
    a_dim1 = *lda;
#line 176 "MB01VD.f"
    a_offset = 1 + a_dim1;
#line 176 "MB01VD.f"
    a -= a_offset;
#line 176 "MB01VD.f"
    b_dim1 = *ldb;
#line 176 "MB01VD.f"
    b_offset = 1 + b_dim1;
#line 176 "MB01VD.f"
    b -= b_offset;
#line 176 "MB01VD.f"
    c_dim1 = *ldc;
#line 176 "MB01VD.f"
    c_offset = 1 + c_dim1;
#line 176 "MB01VD.f"
    c__ -= c_offset;
#line 176 "MB01VD.f"

#line 176 "MB01VD.f"
    /* Function Body */
#line 176 "MB01VD.f"
    transa = lsame_(trana, "T", (ftnlen)1, (ftnlen)1) || lsame_(trana, "C", (
	    ftnlen)1, (ftnlen)1);
#line 177 "MB01VD.f"
    transb = lsame_(tranb, "T", (ftnlen)1, (ftnlen)1) || lsame_(tranb, "C", (
	    ftnlen)1, (ftnlen)1);
#line 178 "MB01VD.f"
    *mc = *ma * *mb;
#line 179 "MB01VD.f"
    *info = 0;
#line 180 "MB01VD.f"
    if (! (transa || lsame_(trana, "N", (ftnlen)1, (ftnlen)1))) {
#line 181 "MB01VD.f"
	*info = -1;
#line 182 "MB01VD.f"
    } else if (! (transb || lsame_(tranb, "N", (ftnlen)1, (ftnlen)1))) {
#line 183 "MB01VD.f"
	*info = -2;
#line 184 "MB01VD.f"
    } else if (*ma < 0) {
#line 185 "MB01VD.f"
	*info = -3;
#line 186 "MB01VD.f"
    } else if (*na < 0) {
#line 187 "MB01VD.f"
	*info = -4;
#line 188 "MB01VD.f"
    } else if (*mb < 0) {
#line 189 "MB01VD.f"
	*info = -5;
#line 190 "MB01VD.f"
    } else if (*nb < 0) {
#line 191 "MB01VD.f"
	*info = -6;
#line 192 "MB01VD.f"
    } else if (transa && *lda < *na || *lda < 1 || ! transa && *lda < *ma) {
#line 194 "MB01VD.f"
	*info = -10;
#line 195 "MB01VD.f"
    } else if (transb && *ldb < *nb || *ldb < 1 || ! transb && *ldb < *mb) {
#line 197 "MB01VD.f"
	*info = -12;
#line 198 "MB01VD.f"
    } else if (*ldc < max(1,*mc)) {
#line 199 "MB01VD.f"
	*info = -14;
#line 200 "MB01VD.f"
    }

#line 202 "MB01VD.f"
    if (*info != 0) {

/*        Error return. */

#line 206 "MB01VD.f"
	i__1 = -(*info);
#line 206 "MB01VD.f"
	xerbla_("MB01VD", &i__1, (ftnlen)6);
#line 207 "MB01VD.f"
	return 0;
#line 208 "MB01VD.f"
    }

/*     Quick return, if possible. */

#line 212 "MB01VD.f"
    *nc = *na * *nb;
#line 213 "MB01VD.f"
    if (*mc == 0 || *nc == 0) {
#line 213 "MB01VD.f"
	return 0;
#line 213 "MB01VD.f"
    }

#line 216 "MB01VD.f"
    if (*alpha == 0.) {
#line 217 "MB01VD.f"
	if (*beta == 0.) {
#line 218 "MB01VD.f"
	    dlaset_("Full", mc, nc, &c_b10, &c_b10, &c__[c_offset], ldc, (
		    ftnlen)4);
#line 219 "MB01VD.f"
	} else if (*beta != 1.) {

#line 221 "MB01VD.f"
	    i__1 = *nc;
#line 221 "MB01VD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 222 "MB01VD.f"
		dscal_(mc, beta, &c__[j * c_dim1 + 1], &c__1);
#line 223 "MB01VD.f"
/* L10: */
#line 223 "MB01VD.f"
	    }

#line 225 "MB01VD.f"
	}
#line 226 "MB01VD.f"
	return 0;
#line 227 "MB01VD.f"
    }

#line 229 "MB01VD.f"
    dum[0] = 0.;
#line 230 "MB01VD.f"
    jc = 1;
#line 231 "MB01VD.f"
    nz = 0;

/*     Compute the Kronecker product of the matrices op(A) and op(B), */
/*        C = alpha*kron( op(A), op(B) ) + beta*C. */
/*     First, check if A is sparse. Here, A is considered as being sparse */
/*     if (number of zeros in A)/(MA*NA) >= SPARST. */

#line 238 "MB01VD.f"
    i__1 = *na;
#line 238 "MB01VD.f"
    for (j = 1; j <= i__1; ++j) {

#line 240 "MB01VD.f"
	i__2 = *ma;
#line 240 "MB01VD.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 241 "MB01VD.f"
	    if (transa) {
#line 242 "MB01VD.f"
		if (a[j + i__ * a_dim1] == 0.) {
#line 242 "MB01VD.f"
		    ++nz;
#line 242 "MB01VD.f"
		}
#line 244 "MB01VD.f"
	    } else {
#line 245 "MB01VD.f"
		if (a[i__ + j * a_dim1] == 0.) {
#line 245 "MB01VD.f"
		    ++nz;
#line 245 "MB01VD.f"
		}
#line 247 "MB01VD.f"
	    }
#line 248 "MB01VD.f"
/* L20: */
#line 248 "MB01VD.f"
	}

#line 250 "MB01VD.f"
/* L30: */
#line 250 "MB01VD.f"
    }

#line 252 "MB01VD.f"
    sparse = (doublereal) nz / (doublereal) (*ma * *na) >= .8;

#line 254 "MB01VD.f"
    if (! transa && ! transb) {

/*        Case op(A) = A and op(B) = B. */

#line 258 "MB01VD.f"
	if (*beta == 0.) {
#line 259 "MB01VD.f"
	    if (*alpha == 1.) {
#line 260 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 0, alpha = 1, A sparse. */

#line 264 "MB01VD.f"
		    i__1 = *na;
#line 264 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 266 "MB01VD.f"
			i__2 = *nb;
#line 266 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 267 "MB01VD.f"
			    ic = 1;

#line 269 "MB01VD.f"
			    i__3 = *ma;
#line 269 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 270 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 271 "MB01VD.f"
				if (aij == 0.) {
#line 272 "MB01VD.f"
				    dcopy_(mb, dum, &c__0, &c__[ic + jc * 
					    c_dim1], &c__1);
#line 273 "MB01VD.f"
				} else if (aij == 1.) {
#line 274 "MB01VD.f"
				    dcopy_(mb, &b[k * b_dim1 + 1], &c__1, &
					    c__[ic + jc * c_dim1], &c__1);
#line 275 "MB01VD.f"
				} else {
#line 276 "MB01VD.f"
				    lc = ic;

#line 278 "MB01VD.f"
				    i__4 = *mb;
#line 278 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 279 "MB01VD.f"
					c__[lc + jc * c_dim1] = aij * b[l + k 
						* b_dim1];
#line 280 "MB01VD.f"
					++lc;
#line 281 "MB01VD.f"
/* L50: */
#line 281 "MB01VD.f"
				    }

#line 283 "MB01VD.f"
				}
#line 284 "MB01VD.f"
				ic += *mb;
#line 285 "MB01VD.f"
/* L60: */
#line 285 "MB01VD.f"
			    }

#line 287 "MB01VD.f"
			    ++jc;
#line 288 "MB01VD.f"
/* L70: */
#line 288 "MB01VD.f"
			}

#line 290 "MB01VD.f"
/* L80: */
#line 290 "MB01VD.f"
		    }

#line 292 "MB01VD.f"
		} else {

/*                 Case beta = 0, alpha = 1, A not sparse. */

#line 296 "MB01VD.f"
		    i__1 = *na;
#line 296 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 298 "MB01VD.f"
			i__2 = *nb;
#line 298 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 299 "MB01VD.f"
			    ic = 1;

#line 301 "MB01VD.f"
			    i__3 = *ma;
#line 301 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 302 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 303 "MB01VD.f"
				lc = ic;

#line 305 "MB01VD.f"
				i__4 = *mb;
#line 305 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 306 "MB01VD.f"
				    c__[lc + jc * c_dim1] = aij * b[l + k * 
					    b_dim1];
#line 307 "MB01VD.f"
				    ++lc;
#line 308 "MB01VD.f"
/* L90: */
#line 308 "MB01VD.f"
				}

#line 310 "MB01VD.f"
				ic += *mb;
#line 311 "MB01VD.f"
/* L100: */
#line 311 "MB01VD.f"
			    }

#line 313 "MB01VD.f"
			    ++jc;
#line 314 "MB01VD.f"
/* L110: */
#line 314 "MB01VD.f"
			}

#line 316 "MB01VD.f"
/* L120: */
#line 316 "MB01VD.f"
		    }

#line 318 "MB01VD.f"
		}
#line 319 "MB01VD.f"
	    } else {
#line 320 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 0, alpha <> 1, A sparse. */

#line 324 "MB01VD.f"
		    i__1 = *na;
#line 324 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 326 "MB01VD.f"
			i__2 = *nb;
#line 326 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 327 "MB01VD.f"
			    ic = 1;

#line 329 "MB01VD.f"
			    i__3 = *ma;
#line 329 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 330 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 331 "MB01VD.f"
				if (aij == 0.) {
#line 332 "MB01VD.f"
				    dcopy_(mb, dum, &c__0, &c__[ic + jc * 
					    c_dim1], &c__1);
#line 333 "MB01VD.f"
				} else {
#line 334 "MB01VD.f"
				    lc = ic;

#line 336 "MB01VD.f"
				    i__4 = *mb;
#line 336 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 337 "MB01VD.f"
					c__[lc + jc * c_dim1] = aij * b[l + k 
						* b_dim1];
#line 338 "MB01VD.f"
					++lc;
#line 339 "MB01VD.f"
/* L130: */
#line 339 "MB01VD.f"
				    }

#line 341 "MB01VD.f"
				}
#line 342 "MB01VD.f"
				ic += *mb;
#line 343 "MB01VD.f"
/* L140: */
#line 343 "MB01VD.f"
			    }

#line 345 "MB01VD.f"
			    ++jc;
#line 346 "MB01VD.f"
/* L150: */
#line 346 "MB01VD.f"
			}

#line 348 "MB01VD.f"
/* L160: */
#line 348 "MB01VD.f"
		    }

#line 350 "MB01VD.f"
		} else {

/*                 Case beta = 0, alpha <> 1, A not sparse. */

#line 354 "MB01VD.f"
		    i__1 = *na;
#line 354 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 356 "MB01VD.f"
			i__2 = *nb;
#line 356 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 357 "MB01VD.f"
			    ic = 1;

#line 359 "MB01VD.f"
			    i__3 = *ma;
#line 359 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 360 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 361 "MB01VD.f"
				lc = ic;

#line 363 "MB01VD.f"
				i__4 = *mb;
#line 363 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 364 "MB01VD.f"
				    c__[lc + jc * c_dim1] = aij * b[l + k * 
					    b_dim1];
#line 365 "MB01VD.f"
				    ++lc;
#line 366 "MB01VD.f"
/* L170: */
#line 366 "MB01VD.f"
				}

#line 368 "MB01VD.f"
				ic += *mb;
#line 369 "MB01VD.f"
/* L180: */
#line 369 "MB01VD.f"
			    }

#line 371 "MB01VD.f"
			    ++jc;
#line 372 "MB01VD.f"
/* L190: */
#line 372 "MB01VD.f"
			}

#line 374 "MB01VD.f"
/* L200: */
#line 374 "MB01VD.f"
		    }

#line 376 "MB01VD.f"
		}
#line 377 "MB01VD.f"
	    }
#line 378 "MB01VD.f"
	} else if (*beta == 1.) {
#line 379 "MB01VD.f"
	    if (*alpha == 1.) {
#line 380 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 1, alpha = 1, A sparse. */

#line 384 "MB01VD.f"
		    i__1 = *na;
#line 384 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 386 "MB01VD.f"
			i__2 = *nb;
#line 386 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 387 "MB01VD.f"
			    ic = 1;

#line 389 "MB01VD.f"
			    i__3 = *ma;
#line 389 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 390 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 391 "MB01VD.f"
				if (aij != 0.) {
#line 392 "MB01VD.f"
				    lc = ic;

#line 394 "MB01VD.f"
				    i__4 = *mb;
#line 394 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 395 "MB01VD.f"
					c__[lc + jc * c_dim1] += aij * b[l + 
						k * b_dim1];
#line 396 "MB01VD.f"
					++lc;
#line 397 "MB01VD.f"
/* L210: */
#line 397 "MB01VD.f"
				    }

#line 399 "MB01VD.f"
				}
#line 400 "MB01VD.f"
				ic += *mb;
#line 401 "MB01VD.f"
/* L220: */
#line 401 "MB01VD.f"
			    }

#line 403 "MB01VD.f"
			    ++jc;
#line 404 "MB01VD.f"
/* L230: */
#line 404 "MB01VD.f"
			}

#line 406 "MB01VD.f"
/* L240: */
#line 406 "MB01VD.f"
		    }

#line 408 "MB01VD.f"
		} else {

/*                 Case beta = 1, alpha = 1, A not sparse. */

#line 412 "MB01VD.f"
		    i__1 = *na;
#line 412 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 414 "MB01VD.f"
			i__2 = *nb;
#line 414 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 415 "MB01VD.f"
			    ic = 1;

#line 417 "MB01VD.f"
			    i__3 = *ma;
#line 417 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 418 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 419 "MB01VD.f"
				lc = ic;

#line 421 "MB01VD.f"
				i__4 = *mb;
#line 421 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 422 "MB01VD.f"
				    c__[lc + jc * c_dim1] += aij * b[l + k * 
					    b_dim1];
#line 423 "MB01VD.f"
				    ++lc;
#line 424 "MB01VD.f"
/* L250: */
#line 424 "MB01VD.f"
				}

#line 426 "MB01VD.f"
				ic += *mb;
#line 427 "MB01VD.f"
/* L260: */
#line 427 "MB01VD.f"
			    }

#line 429 "MB01VD.f"
			    ++jc;
#line 430 "MB01VD.f"
/* L270: */
#line 430 "MB01VD.f"
			}

#line 432 "MB01VD.f"
/* L280: */
#line 432 "MB01VD.f"
		    }

#line 434 "MB01VD.f"
		}
#line 435 "MB01VD.f"
	    } else {
#line 436 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 1, alpha <> 1, A sparse. */

#line 440 "MB01VD.f"
		    i__1 = *na;
#line 440 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 442 "MB01VD.f"
			i__2 = *nb;
#line 442 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 443 "MB01VD.f"
			    ic = 1;

#line 445 "MB01VD.f"
			    i__3 = *ma;
#line 445 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 446 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 447 "MB01VD.f"
				if (aij != 0.) {
#line 448 "MB01VD.f"
				    lc = ic;

#line 450 "MB01VD.f"
				    i__4 = *mb;
#line 450 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 451 "MB01VD.f"
					c__[lc + jc * c_dim1] += aij * b[l + 
						k * b_dim1];
#line 452 "MB01VD.f"
					++lc;
#line 453 "MB01VD.f"
/* L290: */
#line 453 "MB01VD.f"
				    }

#line 455 "MB01VD.f"
				}
#line 456 "MB01VD.f"
				ic += *mb;
#line 457 "MB01VD.f"
/* L300: */
#line 457 "MB01VD.f"
			    }

#line 459 "MB01VD.f"
			    ++jc;
#line 460 "MB01VD.f"
/* L310: */
#line 460 "MB01VD.f"
			}

#line 462 "MB01VD.f"
/* L320: */
#line 462 "MB01VD.f"
		    }

#line 464 "MB01VD.f"
		} else {

/*                 Case beta = 1, alpha <> 1, A not sparse. */

#line 468 "MB01VD.f"
		    i__1 = *na;
#line 468 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 470 "MB01VD.f"
			i__2 = *nb;
#line 470 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 471 "MB01VD.f"
			    ic = 1;

#line 473 "MB01VD.f"
			    i__3 = *ma;
#line 473 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 474 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 475 "MB01VD.f"
				lc = ic;

#line 477 "MB01VD.f"
				i__4 = *mb;
#line 477 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 478 "MB01VD.f"
				    c__[lc + jc * c_dim1] += aij * b[l + k * 
					    b_dim1];
#line 479 "MB01VD.f"
				    ++lc;
#line 480 "MB01VD.f"
/* L330: */
#line 480 "MB01VD.f"
				}

#line 482 "MB01VD.f"
				ic += *mb;
#line 483 "MB01VD.f"
/* L340: */
#line 483 "MB01VD.f"
			    }

#line 485 "MB01VD.f"
			    ++jc;
#line 486 "MB01VD.f"
/* L350: */
#line 486 "MB01VD.f"
			}

#line 488 "MB01VD.f"
/* L360: */
#line 488 "MB01VD.f"
		    }

#line 490 "MB01VD.f"
		}
#line 491 "MB01VD.f"
	    }
#line 492 "MB01VD.f"
	} else {
#line 493 "MB01VD.f"
	    if (*alpha == 1.) {
#line 494 "MB01VD.f"
		if (sparse) {

/*                 Case beta <> 0 or 1, alpha = 1, A sparse. */

#line 498 "MB01VD.f"
		    i__1 = *na;
#line 498 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 500 "MB01VD.f"
			i__2 = *nb;
#line 500 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 501 "MB01VD.f"
			    ic = 1;

#line 503 "MB01VD.f"
			    i__3 = *ma;
#line 503 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 504 "MB01VD.f"
				aij = a[i__ + j * a_dim1];

#line 506 "MB01VD.f"
				if (aij == 0.) {
#line 507 "MB01VD.f"
				    dscal_(mb, beta, &c__[ic + jc * c_dim1], &
					    c__1);
#line 508 "MB01VD.f"
				} else {
#line 509 "MB01VD.f"
				    lc = ic;

#line 511 "MB01VD.f"
				    i__4 = *mb;
#line 511 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 512 "MB01VD.f"
					c__[lc + jc * c_dim1] = *beta * c__[
						lc + jc * c_dim1] + aij * b[l 
						+ k * b_dim1];
#line 513 "MB01VD.f"
					++lc;
#line 514 "MB01VD.f"
/* L370: */
#line 514 "MB01VD.f"
				    }

#line 516 "MB01VD.f"
				}
#line 517 "MB01VD.f"
				ic += *mb;
#line 518 "MB01VD.f"
/* L380: */
#line 518 "MB01VD.f"
			    }

#line 520 "MB01VD.f"
			    ++jc;
#line 521 "MB01VD.f"
/* L390: */
#line 521 "MB01VD.f"
			}

#line 523 "MB01VD.f"
/* L400: */
#line 523 "MB01VD.f"
		    }

#line 525 "MB01VD.f"
		} else {

/*                 Case beta <> 0 or 1, alpha = 1, A not sparse. */

#line 529 "MB01VD.f"
		    i__1 = *na;
#line 529 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 531 "MB01VD.f"
			i__2 = *nb;
#line 531 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 532 "MB01VD.f"
			    ic = 1;

#line 534 "MB01VD.f"
			    i__3 = *ma;
#line 534 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 535 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 536 "MB01VD.f"
				lc = ic;

#line 538 "MB01VD.f"
				i__4 = *mb;
#line 538 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 539 "MB01VD.f"
				    c__[lc + jc * c_dim1] = *beta * c__[lc + 
					    jc * c_dim1] + aij * b[l + k * 
					    b_dim1];
#line 540 "MB01VD.f"
				    ++lc;
#line 541 "MB01VD.f"
/* L410: */
#line 541 "MB01VD.f"
				}

#line 543 "MB01VD.f"
				ic += *mb;
#line 544 "MB01VD.f"
/* L420: */
#line 544 "MB01VD.f"
			    }

#line 546 "MB01VD.f"
			    ++jc;
#line 547 "MB01VD.f"
/* L430: */
#line 547 "MB01VD.f"
			}

#line 549 "MB01VD.f"
/* L440: */
#line 549 "MB01VD.f"
		    }

#line 551 "MB01VD.f"
		}
#line 552 "MB01VD.f"
	    } else {
#line 553 "MB01VD.f"
		if (sparse) {

/*                 Case beta <> 0 or 1, alpha <> 1, A sparse. */

#line 557 "MB01VD.f"
		    i__1 = *na;
#line 557 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 559 "MB01VD.f"
			i__2 = *nb;
#line 559 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 560 "MB01VD.f"
			    ic = 1;

#line 562 "MB01VD.f"
			    i__3 = *ma;
#line 562 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 563 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];

#line 565 "MB01VD.f"
				if (aij == 0.) {
#line 566 "MB01VD.f"
				    dscal_(mb, beta, &c__[ic + jc * c_dim1], &
					    c__1);
#line 567 "MB01VD.f"
				} else {
#line 568 "MB01VD.f"
				    lc = ic;

#line 570 "MB01VD.f"
				    i__4 = *mb;
#line 570 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 571 "MB01VD.f"
					c__[lc + jc * c_dim1] = *beta * c__[
						lc + jc * c_dim1] + aij * b[l 
						+ k * b_dim1];
#line 572 "MB01VD.f"
					++lc;
#line 573 "MB01VD.f"
/* L450: */
#line 573 "MB01VD.f"
				    }

#line 575 "MB01VD.f"
				}
#line 576 "MB01VD.f"
				ic += *mb;
#line 577 "MB01VD.f"
/* L460: */
#line 577 "MB01VD.f"
			    }

#line 579 "MB01VD.f"
			    ++jc;
#line 580 "MB01VD.f"
/* L470: */
#line 580 "MB01VD.f"
			}

#line 582 "MB01VD.f"
/* L480: */
#line 582 "MB01VD.f"
		    }

#line 584 "MB01VD.f"
		} else {

/*                 Case beta <> 0 or 1, alpha <> 1, A not sparse. */

#line 588 "MB01VD.f"
		    i__1 = *na;
#line 588 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 590 "MB01VD.f"
			i__2 = *nb;
#line 590 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 591 "MB01VD.f"
			    ic = 1;

#line 593 "MB01VD.f"
			    i__3 = *ma;
#line 593 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 594 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 595 "MB01VD.f"
				lc = ic;

#line 597 "MB01VD.f"
				i__4 = *mb;
#line 597 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 598 "MB01VD.f"
				    c__[lc + jc * c_dim1] = *beta * c__[lc + 
					    jc * c_dim1] + aij * b[l + k * 
					    b_dim1];
#line 599 "MB01VD.f"
				    ++lc;
#line 600 "MB01VD.f"
/* L490: */
#line 600 "MB01VD.f"
				}

#line 602 "MB01VD.f"
				ic += *mb;
#line 603 "MB01VD.f"
/* L500: */
#line 603 "MB01VD.f"
			    }

#line 605 "MB01VD.f"
			    ++jc;
#line 606 "MB01VD.f"
/* L510: */
#line 606 "MB01VD.f"
			}

#line 608 "MB01VD.f"
/* L520: */
#line 608 "MB01VD.f"
		    }

#line 610 "MB01VD.f"
		}
#line 611 "MB01VD.f"
	    }
#line 612 "MB01VD.f"
	}
#line 613 "MB01VD.f"
    } else if (transa && ! transb) {

/*        Case op(A) = A' and op(B) = B. */

#line 617 "MB01VD.f"
	if (*beta == 0.) {
#line 618 "MB01VD.f"
	    if (*alpha == 1.) {
#line 619 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 0, alpha = 1, A sparse. */

#line 623 "MB01VD.f"
		    i__1 = *na;
#line 623 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 625 "MB01VD.f"
			i__2 = *nb;
#line 625 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 626 "MB01VD.f"
			    ic = 1;

#line 628 "MB01VD.f"
			    i__3 = *ma;
#line 628 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 629 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 630 "MB01VD.f"
				if (aij == 0.) {
#line 631 "MB01VD.f"
				    dcopy_(mb, dum, &c__0, &c__[ic + jc * 
					    c_dim1], &c__1);
#line 632 "MB01VD.f"
				} else if (aij == 1.) {
#line 633 "MB01VD.f"
				    dcopy_(mb, &b[k * b_dim1 + 1], &c__1, &
					    c__[ic + jc * c_dim1], &c__1);
#line 634 "MB01VD.f"
				} else {
#line 635 "MB01VD.f"
				    lc = ic;

#line 637 "MB01VD.f"
				    i__4 = *mb;
#line 637 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 638 "MB01VD.f"
					c__[lc + jc * c_dim1] = aij * b[l + k 
						* b_dim1];
#line 639 "MB01VD.f"
					++lc;
#line 640 "MB01VD.f"
/* L530: */
#line 640 "MB01VD.f"
				    }

#line 642 "MB01VD.f"
				}
#line 643 "MB01VD.f"
				ic += *mb;
#line 644 "MB01VD.f"
/* L540: */
#line 644 "MB01VD.f"
			    }

#line 646 "MB01VD.f"
			    ++jc;
#line 647 "MB01VD.f"
/* L550: */
#line 647 "MB01VD.f"
			}

#line 649 "MB01VD.f"
/* L560: */
#line 649 "MB01VD.f"
		    }

#line 651 "MB01VD.f"
		} else {

/*                 Case beta = 0, alpha = 1, A not sparse. */

#line 655 "MB01VD.f"
		    i__1 = *na;
#line 655 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 657 "MB01VD.f"
			i__2 = *nb;
#line 657 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 658 "MB01VD.f"
			    ic = 1;

#line 660 "MB01VD.f"
			    i__3 = *ma;
#line 660 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 661 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 662 "MB01VD.f"
				lc = ic;

#line 664 "MB01VD.f"
				i__4 = *mb;
#line 664 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 665 "MB01VD.f"
				    c__[lc + jc * c_dim1] = aij * b[l + k * 
					    b_dim1];
#line 666 "MB01VD.f"
				    ++lc;
#line 667 "MB01VD.f"
/* L570: */
#line 667 "MB01VD.f"
				}

#line 669 "MB01VD.f"
				ic += *mb;
#line 670 "MB01VD.f"
/* L580: */
#line 670 "MB01VD.f"
			    }

#line 672 "MB01VD.f"
			    ++jc;
#line 673 "MB01VD.f"
/* L590: */
#line 673 "MB01VD.f"
			}

#line 675 "MB01VD.f"
/* L600: */
#line 675 "MB01VD.f"
		    }

#line 677 "MB01VD.f"
		}
#line 678 "MB01VD.f"
	    } else {
#line 679 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 0, alpha <> 1, A sparse. */

#line 683 "MB01VD.f"
		    i__1 = *na;
#line 683 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 685 "MB01VD.f"
			i__2 = *nb;
#line 685 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 686 "MB01VD.f"
			    ic = 1;

#line 688 "MB01VD.f"
			    i__3 = *ma;
#line 688 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 689 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 690 "MB01VD.f"
				if (aij == 0.) {
#line 691 "MB01VD.f"
				    dcopy_(mb, dum, &c__0, &c__[ic + jc * 
					    c_dim1], &c__1);
#line 692 "MB01VD.f"
				} else {
#line 693 "MB01VD.f"
				    lc = ic;

#line 695 "MB01VD.f"
				    i__4 = *mb;
#line 695 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 696 "MB01VD.f"
					c__[lc + jc * c_dim1] = aij * b[l + k 
						* b_dim1];
#line 697 "MB01VD.f"
					++lc;
#line 698 "MB01VD.f"
/* L610: */
#line 698 "MB01VD.f"
				    }

#line 700 "MB01VD.f"
				}
#line 701 "MB01VD.f"
				ic += *mb;
#line 702 "MB01VD.f"
/* L620: */
#line 702 "MB01VD.f"
			    }

#line 704 "MB01VD.f"
			    ++jc;
#line 705 "MB01VD.f"
/* L630: */
#line 705 "MB01VD.f"
			}

#line 707 "MB01VD.f"
/* L640: */
#line 707 "MB01VD.f"
		    }

#line 709 "MB01VD.f"
		} else {

/*                 Case beta = 0, alpha <> 1, A not sparse. */

#line 713 "MB01VD.f"
		    i__1 = *na;
#line 713 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 715 "MB01VD.f"
			i__2 = *nb;
#line 715 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 716 "MB01VD.f"
			    ic = 1;

#line 718 "MB01VD.f"
			    i__3 = *ma;
#line 718 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 719 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 720 "MB01VD.f"
				lc = ic;

#line 722 "MB01VD.f"
				i__4 = *mb;
#line 722 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 723 "MB01VD.f"
				    c__[lc + jc * c_dim1] = aij * b[l + k * 
					    b_dim1];
#line 724 "MB01VD.f"
				    ++lc;
#line 725 "MB01VD.f"
/* L650: */
#line 725 "MB01VD.f"
				}

#line 727 "MB01VD.f"
				ic += *mb;
#line 728 "MB01VD.f"
/* L660: */
#line 728 "MB01VD.f"
			    }

#line 730 "MB01VD.f"
			    ++jc;
#line 731 "MB01VD.f"
/* L670: */
#line 731 "MB01VD.f"
			}

#line 733 "MB01VD.f"
/* L680: */
#line 733 "MB01VD.f"
		    }

#line 735 "MB01VD.f"
		}
#line 736 "MB01VD.f"
	    }
#line 737 "MB01VD.f"
	} else if (*beta == 1.) {
#line 738 "MB01VD.f"
	    if (*alpha == 1.) {
#line 739 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 1, alpha = 1, A sparse. */

#line 743 "MB01VD.f"
		    i__1 = *na;
#line 743 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 745 "MB01VD.f"
			i__2 = *nb;
#line 745 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 746 "MB01VD.f"
			    ic = 1;

#line 748 "MB01VD.f"
			    i__3 = *ma;
#line 748 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 749 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 750 "MB01VD.f"
				if (aij != 0.) {
#line 751 "MB01VD.f"
				    lc = ic;

#line 753 "MB01VD.f"
				    i__4 = *mb;
#line 753 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 754 "MB01VD.f"
					c__[lc + jc * c_dim1] += aij * b[l + 
						k * b_dim1];
#line 755 "MB01VD.f"
					++lc;
#line 756 "MB01VD.f"
/* L690: */
#line 756 "MB01VD.f"
				    }

#line 758 "MB01VD.f"
				}
#line 759 "MB01VD.f"
				ic += *mb;
#line 760 "MB01VD.f"
/* L700: */
#line 760 "MB01VD.f"
			    }

#line 762 "MB01VD.f"
			    ++jc;
#line 763 "MB01VD.f"
/* L710: */
#line 763 "MB01VD.f"
			}

#line 765 "MB01VD.f"
/* L720: */
#line 765 "MB01VD.f"
		    }

#line 767 "MB01VD.f"
		} else {

/*                 Case beta = 1, alpha = 1, A not sparse. */

#line 771 "MB01VD.f"
		    i__1 = *na;
#line 771 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 773 "MB01VD.f"
			i__2 = *nb;
#line 773 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 774 "MB01VD.f"
			    ic = 1;

#line 776 "MB01VD.f"
			    i__3 = *ma;
#line 776 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 777 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 778 "MB01VD.f"
				lc = ic;

#line 780 "MB01VD.f"
				i__4 = *mb;
#line 780 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 781 "MB01VD.f"
				    c__[lc + jc * c_dim1] += aij * b[l + k * 
					    b_dim1];
#line 782 "MB01VD.f"
				    ++lc;
#line 783 "MB01VD.f"
/* L730: */
#line 783 "MB01VD.f"
				}

#line 785 "MB01VD.f"
				ic += *mb;
#line 786 "MB01VD.f"
/* L740: */
#line 786 "MB01VD.f"
			    }

#line 788 "MB01VD.f"
			    ++jc;
#line 789 "MB01VD.f"
/* L750: */
#line 789 "MB01VD.f"
			}

#line 791 "MB01VD.f"
/* L760: */
#line 791 "MB01VD.f"
		    }

#line 793 "MB01VD.f"
		}
#line 794 "MB01VD.f"
	    } else {
#line 795 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 1, alpha <> 1, A sparse. */

#line 799 "MB01VD.f"
		    i__1 = *na;
#line 799 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 801 "MB01VD.f"
			i__2 = *nb;
#line 801 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 802 "MB01VD.f"
			    ic = 1;

#line 804 "MB01VD.f"
			    i__3 = *ma;
#line 804 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 805 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 806 "MB01VD.f"
				if (aij != 0.) {
#line 807 "MB01VD.f"
				    lc = ic;

#line 809 "MB01VD.f"
				    i__4 = *mb;
#line 809 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 810 "MB01VD.f"
					c__[lc + jc * c_dim1] += aij * b[l + 
						k * b_dim1];
#line 811 "MB01VD.f"
					++lc;
#line 812 "MB01VD.f"
/* L770: */
#line 812 "MB01VD.f"
				    }

#line 814 "MB01VD.f"
				}
#line 815 "MB01VD.f"
				ic += *mb;
#line 816 "MB01VD.f"
/* L780: */
#line 816 "MB01VD.f"
			    }

#line 818 "MB01VD.f"
			    ++jc;
#line 819 "MB01VD.f"
/* L790: */
#line 819 "MB01VD.f"
			}

#line 821 "MB01VD.f"
/* L800: */
#line 821 "MB01VD.f"
		    }

#line 823 "MB01VD.f"
		} else {

/*                 Case beta = 1, alpha <> 1, A not sparse. */

#line 827 "MB01VD.f"
		    i__1 = *na;
#line 827 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 829 "MB01VD.f"
			i__2 = *nb;
#line 829 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 830 "MB01VD.f"
			    ic = 1;

#line 832 "MB01VD.f"
			    i__3 = *ma;
#line 832 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 833 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 834 "MB01VD.f"
				lc = ic;

#line 836 "MB01VD.f"
				i__4 = *mb;
#line 836 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 837 "MB01VD.f"
				    c__[lc + jc * c_dim1] += aij * b[l + k * 
					    b_dim1];
#line 838 "MB01VD.f"
				    ++lc;
#line 839 "MB01VD.f"
/* L810: */
#line 839 "MB01VD.f"
				}

#line 841 "MB01VD.f"
				ic += *mb;
#line 842 "MB01VD.f"
/* L820: */
#line 842 "MB01VD.f"
			    }

#line 844 "MB01VD.f"
			    ++jc;
#line 845 "MB01VD.f"
/* L830: */
#line 845 "MB01VD.f"
			}

#line 847 "MB01VD.f"
/* L840: */
#line 847 "MB01VD.f"
		    }

#line 849 "MB01VD.f"
		}
#line 850 "MB01VD.f"
	    }
#line 851 "MB01VD.f"
	} else {
#line 852 "MB01VD.f"
	    if (*alpha == 1.) {
#line 853 "MB01VD.f"
		if (sparse) {

/*                 Case beta <> 0 or 1, alpha = 1, A sparse. */

#line 857 "MB01VD.f"
		    i__1 = *na;
#line 857 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 859 "MB01VD.f"
			i__2 = *nb;
#line 859 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 860 "MB01VD.f"
			    ic = 1;

#line 862 "MB01VD.f"
			    i__3 = *ma;
#line 862 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 863 "MB01VD.f"
				aij = a[j + i__ * a_dim1];

#line 865 "MB01VD.f"
				if (aij == 0.) {
#line 866 "MB01VD.f"
				    dscal_(mb, beta, &c__[ic + jc * c_dim1], &
					    c__1);
#line 867 "MB01VD.f"
				} else {
#line 868 "MB01VD.f"
				    lc = ic;

#line 870 "MB01VD.f"
				    i__4 = *mb;
#line 870 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 871 "MB01VD.f"
					c__[lc + jc * c_dim1] = *beta * c__[
						lc + jc * c_dim1] + aij * b[l 
						+ k * b_dim1];
#line 872 "MB01VD.f"
					++lc;
#line 873 "MB01VD.f"
/* L850: */
#line 873 "MB01VD.f"
				    }

#line 875 "MB01VD.f"
				}
#line 876 "MB01VD.f"
				ic += *mb;
#line 877 "MB01VD.f"
/* L860: */
#line 877 "MB01VD.f"
			    }

#line 879 "MB01VD.f"
			    ++jc;
#line 880 "MB01VD.f"
/* L870: */
#line 880 "MB01VD.f"
			}

#line 882 "MB01VD.f"
/* L880: */
#line 882 "MB01VD.f"
		    }

#line 884 "MB01VD.f"
		} else {

/*                 Case beta <> 0 or 1, alpha = 1, A not sparse. */

#line 888 "MB01VD.f"
		    i__1 = *na;
#line 888 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 890 "MB01VD.f"
			i__2 = *nb;
#line 890 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 891 "MB01VD.f"
			    ic = 1;

#line 893 "MB01VD.f"
			    i__3 = *ma;
#line 893 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 894 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 895 "MB01VD.f"
				lc = ic;

#line 897 "MB01VD.f"
				i__4 = *mb;
#line 897 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 898 "MB01VD.f"
				    c__[lc + jc * c_dim1] = *beta * c__[lc + 
					    jc * c_dim1] + aij * b[l + k * 
					    b_dim1];
#line 899 "MB01VD.f"
				    ++lc;
#line 900 "MB01VD.f"
/* L890: */
#line 900 "MB01VD.f"
				}

#line 902 "MB01VD.f"
				ic += *mb;
#line 903 "MB01VD.f"
/* L900: */
#line 903 "MB01VD.f"
			    }

#line 905 "MB01VD.f"
			    ++jc;
#line 906 "MB01VD.f"
/* L910: */
#line 906 "MB01VD.f"
			}

#line 908 "MB01VD.f"
/* L920: */
#line 908 "MB01VD.f"
		    }

#line 910 "MB01VD.f"
		}
#line 911 "MB01VD.f"
	    } else {
#line 912 "MB01VD.f"
		if (sparse) {

/*                 Case beta <> 0 or 1, alpha <> 1, A sparse. */

#line 916 "MB01VD.f"
		    i__1 = *na;
#line 916 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 918 "MB01VD.f"
			i__2 = *nb;
#line 918 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 919 "MB01VD.f"
			    ic = 1;

#line 921 "MB01VD.f"
			    i__3 = *ma;
#line 921 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 922 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];

#line 924 "MB01VD.f"
				if (aij == 0.) {
#line 925 "MB01VD.f"
				    dscal_(mb, beta, &c__[ic + jc * c_dim1], &
					    c__1);
#line 926 "MB01VD.f"
				} else {
#line 927 "MB01VD.f"
				    lc = ic;

#line 929 "MB01VD.f"
				    i__4 = *mb;
#line 929 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 930 "MB01VD.f"
					c__[lc + jc * c_dim1] = *beta * c__[
						lc + jc * c_dim1] + aij * b[l 
						+ k * b_dim1];
#line 931 "MB01VD.f"
					++lc;
#line 932 "MB01VD.f"
/* L930: */
#line 932 "MB01VD.f"
				    }

#line 934 "MB01VD.f"
				}
#line 935 "MB01VD.f"
				ic += *mb;
#line 936 "MB01VD.f"
/* L940: */
#line 936 "MB01VD.f"
			    }

#line 938 "MB01VD.f"
			    ++jc;
#line 939 "MB01VD.f"
/* L950: */
#line 939 "MB01VD.f"
			}

#line 941 "MB01VD.f"
/* L960: */
#line 941 "MB01VD.f"
		    }

#line 943 "MB01VD.f"
		} else {

/*                 Case beta <> 0 or 1, alpha <> 1, A not sparse. */

#line 947 "MB01VD.f"
		    i__1 = *na;
#line 947 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 949 "MB01VD.f"
			i__2 = *nb;
#line 949 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 950 "MB01VD.f"
			    ic = 1;

#line 952 "MB01VD.f"
			    i__3 = *ma;
#line 952 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 953 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 954 "MB01VD.f"
				lc = ic;

#line 956 "MB01VD.f"
				i__4 = *mb;
#line 956 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 957 "MB01VD.f"
				    c__[lc + jc * c_dim1] = *beta * c__[lc + 
					    jc * c_dim1] + aij * b[l + k * 
					    b_dim1];
#line 958 "MB01VD.f"
				    ++lc;
#line 959 "MB01VD.f"
/* L970: */
#line 959 "MB01VD.f"
				}

#line 961 "MB01VD.f"
				ic += *mb;
#line 962 "MB01VD.f"
/* L980: */
#line 962 "MB01VD.f"
			    }

#line 964 "MB01VD.f"
			    ++jc;
#line 965 "MB01VD.f"
/* L990: */
#line 965 "MB01VD.f"
			}

#line 967 "MB01VD.f"
/* L1000: */
#line 967 "MB01VD.f"
		    }

#line 969 "MB01VD.f"
		}
#line 970 "MB01VD.f"
	    }
#line 971 "MB01VD.f"
	}
#line 972 "MB01VD.f"
    } else if (transb && ! transa) {

/*        Case op(A) = A and op(B) = B'. */

#line 976 "MB01VD.f"
	if (*beta == 0.) {
#line 977 "MB01VD.f"
	    if (*alpha == 1.) {
#line 978 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 0, alpha = 1, A sparse. */

#line 982 "MB01VD.f"
		    i__1 = *na;
#line 982 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 984 "MB01VD.f"
			i__2 = *nb;
#line 984 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 985 "MB01VD.f"
			    ic = 1;

#line 987 "MB01VD.f"
			    i__3 = *ma;
#line 987 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 988 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 989 "MB01VD.f"
				if (aij == 0.) {
#line 990 "MB01VD.f"
				    dcopy_(mb, dum, &c__0, &c__[ic + jc * 
					    c_dim1], &c__1);
#line 991 "MB01VD.f"
				} else if (aij == 1.) {
#line 992 "MB01VD.f"
				    dcopy_(mb, &b[k + b_dim1], ldb, &c__[ic + 
					    jc * c_dim1], &c__1);
#line 993 "MB01VD.f"
				} else {
#line 994 "MB01VD.f"
				    lc = ic;

#line 996 "MB01VD.f"
				    i__4 = *mb;
#line 996 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 997 "MB01VD.f"
					c__[lc + jc * c_dim1] = aij * b[k + l 
						* b_dim1];
#line 998 "MB01VD.f"
					++lc;
#line 999 "MB01VD.f"
/* L1050: */
#line 999 "MB01VD.f"
				    }

#line 1001 "MB01VD.f"
				}
#line 1002 "MB01VD.f"
				ic += *mb;
#line 1003 "MB01VD.f"
/* L1060: */
#line 1003 "MB01VD.f"
			    }

#line 1005 "MB01VD.f"
			    ++jc;
#line 1006 "MB01VD.f"
/* L1070: */
#line 1006 "MB01VD.f"
			}

#line 1008 "MB01VD.f"
/* L1080: */
#line 1008 "MB01VD.f"
		    }

#line 1010 "MB01VD.f"
		} else {

/*                 Case beta = 0, alpha = 1, A not sparse. */

#line 1014 "MB01VD.f"
		    i__1 = *na;
#line 1014 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1016 "MB01VD.f"
			i__2 = *nb;
#line 1016 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1017 "MB01VD.f"
			    ic = 1;

#line 1019 "MB01VD.f"
			    i__3 = *ma;
#line 1019 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1020 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 1021 "MB01VD.f"
				lc = ic;

#line 1023 "MB01VD.f"
				i__4 = *mb;
#line 1023 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1024 "MB01VD.f"
				    c__[lc + jc * c_dim1] = aij * b[k + l * 
					    b_dim1];
#line 1025 "MB01VD.f"
				    ++lc;
#line 1026 "MB01VD.f"
/* L1090: */
#line 1026 "MB01VD.f"
				}

#line 1028 "MB01VD.f"
				ic += *mb;
#line 1029 "MB01VD.f"
/* L1100: */
#line 1029 "MB01VD.f"
			    }

#line 1031 "MB01VD.f"
			    ++jc;
#line 1032 "MB01VD.f"
/* L1110: */
#line 1032 "MB01VD.f"
			}

#line 1034 "MB01VD.f"
/* L1120: */
#line 1034 "MB01VD.f"
		    }

#line 1036 "MB01VD.f"
		}
#line 1037 "MB01VD.f"
	    } else {
#line 1038 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 0, alpha <> 1, A sparse. */

#line 1042 "MB01VD.f"
		    i__1 = *na;
#line 1042 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1044 "MB01VD.f"
			i__2 = *nb;
#line 1044 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1045 "MB01VD.f"
			    ic = 1;

#line 1047 "MB01VD.f"
			    i__3 = *ma;
#line 1047 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1048 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 1049 "MB01VD.f"
				if (aij == 0.) {
#line 1050 "MB01VD.f"
				    dcopy_(mb, dum, &c__0, &c__[ic + jc * 
					    c_dim1], &c__1);
#line 1051 "MB01VD.f"
				} else {
#line 1052 "MB01VD.f"
				    lc = ic;

#line 1054 "MB01VD.f"
				    i__4 = *mb;
#line 1054 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1055 "MB01VD.f"
					c__[lc + jc * c_dim1] = aij * b[k + l 
						* b_dim1];
#line 1056 "MB01VD.f"
					++lc;
#line 1057 "MB01VD.f"
/* L1130: */
#line 1057 "MB01VD.f"
				    }

#line 1059 "MB01VD.f"
				}
#line 1060 "MB01VD.f"
				ic += *mb;
#line 1061 "MB01VD.f"
/* L1140: */
#line 1061 "MB01VD.f"
			    }

#line 1063 "MB01VD.f"
			    ++jc;
#line 1064 "MB01VD.f"
/* L1150: */
#line 1064 "MB01VD.f"
			}

#line 1066 "MB01VD.f"
/* L1160: */
#line 1066 "MB01VD.f"
		    }

#line 1068 "MB01VD.f"
		} else {

/*                 Case beta = 0, alpha <> 1, A not sparse. */

#line 1072 "MB01VD.f"
		    i__1 = *na;
#line 1072 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1074 "MB01VD.f"
			i__2 = *nb;
#line 1074 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1075 "MB01VD.f"
			    ic = 1;

#line 1077 "MB01VD.f"
			    i__3 = *ma;
#line 1077 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1078 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 1079 "MB01VD.f"
				lc = ic;

#line 1081 "MB01VD.f"
				i__4 = *mb;
#line 1081 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1082 "MB01VD.f"
				    c__[lc + jc * c_dim1] = aij * b[k + l * 
					    b_dim1];
#line 1083 "MB01VD.f"
				    ++lc;
#line 1084 "MB01VD.f"
/* L1170: */
#line 1084 "MB01VD.f"
				}

#line 1086 "MB01VD.f"
				ic += *mb;
#line 1087 "MB01VD.f"
/* L1180: */
#line 1087 "MB01VD.f"
			    }

#line 1089 "MB01VD.f"
			    ++jc;
#line 1090 "MB01VD.f"
/* L1190: */
#line 1090 "MB01VD.f"
			}

#line 1092 "MB01VD.f"
/* L1200: */
#line 1092 "MB01VD.f"
		    }

#line 1094 "MB01VD.f"
		}
#line 1095 "MB01VD.f"
	    }
#line 1096 "MB01VD.f"
	} else if (*beta == 1.) {
#line 1097 "MB01VD.f"
	    if (*alpha == 1.) {
#line 1098 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 1, alpha = 1, A sparse. */

#line 1102 "MB01VD.f"
		    i__1 = *na;
#line 1102 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1104 "MB01VD.f"
			i__2 = *nb;
#line 1104 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1105 "MB01VD.f"
			    ic = 1;

#line 1107 "MB01VD.f"
			    i__3 = *ma;
#line 1107 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1108 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 1109 "MB01VD.f"
				if (aij != 0.) {
#line 1110 "MB01VD.f"
				    lc = ic;

#line 1112 "MB01VD.f"
				    i__4 = *mb;
#line 1112 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1113 "MB01VD.f"
					c__[lc + jc * c_dim1] += aij * b[k + 
						l * b_dim1];
#line 1114 "MB01VD.f"
					++lc;
#line 1115 "MB01VD.f"
/* L1210: */
#line 1115 "MB01VD.f"
				    }

#line 1117 "MB01VD.f"
				}
#line 1118 "MB01VD.f"
				ic += *mb;
#line 1119 "MB01VD.f"
/* L1220: */
#line 1119 "MB01VD.f"
			    }

#line 1121 "MB01VD.f"
			    ++jc;
#line 1122 "MB01VD.f"
/* L1230: */
#line 1122 "MB01VD.f"
			}

#line 1124 "MB01VD.f"
/* L1240: */
#line 1124 "MB01VD.f"
		    }

#line 1126 "MB01VD.f"
		} else {

/*                 Case beta = 1, alpha = 1, A not sparse. */

#line 1130 "MB01VD.f"
		    i__1 = *na;
#line 1130 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1132 "MB01VD.f"
			i__2 = *nb;
#line 1132 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1133 "MB01VD.f"
			    ic = 1;

#line 1135 "MB01VD.f"
			    i__3 = *ma;
#line 1135 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1136 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 1137 "MB01VD.f"
				lc = ic;

#line 1139 "MB01VD.f"
				i__4 = *mb;
#line 1139 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1140 "MB01VD.f"
				    c__[lc + jc * c_dim1] += aij * b[k + l * 
					    b_dim1];
#line 1141 "MB01VD.f"
				    ++lc;
#line 1142 "MB01VD.f"
/* L1250: */
#line 1142 "MB01VD.f"
				}

#line 1144 "MB01VD.f"
				ic += *mb;
#line 1145 "MB01VD.f"
/* L1260: */
#line 1145 "MB01VD.f"
			    }

#line 1147 "MB01VD.f"
			    ++jc;
#line 1148 "MB01VD.f"
/* L1270: */
#line 1148 "MB01VD.f"
			}

#line 1150 "MB01VD.f"
/* L1280: */
#line 1150 "MB01VD.f"
		    }

#line 1152 "MB01VD.f"
		}
#line 1153 "MB01VD.f"
	    } else {
#line 1154 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 1, alpha <> 1, A sparse. */

#line 1158 "MB01VD.f"
		    i__1 = *na;
#line 1158 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1160 "MB01VD.f"
			i__2 = *nb;
#line 1160 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1161 "MB01VD.f"
			    ic = 1;

#line 1163 "MB01VD.f"
			    i__3 = *ma;
#line 1163 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1164 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 1165 "MB01VD.f"
				if (aij != 0.) {
#line 1166 "MB01VD.f"
				    lc = ic;

#line 1168 "MB01VD.f"
				    i__4 = *mb;
#line 1168 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1169 "MB01VD.f"
					c__[lc + jc * c_dim1] += aij * b[k + 
						l * b_dim1];
#line 1170 "MB01VD.f"
					++lc;
#line 1171 "MB01VD.f"
/* L1290: */
#line 1171 "MB01VD.f"
				    }

#line 1173 "MB01VD.f"
				}
#line 1174 "MB01VD.f"
				ic += *mb;
#line 1175 "MB01VD.f"
/* L1300: */
#line 1175 "MB01VD.f"
			    }

#line 1177 "MB01VD.f"
			    ++jc;
#line 1178 "MB01VD.f"
/* L1310: */
#line 1178 "MB01VD.f"
			}

#line 1180 "MB01VD.f"
/* L1320: */
#line 1180 "MB01VD.f"
		    }

#line 1182 "MB01VD.f"
		} else {

/*                 Case beta = 1, alpha <> 1, A not sparse. */

#line 1186 "MB01VD.f"
		    i__1 = *na;
#line 1186 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1188 "MB01VD.f"
			i__2 = *nb;
#line 1188 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1189 "MB01VD.f"
			    ic = 1;

#line 1191 "MB01VD.f"
			    i__3 = *ma;
#line 1191 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1192 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 1193 "MB01VD.f"
				lc = ic;

#line 1195 "MB01VD.f"
				i__4 = *mb;
#line 1195 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1196 "MB01VD.f"
				    c__[lc + jc * c_dim1] += aij * b[k + l * 
					    b_dim1];
#line 1197 "MB01VD.f"
				    ++lc;
#line 1198 "MB01VD.f"
/* L1330: */
#line 1198 "MB01VD.f"
				}

#line 1200 "MB01VD.f"
				ic += *mb;
#line 1201 "MB01VD.f"
/* L1340: */
#line 1201 "MB01VD.f"
			    }

#line 1203 "MB01VD.f"
			    ++jc;
#line 1204 "MB01VD.f"
/* L1350: */
#line 1204 "MB01VD.f"
			}

#line 1206 "MB01VD.f"
/* L1360: */
#line 1206 "MB01VD.f"
		    }

#line 1208 "MB01VD.f"
		}
#line 1209 "MB01VD.f"
	    }
#line 1210 "MB01VD.f"
	} else {
#line 1211 "MB01VD.f"
	    if (*alpha == 1.) {
#line 1212 "MB01VD.f"
		if (sparse) {

/*                 Case beta <> 0 or 1, alpha = 1, A sparse. */

#line 1216 "MB01VD.f"
		    i__1 = *na;
#line 1216 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1218 "MB01VD.f"
			i__2 = *nb;
#line 1218 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1219 "MB01VD.f"
			    ic = 1;

#line 1221 "MB01VD.f"
			    i__3 = *ma;
#line 1221 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1222 "MB01VD.f"
				aij = a[i__ + j * a_dim1];

#line 1224 "MB01VD.f"
				if (aij == 0.) {
#line 1225 "MB01VD.f"
				    dscal_(mb, beta, &c__[ic + jc * c_dim1], &
					    c__1);
#line 1226 "MB01VD.f"
				} else {
#line 1227 "MB01VD.f"
				    lc = ic;

#line 1229 "MB01VD.f"
				    i__4 = *mb;
#line 1229 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1230 "MB01VD.f"
					c__[lc + jc * c_dim1] = *beta * c__[
						lc + jc * c_dim1] + aij * b[k 
						+ l * b_dim1];
#line 1231 "MB01VD.f"
					++lc;
#line 1232 "MB01VD.f"
/* L1370: */
#line 1232 "MB01VD.f"
				    }

#line 1234 "MB01VD.f"
				}
#line 1235 "MB01VD.f"
				ic += *mb;
#line 1236 "MB01VD.f"
/* L1380: */
#line 1236 "MB01VD.f"
			    }

#line 1238 "MB01VD.f"
			    ++jc;
#line 1239 "MB01VD.f"
/* L1390: */
#line 1239 "MB01VD.f"
			}

#line 1241 "MB01VD.f"
/* L1400: */
#line 1241 "MB01VD.f"
		    }

#line 1243 "MB01VD.f"
		} else {

/*                 Case beta <> 0 or 1, alpha = 1, A not sparse. */

#line 1247 "MB01VD.f"
		    i__1 = *na;
#line 1247 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1249 "MB01VD.f"
			i__2 = *nb;
#line 1249 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1250 "MB01VD.f"
			    ic = 1;

#line 1252 "MB01VD.f"
			    i__3 = *ma;
#line 1252 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1253 "MB01VD.f"
				aij = a[i__ + j * a_dim1];
#line 1254 "MB01VD.f"
				lc = ic;

#line 1256 "MB01VD.f"
				i__4 = *mb;
#line 1256 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1257 "MB01VD.f"
				    c__[lc + jc * c_dim1] = *beta * c__[lc + 
					    jc * c_dim1] + aij * b[k + l * 
					    b_dim1];
#line 1258 "MB01VD.f"
				    ++lc;
#line 1259 "MB01VD.f"
/* L1410: */
#line 1259 "MB01VD.f"
				}

#line 1261 "MB01VD.f"
				ic += *mb;
#line 1262 "MB01VD.f"
/* L1420: */
#line 1262 "MB01VD.f"
			    }

#line 1264 "MB01VD.f"
			    ++jc;
#line 1265 "MB01VD.f"
/* L1430: */
#line 1265 "MB01VD.f"
			}

#line 1267 "MB01VD.f"
/* L1440: */
#line 1267 "MB01VD.f"
		    }

#line 1269 "MB01VD.f"
		}
#line 1270 "MB01VD.f"
	    } else {
#line 1271 "MB01VD.f"
		if (sparse) {

/*                 Case beta <> 0 or 1, alpha <> 1, A sparse. */

#line 1275 "MB01VD.f"
		    i__1 = *na;
#line 1275 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1277 "MB01VD.f"
			i__2 = *nb;
#line 1277 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1278 "MB01VD.f"
			    ic = 1;

#line 1280 "MB01VD.f"
			    i__3 = *ma;
#line 1280 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1281 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];

#line 1283 "MB01VD.f"
				if (aij == 0.) {
#line 1284 "MB01VD.f"
				    dscal_(mb, beta, &c__[ic + jc * c_dim1], &
					    c__1);
#line 1285 "MB01VD.f"
				} else {
#line 1286 "MB01VD.f"
				    lc = ic;

#line 1288 "MB01VD.f"
				    i__4 = *mb;
#line 1288 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1289 "MB01VD.f"
					c__[lc + jc * c_dim1] = *beta * c__[
						lc + jc * c_dim1] + aij * b[k 
						+ l * b_dim1];
#line 1290 "MB01VD.f"
					++lc;
#line 1291 "MB01VD.f"
/* L1450: */
#line 1291 "MB01VD.f"
				    }

#line 1293 "MB01VD.f"
				}
#line 1294 "MB01VD.f"
				ic += *mb;
#line 1295 "MB01VD.f"
/* L1460: */
#line 1295 "MB01VD.f"
			    }

#line 1297 "MB01VD.f"
			    ++jc;
#line 1298 "MB01VD.f"
/* L1470: */
#line 1298 "MB01VD.f"
			}

#line 1300 "MB01VD.f"
/* L1480: */
#line 1300 "MB01VD.f"
		    }

#line 1302 "MB01VD.f"
		} else {

/*                 Case beta <> 0 or 1, alpha <> 1, A not sparse. */

#line 1306 "MB01VD.f"
		    i__1 = *na;
#line 1306 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1308 "MB01VD.f"
			i__2 = *nb;
#line 1308 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1309 "MB01VD.f"
			    ic = 1;

#line 1311 "MB01VD.f"
			    i__3 = *ma;
#line 1311 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1312 "MB01VD.f"
				aij = *alpha * a[i__ + j * a_dim1];
#line 1313 "MB01VD.f"
				lc = ic;

#line 1315 "MB01VD.f"
				i__4 = *mb;
#line 1315 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1316 "MB01VD.f"
				    c__[lc + jc * c_dim1] = *beta * c__[lc + 
					    jc * c_dim1] + aij * b[k + l * 
					    b_dim1];
#line 1317 "MB01VD.f"
				    ++lc;
#line 1318 "MB01VD.f"
/* L1490: */
#line 1318 "MB01VD.f"
				}

#line 1320 "MB01VD.f"
				ic += *mb;
#line 1321 "MB01VD.f"
/* L1500: */
#line 1321 "MB01VD.f"
			    }

#line 1323 "MB01VD.f"
			    ++jc;
#line 1324 "MB01VD.f"
/* L1510: */
#line 1324 "MB01VD.f"
			}

#line 1326 "MB01VD.f"
/* L1520: */
#line 1326 "MB01VD.f"
		    }

#line 1328 "MB01VD.f"
		}
#line 1329 "MB01VD.f"
	    }
#line 1330 "MB01VD.f"
	}
#line 1331 "MB01VD.f"
    } else {

/*        Case op(A) = A' and op(B) = B'. */

#line 1335 "MB01VD.f"
	if (*beta == 0.) {
#line 1336 "MB01VD.f"
	    if (*alpha == 1.) {
#line 1337 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 0, alpha = 1, A sparse. */

#line 1341 "MB01VD.f"
		    i__1 = *na;
#line 1341 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1343 "MB01VD.f"
			i__2 = *nb;
#line 1343 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1344 "MB01VD.f"
			    ic = 1;

#line 1346 "MB01VD.f"
			    i__3 = *ma;
#line 1346 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1347 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 1348 "MB01VD.f"
				if (aij == 0.) {
#line 1349 "MB01VD.f"
				    dcopy_(mb, dum, &c__0, &c__[ic + jc * 
					    c_dim1], &c__1);
#line 1350 "MB01VD.f"
				} else if (aij == 1.) {
#line 1351 "MB01VD.f"
				    dcopy_(mb, &b[k + b_dim1], ldb, &c__[ic + 
					    jc * c_dim1], &c__1);
#line 1352 "MB01VD.f"
				} else {
#line 1353 "MB01VD.f"
				    lc = ic;

#line 1355 "MB01VD.f"
				    i__4 = *mb;
#line 1355 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1356 "MB01VD.f"
					c__[lc + jc * c_dim1] = aij * b[k + l 
						* b_dim1];
#line 1357 "MB01VD.f"
					++lc;
#line 1358 "MB01VD.f"
/* L1550: */
#line 1358 "MB01VD.f"
				    }

#line 1360 "MB01VD.f"
				}
#line 1361 "MB01VD.f"
				ic += *mb;
#line 1362 "MB01VD.f"
/* L1560: */
#line 1362 "MB01VD.f"
			    }

#line 1364 "MB01VD.f"
			    ++jc;
#line 1365 "MB01VD.f"
/* L1570: */
#line 1365 "MB01VD.f"
			}

#line 1367 "MB01VD.f"
/* L1580: */
#line 1367 "MB01VD.f"
		    }

#line 1369 "MB01VD.f"
		} else {

/*                 Case beta = 0, alpha = 1, A not sparse. */

#line 1373 "MB01VD.f"
		    i__1 = *na;
#line 1373 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1375 "MB01VD.f"
			i__2 = *nb;
#line 1375 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1376 "MB01VD.f"
			    ic = 1;

#line 1378 "MB01VD.f"
			    i__3 = *ma;
#line 1378 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1379 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 1380 "MB01VD.f"
				lc = ic;

#line 1382 "MB01VD.f"
				i__4 = *mb;
#line 1382 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1383 "MB01VD.f"
				    c__[lc + jc * c_dim1] = aij * b[k + l * 
					    b_dim1];
#line 1384 "MB01VD.f"
				    ++lc;
#line 1385 "MB01VD.f"
/* L1590: */
#line 1385 "MB01VD.f"
				}

#line 1387 "MB01VD.f"
				ic += *mb;
#line 1388 "MB01VD.f"
/* L1600: */
#line 1388 "MB01VD.f"
			    }

#line 1390 "MB01VD.f"
			    ++jc;
#line 1391 "MB01VD.f"
/* L1610: */
#line 1391 "MB01VD.f"
			}

#line 1393 "MB01VD.f"
/* L1620: */
#line 1393 "MB01VD.f"
		    }

#line 1395 "MB01VD.f"
		}
#line 1396 "MB01VD.f"
	    } else {
#line 1397 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 0, alpha <> 1, A sparse. */

#line 1401 "MB01VD.f"
		    i__1 = *na;
#line 1401 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1403 "MB01VD.f"
			i__2 = *nb;
#line 1403 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1404 "MB01VD.f"
			    ic = 1;

#line 1406 "MB01VD.f"
			    i__3 = *ma;
#line 1406 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1407 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 1408 "MB01VD.f"
				if (aij == 0.) {
#line 1409 "MB01VD.f"
				    dcopy_(mb, dum, &c__0, &c__[ic + jc * 
					    c_dim1], &c__1);
#line 1410 "MB01VD.f"
				} else {
#line 1411 "MB01VD.f"
				    lc = ic;

#line 1413 "MB01VD.f"
				    i__4 = *mb;
#line 1413 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1414 "MB01VD.f"
					c__[lc + jc * c_dim1] = aij * b[k + l 
						* b_dim1];
#line 1415 "MB01VD.f"
					++lc;
#line 1416 "MB01VD.f"
/* L1630: */
#line 1416 "MB01VD.f"
				    }

#line 1418 "MB01VD.f"
				}
#line 1419 "MB01VD.f"
				ic += *mb;
#line 1420 "MB01VD.f"
/* L1640: */
#line 1420 "MB01VD.f"
			    }

#line 1422 "MB01VD.f"
			    ++jc;
#line 1423 "MB01VD.f"
/* L1650: */
#line 1423 "MB01VD.f"
			}

#line 1425 "MB01VD.f"
/* L1660: */
#line 1425 "MB01VD.f"
		    }

#line 1427 "MB01VD.f"
		} else {

/*                 Case beta = 0, alpha <> 1, A not sparse. */

#line 1431 "MB01VD.f"
		    i__1 = *na;
#line 1431 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1433 "MB01VD.f"
			i__2 = *nb;
#line 1433 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1434 "MB01VD.f"
			    ic = 1;

#line 1436 "MB01VD.f"
			    i__3 = *ma;
#line 1436 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1437 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 1438 "MB01VD.f"
				lc = ic;

#line 1440 "MB01VD.f"
				i__4 = *mb;
#line 1440 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1441 "MB01VD.f"
				    c__[lc + jc * c_dim1] = aij * b[k + l * 
					    b_dim1];
#line 1442 "MB01VD.f"
				    ++lc;
#line 1443 "MB01VD.f"
/* L1670: */
#line 1443 "MB01VD.f"
				}

#line 1445 "MB01VD.f"
				ic += *mb;
#line 1446 "MB01VD.f"
/* L1680: */
#line 1446 "MB01VD.f"
			    }

#line 1448 "MB01VD.f"
			    ++jc;
#line 1449 "MB01VD.f"
/* L1690: */
#line 1449 "MB01VD.f"
			}

#line 1451 "MB01VD.f"
/* L1700: */
#line 1451 "MB01VD.f"
		    }

#line 1453 "MB01VD.f"
		}
#line 1454 "MB01VD.f"
	    }
#line 1455 "MB01VD.f"
	} else if (*beta == 1.) {
#line 1456 "MB01VD.f"
	    if (*alpha == 1.) {
#line 1457 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 1, alpha = 1, A sparse. */

#line 1461 "MB01VD.f"
		    i__1 = *na;
#line 1461 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1463 "MB01VD.f"
			i__2 = *nb;
#line 1463 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1464 "MB01VD.f"
			    ic = 1;

#line 1466 "MB01VD.f"
			    i__3 = *ma;
#line 1466 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1467 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 1468 "MB01VD.f"
				if (aij != 0.) {
#line 1469 "MB01VD.f"
				    lc = ic;

#line 1471 "MB01VD.f"
				    i__4 = *mb;
#line 1471 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1472 "MB01VD.f"
					c__[lc + jc * c_dim1] += aij * b[k + 
						l * b_dim1];
#line 1473 "MB01VD.f"
					++lc;
#line 1474 "MB01VD.f"
/* L1710: */
#line 1474 "MB01VD.f"
				    }

#line 1476 "MB01VD.f"
				}
#line 1477 "MB01VD.f"
				ic += *mb;
#line 1478 "MB01VD.f"
/* L1720: */
#line 1478 "MB01VD.f"
			    }

#line 1480 "MB01VD.f"
			    ++jc;
#line 1481 "MB01VD.f"
/* L1730: */
#line 1481 "MB01VD.f"
			}

#line 1483 "MB01VD.f"
/* L1740: */
#line 1483 "MB01VD.f"
		    }

#line 1485 "MB01VD.f"
		} else {

/*                 Case beta = 1, alpha = 1, A not sparse. */

#line 1489 "MB01VD.f"
		    i__1 = *na;
#line 1489 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1491 "MB01VD.f"
			i__2 = *nb;
#line 1491 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1492 "MB01VD.f"
			    ic = 1;

#line 1494 "MB01VD.f"
			    i__3 = *ma;
#line 1494 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1495 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 1496 "MB01VD.f"
				lc = ic;

#line 1498 "MB01VD.f"
				i__4 = *mb;
#line 1498 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1499 "MB01VD.f"
				    c__[lc + jc * c_dim1] += aij * b[k + l * 
					    b_dim1];
#line 1500 "MB01VD.f"
				    ++lc;
#line 1501 "MB01VD.f"
/* L1750: */
#line 1501 "MB01VD.f"
				}

#line 1503 "MB01VD.f"
				ic += *mb;
#line 1504 "MB01VD.f"
/* L1760: */
#line 1504 "MB01VD.f"
			    }

#line 1506 "MB01VD.f"
			    ++jc;
#line 1507 "MB01VD.f"
/* L1770: */
#line 1507 "MB01VD.f"
			}

#line 1509 "MB01VD.f"
/* L1780: */
#line 1509 "MB01VD.f"
		    }

#line 1511 "MB01VD.f"
		}
#line 1512 "MB01VD.f"
	    } else {
#line 1513 "MB01VD.f"
		if (sparse) {

/*                 Case beta = 1, alpha <> 1, A sparse. */

#line 1517 "MB01VD.f"
		    i__1 = *na;
#line 1517 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1519 "MB01VD.f"
			i__2 = *nb;
#line 1519 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1520 "MB01VD.f"
			    ic = 1;

#line 1522 "MB01VD.f"
			    i__3 = *ma;
#line 1522 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1523 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 1524 "MB01VD.f"
				if (aij != 0.) {
#line 1525 "MB01VD.f"
				    lc = ic;

#line 1527 "MB01VD.f"
				    i__4 = *mb;
#line 1527 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1528 "MB01VD.f"
					c__[lc + jc * c_dim1] += aij * b[k + 
						l * b_dim1];
#line 1529 "MB01VD.f"
					++lc;
#line 1530 "MB01VD.f"
/* L1790: */
#line 1530 "MB01VD.f"
				    }

#line 1532 "MB01VD.f"
				}
#line 1533 "MB01VD.f"
				ic += *mb;
#line 1534 "MB01VD.f"
/* L1800: */
#line 1534 "MB01VD.f"
			    }

#line 1536 "MB01VD.f"
			    ++jc;
#line 1537 "MB01VD.f"
/* L1810: */
#line 1537 "MB01VD.f"
			}

#line 1539 "MB01VD.f"
/* L1820: */
#line 1539 "MB01VD.f"
		    }

#line 1541 "MB01VD.f"
		} else {

/*                 Case beta = 1, alpha <> 1, A not sparse. */

#line 1545 "MB01VD.f"
		    i__1 = *na;
#line 1545 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1547 "MB01VD.f"
			i__2 = *nb;
#line 1547 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1548 "MB01VD.f"
			    ic = 1;

#line 1550 "MB01VD.f"
			    i__3 = *ma;
#line 1550 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1551 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 1552 "MB01VD.f"
				lc = ic;

#line 1554 "MB01VD.f"
				i__4 = *mb;
#line 1554 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1555 "MB01VD.f"
				    c__[lc + jc * c_dim1] += aij * b[k + l * 
					    b_dim1];
#line 1556 "MB01VD.f"
				    ++lc;
#line 1557 "MB01VD.f"
/* L1830: */
#line 1557 "MB01VD.f"
				}

#line 1559 "MB01VD.f"
				ic += *mb;
#line 1560 "MB01VD.f"
/* L1840: */
#line 1560 "MB01VD.f"
			    }

#line 1562 "MB01VD.f"
			    ++jc;
#line 1563 "MB01VD.f"
/* L1850: */
#line 1563 "MB01VD.f"
			}

#line 1565 "MB01VD.f"
/* L1860: */
#line 1565 "MB01VD.f"
		    }

#line 1567 "MB01VD.f"
		}
#line 1568 "MB01VD.f"
	    }
#line 1569 "MB01VD.f"
	} else {
#line 1570 "MB01VD.f"
	    if (*alpha == 1.) {
#line 1571 "MB01VD.f"
		if (sparse) {

/*                 Case beta <> 0 or 1, alpha = 1, A sparse. */

#line 1575 "MB01VD.f"
		    i__1 = *na;
#line 1575 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1577 "MB01VD.f"
			i__2 = *nb;
#line 1577 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1578 "MB01VD.f"
			    ic = 1;

#line 1580 "MB01VD.f"
			    i__3 = *ma;
#line 1580 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1581 "MB01VD.f"
				aij = a[j + i__ * a_dim1];

#line 1583 "MB01VD.f"
				if (aij == 0.) {
#line 1584 "MB01VD.f"
				    dscal_(mb, beta, &c__[ic + jc * c_dim1], &
					    c__1);
#line 1585 "MB01VD.f"
				} else {
#line 1586 "MB01VD.f"
				    lc = ic;

#line 1588 "MB01VD.f"
				    i__4 = *mb;
#line 1588 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1589 "MB01VD.f"
					c__[lc + jc * c_dim1] = *beta * c__[
						lc + jc * c_dim1] + aij * b[k 
						+ l * b_dim1];
#line 1590 "MB01VD.f"
					++lc;
#line 1591 "MB01VD.f"
/* L1870: */
#line 1591 "MB01VD.f"
				    }

#line 1593 "MB01VD.f"
				}
#line 1594 "MB01VD.f"
				ic += *mb;
#line 1595 "MB01VD.f"
/* L1880: */
#line 1595 "MB01VD.f"
			    }

#line 1597 "MB01VD.f"
			    ++jc;
#line 1598 "MB01VD.f"
/* L1890: */
#line 1598 "MB01VD.f"
			}

#line 1600 "MB01VD.f"
/* L1900: */
#line 1600 "MB01VD.f"
		    }

#line 1602 "MB01VD.f"
		} else {

/*                 Case beta <> 0 or 1, alpha = 1, A not sparse. */

#line 1606 "MB01VD.f"
		    i__1 = *na;
#line 1606 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1608 "MB01VD.f"
			i__2 = *nb;
#line 1608 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1609 "MB01VD.f"
			    ic = 1;

#line 1611 "MB01VD.f"
			    i__3 = *ma;
#line 1611 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1612 "MB01VD.f"
				aij = a[j + i__ * a_dim1];
#line 1613 "MB01VD.f"
				lc = ic;

#line 1615 "MB01VD.f"
				i__4 = *mb;
#line 1615 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1616 "MB01VD.f"
				    c__[lc + jc * c_dim1] = *beta * c__[lc + 
					    jc * c_dim1] + aij * b[k + l * 
					    b_dim1];
#line 1617 "MB01VD.f"
				    ++lc;
#line 1618 "MB01VD.f"
/* L1910: */
#line 1618 "MB01VD.f"
				}

#line 1620 "MB01VD.f"
				ic += *mb;
#line 1621 "MB01VD.f"
/* L1920: */
#line 1621 "MB01VD.f"
			    }

#line 1623 "MB01VD.f"
			    ++jc;
#line 1624 "MB01VD.f"
/* L1930: */
#line 1624 "MB01VD.f"
			}

#line 1626 "MB01VD.f"
/* L1940: */
#line 1626 "MB01VD.f"
		    }

#line 1628 "MB01VD.f"
		}
#line 1629 "MB01VD.f"
	    } else {
#line 1630 "MB01VD.f"
		if (sparse) {

/*                 Case beta <> 0 or 1, alpha <> 1, A sparse. */

#line 1634 "MB01VD.f"
		    i__1 = *na;
#line 1634 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1636 "MB01VD.f"
			i__2 = *nb;
#line 1636 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1637 "MB01VD.f"
			    ic = 1;

#line 1639 "MB01VD.f"
			    i__3 = *ma;
#line 1639 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1640 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];

#line 1642 "MB01VD.f"
				if (aij == 0.) {
#line 1643 "MB01VD.f"
				    dscal_(mb, beta, &c__[ic + jc * c_dim1], &
					    c__1);
#line 1644 "MB01VD.f"
				} else {
#line 1645 "MB01VD.f"
				    lc = ic;

#line 1647 "MB01VD.f"
				    i__4 = *mb;
#line 1647 "MB01VD.f"
				    for (l = 1; l <= i__4; ++l) {
#line 1648 "MB01VD.f"
					c__[lc + jc * c_dim1] = *beta * c__[
						lc + jc * c_dim1] + aij * b[k 
						+ l * b_dim1];
#line 1649 "MB01VD.f"
					++lc;
#line 1650 "MB01VD.f"
/* L1950: */
#line 1650 "MB01VD.f"
				    }

#line 1652 "MB01VD.f"
				}
#line 1653 "MB01VD.f"
				ic += *mb;
#line 1654 "MB01VD.f"
/* L1960: */
#line 1654 "MB01VD.f"
			    }

#line 1656 "MB01VD.f"
			    ++jc;
#line 1657 "MB01VD.f"
/* L1970: */
#line 1657 "MB01VD.f"
			}

#line 1659 "MB01VD.f"
/* L1980: */
#line 1659 "MB01VD.f"
		    }

#line 1661 "MB01VD.f"
		} else {

/*                 Case beta <> 0 or 1, alpha <> 1, A not sparse. */

#line 1665 "MB01VD.f"
		    i__1 = *na;
#line 1665 "MB01VD.f"
		    for (j = 1; j <= i__1; ++j) {

#line 1667 "MB01VD.f"
			i__2 = *nb;
#line 1667 "MB01VD.f"
			for (k = 1; k <= i__2; ++k) {
#line 1668 "MB01VD.f"
			    ic = 1;

#line 1670 "MB01VD.f"
			    i__3 = *ma;
#line 1670 "MB01VD.f"
			    for (i__ = 1; i__ <= i__3; ++i__) {
#line 1671 "MB01VD.f"
				aij = *alpha * a[j + i__ * a_dim1];
#line 1672 "MB01VD.f"
				lc = ic;

#line 1674 "MB01VD.f"
				i__4 = *mb;
#line 1674 "MB01VD.f"
				for (l = 1; l <= i__4; ++l) {
#line 1675 "MB01VD.f"
				    c__[lc + jc * c_dim1] = *beta * c__[lc + 
					    jc * c_dim1] + aij * b[k + l * 
					    b_dim1];
#line 1676 "MB01VD.f"
				    ++lc;
#line 1677 "MB01VD.f"
/* L1990: */
#line 1677 "MB01VD.f"
				}

#line 1679 "MB01VD.f"
				ic += *mb;
#line 1680 "MB01VD.f"
/* L2000: */
#line 1680 "MB01VD.f"
			    }

#line 1682 "MB01VD.f"
			    ++jc;
#line 1683 "MB01VD.f"
/* L2010: */
#line 1683 "MB01VD.f"
			}

#line 1685 "MB01VD.f"
/* L2020: */
#line 1685 "MB01VD.f"
		    }

#line 1687 "MB01VD.f"
		}
#line 1688 "MB01VD.f"
	    }
#line 1689 "MB01VD.f"
	}
#line 1690 "MB01VD.f"
    }
#line 1691 "MB01VD.f"
    return 0;
/* *** Last line of MB01VD *** */
} /* mb01vd_ */

