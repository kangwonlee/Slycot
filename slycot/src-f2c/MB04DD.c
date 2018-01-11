#line 1 "MB04DD.f"
/* MB04DD.f -- translated by f2c (version 20100827).
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

#line 1 "MB04DD.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b29 = -1.;

/* Subroutine */ int mb04dd_(char *job, integer *n, doublereal *a, integer *
	lda, doublereal *qg, integer *ldqg, integer *ilo, doublereal *scale, 
	integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2, i__3, i__4, 
	    i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6;

    /* Local variables */
    static doublereal c__, f;
    static integer i__, j;
    static doublereal r__;
    static integer ic;
    static doublereal gii, qii, maxc;
    static logical conv;
    static doublereal temp, maxr;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static logical lscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical lperm;
    static doublereal sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal sclfac;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer iloold;


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

/*     To balance a real Hamiltonian matrix, */

/*                   [  A   G  ] */
/*              H =  [       T ] , */
/*                   [  Q  -A  ] */

/*     where A is an N-by-N matrix and G, Q are N-by-N symmetric */
/*     matrices. This involves, first, permuting H by a symplectic */
/*     similarity transformation to isolate eigenvalues in the first */
/*     1:ILO-1 elements on the diagonal of A; and second, applying a */
/*     diagonal similarity transformation to rows and columns */
/*     ILO:2*N-ILO+1 to make the rows and columns as close in 1-norm */
/*     as possible. Both steps are optional. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the operations to be performed on H: */
/*             = 'N':  none, set ILO = 1, SCALE(I) = 1.0, I = 1 .. N; */
/*             = 'P':  permute only; */
/*             = 'S':  scale only; */
/*             = 'B':  both permute and scale. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix A of the balanced Hamiltonian. In particular, */
/*             the lower triangular part of the first ILO-1 columns of A */
/*             is zero. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*                            (LDQG,N+1) */
/*             On entry, the leading N-by-N+1 part of this array must */
/*             contain the lower triangular part of the matrix Q and */
/*             the upper triangular part of the matrix G. */
/*             On exit, the leading N-by-N+1 part of this array contains */
/*             the lower and upper triangular parts of the matrices Q and */
/*             G, respectively, of the balanced Hamiltonian. In */
/*             particular, the lower triangular and diagonal part of the */
/*             first ILO-1 columns of QG is zero. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= MAX(1,N). */

/*     ILO     (output) INTEGER */
/*             ILO-1 is the number of deflated eigenvalues in the */
/*             balanced Hamiltonian matrix. */

/*     SCALE   (output) DOUBLE PRECISION array of dimension (N) */
/*             Details of the permutations and scaling factors applied to */
/*             H.  For j = 1,...,ILO-1 let P(j) = SCALE(j). If P(j) <= N, */
/*             then rows and columns P(j) and P(j)+N are interchanged */
/*             with rows and columns j and j+N, respectively. If */
/*             P(j) > N, then row and column P(j)-N are interchanged with */
/*             row and column j+N by a generalized symplectic */
/*             permutation. For j = ILO,...,N the j-th element of SCALE */
/*             contains the factor of the scaling applied to row and */
/*             column j. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     REFERENCES */

/*     [1] Benner, P. */
/*         Symplectic balancing of Hamiltonian matrices. */
/*         SIAM J. Sci. Comput., 22 (5), pp. 1885-1904, 2000. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DHABAL). */

/*     KEYWORDS */

/*     Balancing, Hamiltonian matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 148 "MB04DD.f"
    /* Parameter adjustments */
#line 148 "MB04DD.f"
    a_dim1 = *lda;
#line 148 "MB04DD.f"
    a_offset = 1 + a_dim1;
#line 148 "MB04DD.f"
    a -= a_offset;
#line 148 "MB04DD.f"
    qg_dim1 = *ldqg;
#line 148 "MB04DD.f"
    qg_offset = 1 + qg_dim1;
#line 148 "MB04DD.f"
    qg -= qg_offset;
#line 148 "MB04DD.f"
    --scale;
#line 148 "MB04DD.f"

#line 148 "MB04DD.f"
    /* Function Body */
#line 148 "MB04DD.f"
    *info = 0;
#line 149 "MB04DD.f"
    lperm = lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
#line 150 "MB04DD.f"
    lscal = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);

#line 152 "MB04DD.f"
    if (! lperm && ! lscal && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 154 "MB04DD.f"
	*info = -1;
#line 155 "MB04DD.f"
    } else if (*n < 0) {
#line 156 "MB04DD.f"
	*info = -2;
#line 157 "MB04DD.f"
    } else if (*lda < max(1,*n)) {
#line 158 "MB04DD.f"
	*info = -4;
#line 159 "MB04DD.f"
    } else if (*ldqg < max(1,*n)) {
#line 160 "MB04DD.f"
	*info = -6;
#line 161 "MB04DD.f"
    }

/*     Return if there were illegal values. */

#line 165 "MB04DD.f"
    if (*info != 0) {
#line 166 "MB04DD.f"
	i__1 = -(*info);
#line 166 "MB04DD.f"
	xerbla_("MB04DD", &i__1, (ftnlen)6);
#line 167 "MB04DD.f"
	return 0;
#line 168 "MB04DD.f"
    }

#line 170 "MB04DD.f"
    *ilo = 1;

/*     Quick return if possible. */

#line 174 "MB04DD.f"
    if (*n == 0) {
#line 174 "MB04DD.f"
	return 0;
#line 174 "MB04DD.f"
    }
#line 176 "MB04DD.f"
    if (! lperm && ! lscal) {
#line 177 "MB04DD.f"
	i__1 = *n;
#line 177 "MB04DD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 178 "MB04DD.f"
	    scale[i__] = 1.;
#line 179 "MB04DD.f"
/* L10: */
#line 179 "MB04DD.f"
	}
#line 180 "MB04DD.f"
	return 0;
#line 181 "MB04DD.f"
    }

/*     Permutations to isolate eigenvalues if possible. */

#line 185 "MB04DD.f"
    if (lperm) {
#line 186 "MB04DD.f"
	iloold = 0;
/*        WHILE ( ILO.NE.ILOOLD ) */
#line 188 "MB04DD.f"
L20:
#line 188 "MB04DD.f"
	if (*ilo != iloold) {
#line 189 "MB04DD.f"
	    iloold = *ilo;

/*           Scan columns ILO .. N. */

#line 193 "MB04DD.f"
	    i__ = *ilo;
/*           WHILE ( I.LE.N .AND. ILO.EQ.ILOOLD ) */
#line 195 "MB04DD.f"
L30:
#line 195 "MB04DD.f"
	    if (i__ <= *n && *ilo == iloold) {
#line 196 "MB04DD.f"
		i__1 = i__ - 1;
#line 196 "MB04DD.f"
		for (j = *ilo; j <= i__1; ++j) {
#line 197 "MB04DD.f"
		    if (a[j + i__ * a_dim1] != 0.) {
#line 198 "MB04DD.f"
			++i__;
#line 199 "MB04DD.f"
			goto L30;
#line 200 "MB04DD.f"
		    }
#line 201 "MB04DD.f"
/* L40: */
#line 201 "MB04DD.f"
		}
#line 202 "MB04DD.f"
		i__1 = *n;
#line 202 "MB04DD.f"
		for (j = i__ + 1; j <= i__1; ++j) {
#line 203 "MB04DD.f"
		    if (a[j + i__ * a_dim1] != 0.) {
#line 204 "MB04DD.f"
			++i__;
#line 205 "MB04DD.f"
			goto L30;
#line 206 "MB04DD.f"
		    }
#line 207 "MB04DD.f"
/* L50: */
#line 207 "MB04DD.f"
		}
#line 208 "MB04DD.f"
		i__1 = i__;
#line 208 "MB04DD.f"
		for (j = *ilo; j <= i__1; ++j) {
#line 209 "MB04DD.f"
		    if (qg[i__ + j * qg_dim1] != 0.) {
#line 210 "MB04DD.f"
			++i__;
#line 211 "MB04DD.f"
			goto L30;
#line 212 "MB04DD.f"
		    }
#line 213 "MB04DD.f"
/* L60: */
#line 213 "MB04DD.f"
		}
#line 214 "MB04DD.f"
		i__1 = *n;
#line 214 "MB04DD.f"
		for (j = i__ + 1; j <= i__1; ++j) {
#line 215 "MB04DD.f"
		    if (qg[j + i__ * qg_dim1] != 0.) {
#line 216 "MB04DD.f"
			++i__;
#line 217 "MB04DD.f"
			goto L30;
#line 218 "MB04DD.f"
		    }
#line 219 "MB04DD.f"
/* L70: */
#line 219 "MB04DD.f"
		}

/*              Exchange columns/rows ILO <-> I. */

#line 223 "MB04DD.f"
		scale[*ilo] = (doublereal) i__;
#line 224 "MB04DD.f"
		if (*ilo != i__) {

#line 226 "MB04DD.f"
		    dswap_(n, &a[*ilo * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 
			    1], &c__1);
#line 227 "MB04DD.f"
		    i__1 = *n - *ilo + 1;
#line 227 "MB04DD.f"
		    dswap_(&i__1, &a[*ilo + *ilo * a_dim1], lda, &a[i__ + *
			    ilo * a_dim1], lda);

#line 229 "MB04DD.f"
		    dswap_(&c__1, &qg[i__ + *ilo * qg_dim1], ldqg, &qg[*ilo + 
			    *ilo * qg_dim1], ldqg);
#line 230 "MB04DD.f"
		    i__1 = *n - i__ + 1;
#line 230 "MB04DD.f"
		    dswap_(&i__1, &qg[i__ + i__ * qg_dim1], &c__1, &qg[i__ + *
			    ilo * qg_dim1], &c__1);
#line 231 "MB04DD.f"
		    i__1 = i__ - *ilo;
#line 231 "MB04DD.f"
		    dswap_(&i__1, &qg[*ilo + *ilo * qg_dim1], &c__1, &qg[i__ 
			    + *ilo * qg_dim1], ldqg);

#line 233 "MB04DD.f"
		    dswap_(ilo, &qg[(i__ + 1) * qg_dim1 + 1], &c__1, &qg[(*
			    ilo + 1) * qg_dim1 + 1], &c__1);
#line 234 "MB04DD.f"
		    i__1 = *n - i__ + 1;
#line 234 "MB04DD.f"
		    dswap_(&i__1, &qg[i__ + (i__ + 1) * qg_dim1], ldqg, &qg[*
			    ilo + (i__ + 1) * qg_dim1], ldqg);
#line 236 "MB04DD.f"
		    i__1 = i__ - *ilo;
#line 236 "MB04DD.f"
		    dswap_(&i__1, &qg[*ilo + (*ilo + 1) * qg_dim1], ldqg, &qg[
			    *ilo + (i__ + 1) * qg_dim1], &c__1);
#line 238 "MB04DD.f"
		}
#line 239 "MB04DD.f"
		++(*ilo);
#line 240 "MB04DD.f"
	    }
/*           END WHILE 30 */

/*           Scan columns N+ILO .. 2*N. */

#line 245 "MB04DD.f"
	    i__ = *ilo;
/*           WHILE ( I.LE.N .AND. ILO.EQ.ILOOLD ) */
#line 247 "MB04DD.f"
L80:
#line 247 "MB04DD.f"
	    if (i__ <= *n && *ilo == iloold) {
#line 248 "MB04DD.f"
		i__1 = i__ - 1;
#line 248 "MB04DD.f"
		for (j = *ilo; j <= i__1; ++j) {
#line 249 "MB04DD.f"
		    if (a[i__ + j * a_dim1] != 0.) {
#line 250 "MB04DD.f"
			++i__;
#line 251 "MB04DD.f"
			goto L80;
#line 252 "MB04DD.f"
		    }
#line 253 "MB04DD.f"
/* L90: */
#line 253 "MB04DD.f"
		}
#line 254 "MB04DD.f"
		i__1 = *n;
#line 254 "MB04DD.f"
		for (j = i__ + 1; j <= i__1; ++j) {
#line 255 "MB04DD.f"
		    if (a[i__ + j * a_dim1] != 0.) {
#line 256 "MB04DD.f"
			++i__;
#line 257 "MB04DD.f"
			goto L80;
#line 258 "MB04DD.f"
		    }
#line 259 "MB04DD.f"
/* L100: */
#line 259 "MB04DD.f"
		}
#line 260 "MB04DD.f"
		i__1 = i__;
#line 260 "MB04DD.f"
		for (j = *ilo; j <= i__1; ++j) {
#line 261 "MB04DD.f"
		    if (qg[j + (i__ + 1) * qg_dim1] != 0.) {
#line 262 "MB04DD.f"
			++i__;
#line 263 "MB04DD.f"
			goto L80;
#line 264 "MB04DD.f"
		    }
#line 265 "MB04DD.f"
/* L110: */
#line 265 "MB04DD.f"
		}
#line 266 "MB04DD.f"
		i__1 = *n;
#line 266 "MB04DD.f"
		for (j = i__ + 1; j <= i__1; ++j) {
#line 267 "MB04DD.f"
		    if (qg[i__ + (j + 1) * qg_dim1] != 0.) {
#line 268 "MB04DD.f"
			++i__;
#line 269 "MB04DD.f"
			goto L80;
#line 270 "MB04DD.f"
		    }
#line 271 "MB04DD.f"
/* L120: */
#line 271 "MB04DD.f"
		}
#line 272 "MB04DD.f"
		scale[*ilo] = (doublereal) (*n + i__);

/*              Exchange columns/rows I <-> I+N with a symplectic */
/*              generalized permutation. */

#line 277 "MB04DD.f"
		i__1 = i__ - *ilo;
#line 277 "MB04DD.f"
		dswap_(&i__1, &a[i__ + *ilo * a_dim1], lda, &qg[i__ + *ilo * 
			qg_dim1], ldqg);
#line 278 "MB04DD.f"
		i__1 = i__ - *ilo;
#line 278 "MB04DD.f"
		dscal_(&i__1, &c_b29, &a[i__ + *ilo * a_dim1], lda);
#line 279 "MB04DD.f"
		i__1 = *n - i__;
#line 279 "MB04DD.f"
		dswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &qg[i__ + 1 
			+ i__ * qg_dim1], &c__1);
#line 280 "MB04DD.f"
		i__1 = *n - i__;
#line 280 "MB04DD.f"
		dscal_(&i__1, &c_b29, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 281 "MB04DD.f"
		i__1 = i__ - 1;
#line 281 "MB04DD.f"
		dswap_(&i__1, &a[i__ * a_dim1 + 1], &c__1, &qg[(i__ + 1) * 
			qg_dim1 + 1], &c__1);
#line 282 "MB04DD.f"
		i__1 = i__ - 1;
#line 282 "MB04DD.f"
		dscal_(&i__1, &c_b29, &a[i__ * a_dim1 + 1], &c__1);
#line 283 "MB04DD.f"
		i__1 = *n - i__;
#line 283 "MB04DD.f"
		dswap_(&i__1, &a[i__ + 1 + i__ * a_dim1], &c__1, &qg[i__ + (
			i__ + 2) * qg_dim1], ldqg);
#line 284 "MB04DD.f"
		i__1 = *n - i__;
#line 284 "MB04DD.f"
		dscal_(&i__1, &c_b29, &a[i__ + 1 + i__ * a_dim1], &c__1);
#line 285 "MB04DD.f"
		a[i__ + i__ * a_dim1] = -a[i__ + i__ * a_dim1];
#line 286 "MB04DD.f"
		temp = qg[i__ + i__ * qg_dim1];
#line 287 "MB04DD.f"
		qg[i__ + i__ * qg_dim1] = -qg[i__ + (i__ + 1) * qg_dim1];
#line 288 "MB04DD.f"
		qg[i__ + (i__ + 1) * qg_dim1] = -temp;

/*              Exchange columns/rows ILO <-> I. */

#line 292 "MB04DD.f"
		if (*ilo != i__) {

#line 294 "MB04DD.f"
		    dswap_(n, &a[*ilo * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 
			    1], &c__1);
#line 295 "MB04DD.f"
		    i__1 = *n - *ilo + 1;
#line 295 "MB04DD.f"
		    dswap_(&i__1, &a[*ilo + *ilo * a_dim1], lda, &a[i__ + *
			    ilo * a_dim1], lda);

#line 297 "MB04DD.f"
		    dswap_(&c__1, &qg[i__ + *ilo * qg_dim1], ldqg, &qg[*ilo + 
			    *ilo * qg_dim1], ldqg);
#line 298 "MB04DD.f"
		    i__1 = *n - i__ + 1;
#line 298 "MB04DD.f"
		    dswap_(&i__1, &qg[i__ + i__ * qg_dim1], &c__1, &qg[i__ + *
			    ilo * qg_dim1], &c__1);
#line 299 "MB04DD.f"
		    i__1 = i__ - *ilo;
#line 299 "MB04DD.f"
		    dswap_(&i__1, &qg[*ilo + *ilo * qg_dim1], &c__1, &qg[i__ 
			    + *ilo * qg_dim1], ldqg);

#line 301 "MB04DD.f"
		    dswap_(ilo, &qg[(i__ + 1) * qg_dim1 + 1], &c__1, &qg[(*
			    ilo + 1) * qg_dim1 + 1], &c__1);
#line 302 "MB04DD.f"
		    i__1 = *n - i__ + 1;
#line 302 "MB04DD.f"
		    dswap_(&i__1, &qg[i__ + (i__ + 1) * qg_dim1], ldqg, &qg[*
			    ilo + (i__ + 1) * qg_dim1], ldqg);
#line 304 "MB04DD.f"
		    i__1 = i__ - *ilo;
#line 304 "MB04DD.f"
		    dswap_(&i__1, &qg[*ilo + (*ilo + 1) * qg_dim1], ldqg, &qg[
			    *ilo + (i__ + 1) * qg_dim1], &c__1);
#line 306 "MB04DD.f"
		}
#line 307 "MB04DD.f"
		++(*ilo);
#line 308 "MB04DD.f"
	    }
/*           END WHILE 80 */
#line 310 "MB04DD.f"
	    goto L20;
#line 311 "MB04DD.f"
	}
/*        END WHILE 20 */
#line 313 "MB04DD.f"
    }

#line 315 "MB04DD.f"
    i__1 = *n;
#line 315 "MB04DD.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 316 "MB04DD.f"
	scale[i__] = 1.;
#line 317 "MB04DD.f"
/* L130: */
#line 317 "MB04DD.f"
    }

/*     Scale to reduce the 1-norm of the remaining blocks. */

#line 321 "MB04DD.f"
    if (lscal) {
#line 322 "MB04DD.f"
	sclfac = dlamch_("B", (ftnlen)1);
#line 323 "MB04DD.f"
	sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 324 "MB04DD.f"
	sfmax1 = 1. / sfmin1;
#line 325 "MB04DD.f"
	sfmin2 = sfmin1 * sclfac;
#line 326 "MB04DD.f"
	sfmax2 = 1. / sfmin2;

/*        Scale the rows and columns one at a time to minimize the */
/*        1-norm of the remaining Hamiltonian submatrix. */
/*        Stop when the 1-norm is very roughly minimal. */

#line 332 "MB04DD.f"
L140:
#line 333 "MB04DD.f"
	conv = TRUE_;
#line 334 "MB04DD.f"
	i__1 = *n;
#line 334 "MB04DD.f"
	for (i__ = *ilo; i__ <= i__1; ++i__) {

/*              Compute 1-norm of row and column I without diagonal */
/*              elements. */

#line 339 "MB04DD.f"
	    i__2 = i__ - *ilo;
#line 339 "MB04DD.f"
	    i__3 = *n - i__;
#line 339 "MB04DD.f"
	    i__4 = i__ - *ilo;
#line 339 "MB04DD.f"
	    i__5 = *n - i__;
#line 339 "MB04DD.f"
	    r__ = dasum_(&i__2, &a[i__ + *ilo * a_dim1], lda) + dasum_(&i__3, 
		    &a[i__ + (i__ + 1) * a_dim1], lda) + dasum_(&i__4, &qg[*
		    ilo + (i__ + 1) * qg_dim1], &c__1) + dasum_(&i__5, &qg[
		    i__ + (i__ + 2) * qg_dim1], ldqg);
#line 343 "MB04DD.f"
	    i__2 = i__ - *ilo;
#line 343 "MB04DD.f"
	    i__3 = *n - i__;
#line 343 "MB04DD.f"
	    i__4 = i__ - *ilo;
#line 343 "MB04DD.f"
	    i__5 = *n - i__;
#line 343 "MB04DD.f"
	    c__ = dasum_(&i__2, &a[*ilo + i__ * a_dim1], &c__1) + dasum_(&
		    i__3, &a[i__ + 1 + i__ * a_dim1], &c__1) + dasum_(&i__4, &
		    qg[i__ + *ilo * qg_dim1], ldqg) + dasum_(&i__5, &qg[i__ + 
		    1 + i__ * qg_dim1], &c__1);
#line 347 "MB04DD.f"
	    qii = (d__1 = qg[i__ + i__ * qg_dim1], abs(d__1));
#line 348 "MB04DD.f"
	    gii = (d__1 = qg[i__ + (i__ + 1) * qg_dim1], abs(d__1));

/*              Compute inf-norms of row and column I. */

#line 352 "MB04DD.f"
	    i__2 = *n - *ilo + 1;
#line 352 "MB04DD.f"
	    ic = idamax_(&i__2, &a[i__ + *ilo * a_dim1], lda);
#line 353 "MB04DD.f"
	    maxr = (d__1 = a[i__ + (ic + *ilo - 1) * a_dim1], abs(d__1));
#line 354 "MB04DD.f"
	    if (i__ > 1) {
#line 355 "MB04DD.f"
		i__2 = i__ - 1;
#line 355 "MB04DD.f"
		ic = idamax_(&i__2, &qg[(i__ + 1) * qg_dim1 + 1], &c__1);
/* Computing MAX */
#line 356 "MB04DD.f"
		d__2 = maxr, d__3 = (d__1 = qg[ic + (i__ + 1) * qg_dim1], abs(
			d__1));
#line 356 "MB04DD.f"
		maxr = max(d__2,d__3);
#line 357 "MB04DD.f"
	    }
#line 358 "MB04DD.f"
	    if (*n > i__) {
#line 359 "MB04DD.f"
		i__2 = *n - i__;
#line 359 "MB04DD.f"
		ic = idamax_(&i__2, &qg[i__ + (i__ + 2) * qg_dim1], ldqg);
/* Computing MAX */
#line 360 "MB04DD.f"
		d__2 = maxr, d__3 = (d__1 = qg[i__ + (ic + i__ + 1) * qg_dim1]
			, abs(d__1));
#line 360 "MB04DD.f"
		maxr = max(d__2,d__3);
#line 361 "MB04DD.f"
	    }
#line 362 "MB04DD.f"
	    ic = idamax_(n, &a[i__ * a_dim1 + 1], &c__1);
#line 363 "MB04DD.f"
	    maxc = (d__1 = a[ic + i__ * a_dim1], abs(d__1));
#line 364 "MB04DD.f"
	    if (i__ > *ilo) {
#line 365 "MB04DD.f"
		i__2 = i__ - *ilo;
#line 365 "MB04DD.f"
		ic = idamax_(&i__2, &qg[i__ + *ilo * qg_dim1], ldqg);
/* Computing MAX */
#line 366 "MB04DD.f"
		d__2 = maxc, d__3 = (d__1 = qg[i__ + (ic + *ilo - 1) * 
			qg_dim1], abs(d__1));
#line 366 "MB04DD.f"
		maxc = max(d__2,d__3);
#line 367 "MB04DD.f"
	    }
#line 368 "MB04DD.f"
	    if (*n > i__) {
#line 369 "MB04DD.f"
		i__2 = *n - i__;
#line 369 "MB04DD.f"
		ic = idamax_(&i__2, &qg[i__ + 1 + i__ * qg_dim1], &c__1);
/* Computing MAX */
#line 370 "MB04DD.f"
		d__2 = maxc, d__3 = (d__1 = qg[ic + i__ + i__ * qg_dim1], abs(
			d__1));
#line 370 "MB04DD.f"
		maxc = max(d__2,d__3);
#line 371 "MB04DD.f"
	    }
#line 372 "MB04DD.f"
	    if (c__ + qii == 0. || r__ + gii == 0.) {
#line 372 "MB04DD.f"
		goto L170;
#line 372 "MB04DD.f"
	    }

#line 375 "MB04DD.f"
	    f = 1.;
#line 376 "MB04DD.f"
L150:
/* Computing MAX */
#line 377 "MB04DD.f"
	    d__1 = f * sclfac, d__2 = c__ * sclfac, d__1 = max(d__1,d__2), 
		    d__2 = maxc * sclfac, d__1 = max(d__1,d__2), d__2 = qii * 
		    sclfac * sclfac;
/* Computing MIN */
/* Computing MAX */
#line 377 "MB04DD.f"
	    d__5 = maxr / sclfac, d__6 = gii / sclfac / sclfac;
#line 377 "MB04DD.f"
	    d__3 = (r__ + gii / sclfac) / sclfac, d__4 = max(d__5,d__6);
#line 377 "MB04DD.f"
	    if ((r__ + gii / sclfac) / sclfac >= (c__ + qii * sclfac) * 
		    sclfac && max(d__1,d__2) < sfmax2 && min(d__3,d__4) > 
		    sfmin2) {
#line 383 "MB04DD.f"
		f *= sclfac;
#line 384 "MB04DD.f"
		c__ *= sclfac;
#line 385 "MB04DD.f"
		qii = qii * sclfac * sclfac;
#line 386 "MB04DD.f"
		r__ /= sclfac;
#line 387 "MB04DD.f"
		gii = gii / sclfac / sclfac;
#line 388 "MB04DD.f"
		maxc *= sclfac;
#line 389 "MB04DD.f"
		maxr /= sclfac;
#line 390 "MB04DD.f"
		goto L150;
#line 391 "MB04DD.f"
	    }

#line 393 "MB04DD.f"
L160:
/* Computing MAX */
#line 394 "MB04DD.f"
	    d__1 = r__ * sclfac, d__2 = maxr * sclfac, d__1 = max(d__1,d__2), 
		    d__2 = gii * sclfac * sclfac;
/* Computing MIN */
/* Computing MAX */
#line 394 "MB04DD.f"
	    d__5 = maxc / sclfac, d__6 = qii / sclfac / sclfac;
#line 394 "MB04DD.f"
	    d__3 = f / sclfac, d__4 = (c__ + qii / sclfac) / sclfac, d__3 = 
		    min(d__3,d__4), d__4 = max(d__5,d__6);
#line 394 "MB04DD.f"
	    if ((r__ + gii * sclfac) * sclfac <= (c__ + qii / sclfac) / 
		    sclfac && max(d__1,d__2) < sfmax2 && min(d__3,d__4) > 
		    sfmin2) {
#line 401 "MB04DD.f"
		f /= sclfac;
#line 402 "MB04DD.f"
		c__ /= sclfac;
#line 403 "MB04DD.f"
		qii = qii / sclfac / sclfac;
#line 404 "MB04DD.f"
		r__ *= sclfac;
#line 405 "MB04DD.f"
		gii = gii * sclfac * sclfac;
#line 406 "MB04DD.f"
		maxc /= sclfac;
#line 407 "MB04DD.f"
		maxr *= sclfac;
#line 408 "MB04DD.f"
		goto L160;
#line 409 "MB04DD.f"
	    }

/*              Now balance if necessary. */

#line 413 "MB04DD.f"
	    if (f != 1.) {
#line 414 "MB04DD.f"
		if (f < 1. && scale[i__] < 1.) {
#line 415 "MB04DD.f"
		    if (f * scale[i__] <= sfmin1) {
#line 415 "MB04DD.f"
			goto L170;
#line 415 "MB04DD.f"
		    }
#line 417 "MB04DD.f"
		}
#line 418 "MB04DD.f"
		if (f > 1. && scale[i__] > 1.) {
#line 419 "MB04DD.f"
		    if (scale[i__] >= sfmax1 / f) {
#line 419 "MB04DD.f"
			goto L170;
#line 419 "MB04DD.f"
		    }
#line 421 "MB04DD.f"
		}
#line 422 "MB04DD.f"
		conv = FALSE_;
#line 423 "MB04DD.f"
		scale[i__] *= f;
#line 424 "MB04DD.f"
		i__2 = i__ - *ilo;
#line 424 "MB04DD.f"
		drscl_(&i__2, &f, &a[i__ + *ilo * a_dim1], lda);
#line 425 "MB04DD.f"
		i__2 = *n - i__;
#line 425 "MB04DD.f"
		drscl_(&i__2, &f, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 426 "MB04DD.f"
		i__2 = i__ - 1;
#line 426 "MB04DD.f"
		dscal_(&i__2, &f, &a[i__ * a_dim1 + 1], &c__1);
#line 427 "MB04DD.f"
		i__2 = *n - i__;
#line 427 "MB04DD.f"
		dscal_(&i__2, &f, &a[i__ + 1 + i__ * a_dim1], &c__1);
#line 428 "MB04DD.f"
		i__2 = i__ - 1;
#line 428 "MB04DD.f"
		drscl_(&i__2, &f, &qg[(i__ + 1) * qg_dim1 + 1], &c__1);
#line 429 "MB04DD.f"
		qg[i__ + (i__ + 1) * qg_dim1] = qg[i__ + (i__ + 1) * qg_dim1] 
			/ f / f;
#line 430 "MB04DD.f"
		i__2 = *n - i__;
#line 430 "MB04DD.f"
		drscl_(&i__2, &f, &qg[i__ + (i__ + 2) * qg_dim1], ldqg);
#line 431 "MB04DD.f"
		i__2 = i__ - *ilo;
#line 431 "MB04DD.f"
		dscal_(&i__2, &f, &qg[i__ + *ilo * qg_dim1], ldqg);
#line 432 "MB04DD.f"
		qg[i__ + i__ * qg_dim1] = qg[i__ + i__ * qg_dim1] * f * f;
#line 433 "MB04DD.f"
		i__2 = *n - i__;
#line 433 "MB04DD.f"
		dscal_(&i__2, &f, &qg[i__ + 1 + i__ * qg_dim1], &c__1);
#line 434 "MB04DD.f"
	    }
#line 435 "MB04DD.f"
L170:
#line 435 "MB04DD.f"
	    ;
#line 435 "MB04DD.f"
	}
#line 436 "MB04DD.f"
	if (! conv) {
#line 436 "MB04DD.f"
	    goto L140;
#line 436 "MB04DD.f"
	}
#line 437 "MB04DD.f"
    }
#line 438 "MB04DD.f"
    return 0;
/* *** Last line of MB04DD *** */
} /* mb04dd_ */

