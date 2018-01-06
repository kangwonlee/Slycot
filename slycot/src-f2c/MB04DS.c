#line 1 "MB04DS.f"
/* MB04DS.f -- translated by f2c (version 20100827).
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

#line 1 "MB04DS.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b19 = -1.;

/* Subroutine */ int mb04ds_(char *job, integer *n, doublereal *a, integer *
	lda, doublereal *qg, integer *ldqg, integer *ilo, doublereal *scale, 
	integer *info, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2, i__3, i__4, 
	    i__5;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal c__, f, g;
    static integer i__, j;
    static doublereal r__, s;
    static integer ic;
    static doublereal maxc;
    static logical conv;
    static doublereal maxr;
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

/*     To balance a real skew-Hamiltonian matrix */

/*                   [  A   G  ] */
/*              S =  [       T ] , */
/*                   [  Q   A  ] */

/*     where A is an N-by-N matrix and G, Q are N-by-N skew-symmetric */
/*     matrices. This involves, first, permuting S by a symplectic */
/*     similarity transformation to isolate eigenvalues in the first */
/*     1:ILO-1 elements on the diagonal of A; and second, applying a */
/*     diagonal similarity transformation to rows and columns */
/*     ILO:2*N-ILO+1 to make the rows and columns as close in 1-norm */
/*     as possible. Both steps are optional. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies the operations to be performed on S: */
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
/*             the matrix A of the balanced skew-Hamiltonian. In */
/*             particular, the lower triangular part of the first ILO-1 */
/*             columns of A is zero. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     QG      (input/output) DOUBLE PRECISION array, dimension */
/*                            (LDQG,N) */
/*             On entry, the leading N-by-N+1 part of this array must */
/*             contain in columns 1:N the strictly lower triangular part */
/*             of the matrix Q and in columns 2:N+1 the strictly upper */
/*             triangular part of the matrix G. The parts containing the */
/*             diagonal and the first supdiagonal of this array are not */
/*             referenced. */
/*             On exit, the leading N-by-N+1 part of this array contains */
/*             the strictly lower and strictly upper triangular parts of */
/*             the matrices Q and G, respectively, of the balanced */
/*             skew-Hamiltonian. In particular, the strictly lower */
/*             triangular part of the first ILO-1 columns of QG is zero. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= MAX(1,N). */

/*     ILO     (output) INTEGER */
/*             ILO-1 is the number of deflated eigenvalues in the */
/*             balanced skew-Hamiltonian matrix. */

/*     SCALE   (output) DOUBLE PRECISION array of dimension (N) */
/*             Details of the permutations and scaling factors applied to */
/*             S.  For j = 1,...,ILO-1 let P(j) = SCALE(j). If P(j) <= N, */
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

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DSHBAL). */

/*     KEYWORDS */

/*     Balancing, skew-Hamiltonian matrix. */

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

#line 153 "MB04DS.f"
    /* Parameter adjustments */
#line 153 "MB04DS.f"
    a_dim1 = *lda;
#line 153 "MB04DS.f"
    a_offset = 1 + a_dim1;
#line 153 "MB04DS.f"
    a -= a_offset;
#line 153 "MB04DS.f"
    qg_dim1 = *ldqg;
#line 153 "MB04DS.f"
    qg_offset = 1 + qg_dim1;
#line 153 "MB04DS.f"
    qg -= qg_offset;
#line 153 "MB04DS.f"
    --scale;
#line 153 "MB04DS.f"

#line 153 "MB04DS.f"
    /* Function Body */
#line 153 "MB04DS.f"
    *info = 0;
#line 154 "MB04DS.f"
    lperm = lsame_(job, "P", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);
#line 155 "MB04DS.f"
    lscal = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "B", (
	    ftnlen)1, (ftnlen)1);

#line 157 "MB04DS.f"
    if (! lperm && ! lscal && ! lsame_(job, "N", (ftnlen)1, (ftnlen)1)) {
#line 159 "MB04DS.f"
	*info = -1;
#line 160 "MB04DS.f"
    } else if (*n < 0) {
#line 161 "MB04DS.f"
	*info = -2;
#line 162 "MB04DS.f"
    } else if (*lda < max(1,*n)) {
#line 163 "MB04DS.f"
	*info = -4;
#line 164 "MB04DS.f"
    } else if (*ldqg < max(1,*n)) {
#line 165 "MB04DS.f"
	*info = -6;
#line 166 "MB04DS.f"
    }

/*     Return if there were illegal values. */

#line 170 "MB04DS.f"
    if (*info != 0) {
#line 171 "MB04DS.f"
	i__1 = -(*info);
#line 171 "MB04DS.f"
	xerbla_("MB04DS", &i__1, (ftnlen)6);
#line 172 "MB04DS.f"
	return 0;
#line 173 "MB04DS.f"
    }

#line 175 "MB04DS.f"
    *ilo = 1;

/*     Quick return if possible. */

#line 179 "MB04DS.f"
    if (*n == 0) {
#line 179 "MB04DS.f"
	return 0;
#line 179 "MB04DS.f"
    }
#line 181 "MB04DS.f"
    if (! lperm && ! lscal) {
#line 182 "MB04DS.f"
	i__1 = *n;
#line 182 "MB04DS.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 183 "MB04DS.f"
	    scale[i__] = 1.;
#line 184 "MB04DS.f"
/* L10: */
#line 184 "MB04DS.f"
	}
#line 185 "MB04DS.f"
	return 0;
#line 186 "MB04DS.f"
    }

/*     Permutations to isolate eigenvalues if possible. */

#line 190 "MB04DS.f"
    if (lperm) {
#line 191 "MB04DS.f"
	iloold = 0;
/*        WHILE ( ILO.NE.ILOOLD ) */
#line 193 "MB04DS.f"
L20:
#line 193 "MB04DS.f"
	if (*ilo != iloold) {
#line 194 "MB04DS.f"
	    iloold = *ilo;

/*           Scan columns ILO .. N. */

#line 198 "MB04DS.f"
	    i__ = *ilo;
/*           WHILE ( I.LE.N .AND. ILO.EQ.ILOOLD ) */
#line 200 "MB04DS.f"
L30:
#line 200 "MB04DS.f"
	    if (i__ <= *n && *ilo == iloold) {
#line 201 "MB04DS.f"
		i__1 = i__ - 1;
#line 201 "MB04DS.f"
		for (j = *ilo; j <= i__1; ++j) {
#line 202 "MB04DS.f"
		    if (a[j + i__ * a_dim1] != 0.) {
#line 203 "MB04DS.f"
			++i__;
#line 204 "MB04DS.f"
			goto L30;
#line 205 "MB04DS.f"
		    }
#line 206 "MB04DS.f"
/* L40: */
#line 206 "MB04DS.f"
		}
#line 207 "MB04DS.f"
		i__1 = *n;
#line 207 "MB04DS.f"
		for (j = i__ + 1; j <= i__1; ++j) {
#line 208 "MB04DS.f"
		    if (a[j + i__ * a_dim1] != 0.) {
#line 209 "MB04DS.f"
			++i__;
#line 210 "MB04DS.f"
			goto L30;
#line 211 "MB04DS.f"
		    }
#line 212 "MB04DS.f"
/* L50: */
#line 212 "MB04DS.f"
		}
#line 213 "MB04DS.f"
		i__1 = i__ - 1;
#line 213 "MB04DS.f"
		for (j = *ilo; j <= i__1; ++j) {
#line 214 "MB04DS.f"
		    if (qg[i__ + j * qg_dim1] != 0.) {
#line 215 "MB04DS.f"
			++i__;
#line 216 "MB04DS.f"
			goto L30;
#line 217 "MB04DS.f"
		    }
#line 218 "MB04DS.f"
/* L60: */
#line 218 "MB04DS.f"
		}
#line 219 "MB04DS.f"
		i__1 = *n;
#line 219 "MB04DS.f"
		for (j = i__ + 1; j <= i__1; ++j) {
#line 220 "MB04DS.f"
		    if (qg[j + i__ * qg_dim1] != 0.) {
#line 221 "MB04DS.f"
			++i__;
#line 222 "MB04DS.f"
			goto L30;
#line 223 "MB04DS.f"
		    }
#line 224 "MB04DS.f"
/* L70: */
#line 224 "MB04DS.f"
		}

/*              Exchange columns/rows ILO <-> I. */

#line 228 "MB04DS.f"
		scale[*ilo] = (doublereal) i__;
#line 229 "MB04DS.f"
		if (*ilo != i__) {

#line 231 "MB04DS.f"
		    dswap_(n, &a[*ilo * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 
			    1], &c__1);
#line 232 "MB04DS.f"
		    i__1 = *n - *ilo + 1;
#line 232 "MB04DS.f"
		    dswap_(&i__1, &a[*ilo + *ilo * a_dim1], lda, &a[i__ + *
			    ilo * a_dim1], lda);

#line 234 "MB04DS.f"
		    if (i__ < *n) {
#line 234 "MB04DS.f"
			i__1 = *n - i__;
#line 234 "MB04DS.f"
			dswap_(&i__1, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &
				qg[i__ + 1 + *ilo * qg_dim1], &c__1);
#line 234 "MB04DS.f"
		    }
#line 236 "MB04DS.f"
		    if (i__ > *ilo + 1) {
#line 237 "MB04DS.f"
			i__1 = i__ - *ilo - 1;
#line 237 "MB04DS.f"
			dscal_(&i__1, &c_b19, &qg[*ilo + 1 + *ilo * qg_dim1], 
				&c__1);
#line 238 "MB04DS.f"
			i__1 = i__ - *ilo - 1;
#line 238 "MB04DS.f"
			dswap_(&i__1, &qg[*ilo + 1 + *ilo * qg_dim1], &c__1, &
				qg[i__ + (*ilo + 1) * qg_dim1], ldqg);
#line 240 "MB04DS.f"
		    }

#line 242 "MB04DS.f"
		    i__1 = *ilo - 1;
#line 242 "MB04DS.f"
		    dswap_(&i__1, &qg[(i__ + 1) * qg_dim1 + 1], &c__1, &qg[(*
			    ilo + 1) * qg_dim1 + 1], &c__1);
#line 243 "MB04DS.f"
		    if (*n > i__) {
#line 243 "MB04DS.f"
			i__1 = *n - i__;
#line 243 "MB04DS.f"
			dswap_(&i__1, &qg[i__ + (i__ + 2) * qg_dim1], ldqg, &
				qg[*ilo + (i__ + 2) * qg_dim1], ldqg);
#line 243 "MB04DS.f"
		    }
#line 246 "MB04DS.f"
		    if (i__ > *ilo + 1) {
#line 247 "MB04DS.f"
			i__1 = i__ - *ilo - 1;
#line 247 "MB04DS.f"
			dscal_(&i__1, &c_b19, &qg[*ilo + 1 + (i__ + 1) * 
				qg_dim1], &c__1);
#line 248 "MB04DS.f"
			i__1 = i__ - *ilo - 1;
#line 248 "MB04DS.f"
			dswap_(&i__1, &qg[*ilo + (*ilo + 2) * qg_dim1], ldqg, 
				&qg[*ilo + 1 + (i__ + 1) * qg_dim1], &c__1);
#line 250 "MB04DS.f"
		    }
#line 251 "MB04DS.f"
		    i__1 = i__ - *ilo;
#line 251 "MB04DS.f"
		    dscal_(&i__1, &c_b19, &qg[*ilo + (i__ + 1) * qg_dim1], &
			    c__1);
#line 252 "MB04DS.f"
		}
#line 253 "MB04DS.f"
		++(*ilo);
#line 254 "MB04DS.f"
	    }
/*           END WHILE 30 */

/*           Scan columns N+ILO .. 2*N. */

#line 259 "MB04DS.f"
	    i__ = *ilo;
/*           WHILE ( I.LE.N .AND. ILO.EQ.ILOOLD ) */
#line 261 "MB04DS.f"
L80:
#line 261 "MB04DS.f"
	    if (i__ <= *n && *ilo == iloold) {
#line 262 "MB04DS.f"
		i__1 = i__ - 1;
#line 262 "MB04DS.f"
		for (j = *ilo; j <= i__1; ++j) {
#line 263 "MB04DS.f"
		    if (a[i__ + j * a_dim1] != 0.) {
#line 264 "MB04DS.f"
			++i__;
#line 265 "MB04DS.f"
			goto L80;
#line 266 "MB04DS.f"
		    }
#line 267 "MB04DS.f"
/* L90: */
#line 267 "MB04DS.f"
		}
#line 268 "MB04DS.f"
		i__1 = *n;
#line 268 "MB04DS.f"
		for (j = i__ + 1; j <= i__1; ++j) {
#line 269 "MB04DS.f"
		    if (a[i__ + j * a_dim1] != 0.) {
#line 270 "MB04DS.f"
			++i__;
#line 271 "MB04DS.f"
			goto L80;
#line 272 "MB04DS.f"
		    }
#line 273 "MB04DS.f"
/* L100: */
#line 273 "MB04DS.f"
		}
#line 274 "MB04DS.f"
		i__1 = i__ - 1;
#line 274 "MB04DS.f"
		for (j = *ilo; j <= i__1; ++j) {
#line 275 "MB04DS.f"
		    if (qg[j + (i__ + 1) * qg_dim1] != 0.) {
#line 276 "MB04DS.f"
			++i__;
#line 277 "MB04DS.f"
			goto L80;
#line 278 "MB04DS.f"
		    }
#line 279 "MB04DS.f"
/* L110: */
#line 279 "MB04DS.f"
		}
#line 280 "MB04DS.f"
		i__1 = *n;
#line 280 "MB04DS.f"
		for (j = i__ + 1; j <= i__1; ++j) {
#line 281 "MB04DS.f"
		    if (qg[i__ + (j + 1) * qg_dim1] != 0.) {
#line 282 "MB04DS.f"
			++i__;
#line 283 "MB04DS.f"
			goto L80;
#line 284 "MB04DS.f"
		    }
#line 285 "MB04DS.f"
/* L120: */
#line 285 "MB04DS.f"
		}
#line 286 "MB04DS.f"
		scale[*ilo] = (doublereal) (*n + i__);

/*              Exchange columns/rows I <-> I+N with a symplectic */
/*              generalized permutation. */

#line 291 "MB04DS.f"
		i__1 = i__ - *ilo;
#line 291 "MB04DS.f"
		dswap_(&i__1, &a[i__ + *ilo * a_dim1], lda, &qg[i__ + *ilo * 
			qg_dim1], ldqg);
#line 292 "MB04DS.f"
		i__1 = i__ - *ilo;
#line 292 "MB04DS.f"
		dscal_(&i__1, &c_b19, &a[i__ + *ilo * a_dim1], lda);
#line 293 "MB04DS.f"
		i__1 = *n - i__;
#line 293 "MB04DS.f"
		dswap_(&i__1, &a[i__ + (i__ + 1) * a_dim1], lda, &qg[i__ + 1 
			+ i__ * qg_dim1], &c__1);
#line 294 "MB04DS.f"
		i__1 = *n - i__;
#line 294 "MB04DS.f"
		dscal_(&i__1, &c_b19, &qg[i__ + 1 + i__ * qg_dim1], &c__1);
#line 295 "MB04DS.f"
		i__1 = i__ - 1;
#line 295 "MB04DS.f"
		dswap_(&i__1, &a[i__ * a_dim1 + 1], &c__1, &qg[(i__ + 1) * 
			qg_dim1 + 1], &c__1);
#line 296 "MB04DS.f"
		i__1 = i__ - 1;
#line 296 "MB04DS.f"
		dscal_(&i__1, &c_b19, &a[i__ * a_dim1 + 1], &c__1);
#line 297 "MB04DS.f"
		i__1 = *n - i__;
#line 297 "MB04DS.f"
		dscal_(&i__1, &c_b19, &a[i__ + 1 + i__ * a_dim1], &c__1);
#line 298 "MB04DS.f"
		i__1 = *n - i__;
#line 298 "MB04DS.f"
		dswap_(&i__1, &a[i__ + 1 + i__ * a_dim1], &c__1, &qg[i__ + (
			i__ + 2) * qg_dim1], ldqg);

/*              Exchange columns/rows ILO <-> I. */

#line 302 "MB04DS.f"
		if (*ilo != i__) {

#line 304 "MB04DS.f"
		    dswap_(n, &a[*ilo * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 
			    1], &c__1);
#line 305 "MB04DS.f"
		    i__1 = *n - *ilo + 1;
#line 305 "MB04DS.f"
		    dswap_(&i__1, &a[*ilo + *ilo * a_dim1], lda, &a[i__ + *
			    ilo * a_dim1], lda);

#line 307 "MB04DS.f"
		    if (i__ < *n) {
#line 307 "MB04DS.f"
			i__1 = *n - i__;
#line 307 "MB04DS.f"
			dswap_(&i__1, &qg[i__ + 1 + i__ * qg_dim1], &c__1, &
				qg[i__ + 1 + *ilo * qg_dim1], &c__1);
#line 307 "MB04DS.f"
		    }
#line 309 "MB04DS.f"
		    if (i__ > *ilo + 1) {
#line 310 "MB04DS.f"
			i__1 = i__ - *ilo - 1;
#line 310 "MB04DS.f"
			dscal_(&i__1, &c_b19, &qg[*ilo + 1 + *ilo * qg_dim1], 
				&c__1);
#line 311 "MB04DS.f"
			i__1 = i__ - *ilo - 1;
#line 311 "MB04DS.f"
			dswap_(&i__1, &qg[*ilo + 1 + *ilo * qg_dim1], &c__1, &
				qg[i__ + (*ilo + 1) * qg_dim1], ldqg);
#line 313 "MB04DS.f"
		    }

#line 315 "MB04DS.f"
		    i__1 = *ilo - 1;
#line 315 "MB04DS.f"
		    dswap_(&i__1, &qg[(i__ + 1) * qg_dim1 + 1], &c__1, &qg[(*
			    ilo + 1) * qg_dim1 + 1], &c__1);
#line 316 "MB04DS.f"
		    if (*n > i__) {
#line 316 "MB04DS.f"
			i__1 = *n - i__;
#line 316 "MB04DS.f"
			dswap_(&i__1, &qg[i__ + (i__ + 2) * qg_dim1], ldqg, &
				qg[*ilo + (i__ + 2) * qg_dim1], ldqg);
#line 316 "MB04DS.f"
		    }
#line 319 "MB04DS.f"
		    if (i__ > *ilo + 1) {
#line 320 "MB04DS.f"
			i__1 = i__ - *ilo - 1;
#line 320 "MB04DS.f"
			dscal_(&i__1, &c_b19, &qg[*ilo + 1 + (i__ + 1) * 
				qg_dim1], &c__1);
#line 321 "MB04DS.f"
			i__1 = i__ - *ilo - 1;
#line 321 "MB04DS.f"
			dswap_(&i__1, &qg[*ilo + (*ilo + 2) * qg_dim1], ldqg, 
				&qg[*ilo + 1 + (i__ + 1) * qg_dim1], &c__1);
#line 323 "MB04DS.f"
		    }
#line 324 "MB04DS.f"
		    i__1 = i__ - *ilo;
#line 324 "MB04DS.f"
		    dscal_(&i__1, &c_b19, &qg[*ilo + (i__ + 1) * qg_dim1], &
			    c__1);
#line 325 "MB04DS.f"
		}
#line 326 "MB04DS.f"
		++(*ilo);
#line 327 "MB04DS.f"
	    }
/*           END WHILE 80 */
#line 329 "MB04DS.f"
	    goto L20;
#line 330 "MB04DS.f"
	}
/*        END WHILE 20 */
#line 332 "MB04DS.f"
    }

#line 334 "MB04DS.f"
    i__1 = *n;
#line 334 "MB04DS.f"
    for (i__ = *ilo; i__ <= i__1; ++i__) {
#line 335 "MB04DS.f"
	scale[i__] = 1.;
#line 336 "MB04DS.f"
/* L130: */
#line 336 "MB04DS.f"
    }

/*     Scale to reduce the 1-norm of the remaining blocks. */

#line 340 "MB04DS.f"
    if (lscal) {
#line 341 "MB04DS.f"
	sclfac = dlamch_("B", (ftnlen)1);
#line 342 "MB04DS.f"
	sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 343 "MB04DS.f"
	sfmax1 = 1. / sfmin1;
#line 344 "MB04DS.f"
	sfmin2 = sfmin1 * sclfac;
#line 345 "MB04DS.f"
	sfmax2 = 1. / sfmin2;

/*        Scale the rows and columns one at a time to minimize the */
/*        1-norm of the skew-Hamiltonian submatrix. */
/*        Stop when the 1-norm is very roughly minimal. */

#line 351 "MB04DS.f"
L140:
#line 352 "MB04DS.f"
	conv = TRUE_;
#line 353 "MB04DS.f"
	i__1 = *n;
#line 353 "MB04DS.f"
	for (i__ = *ilo; i__ <= i__1; ++i__) {

/*              Compute 1-norm of row and column I without diagonal */
/*              elements. */

#line 358 "MB04DS.f"
	    i__2 = i__ - *ilo;
#line 358 "MB04DS.f"
	    i__3 = *n - i__;
#line 358 "MB04DS.f"
	    i__4 = i__ - *ilo;
#line 358 "MB04DS.f"
	    i__5 = *n - i__;
#line 358 "MB04DS.f"
	    r__ = dasum_(&i__2, &a[i__ + *ilo * a_dim1], lda) + dasum_(&i__3, 
		    &a[i__ + (i__ + 1) * a_dim1], lda) + dasum_(&i__4, &qg[*
		    ilo + (i__ + 1) * qg_dim1], &c__1) + dasum_(&i__5, &qg[
		    i__ + (i__ + 2) * qg_dim1], ldqg);
#line 362 "MB04DS.f"
	    i__2 = i__ - *ilo;
#line 362 "MB04DS.f"
	    i__3 = *n - i__;
#line 362 "MB04DS.f"
	    i__4 = i__ - *ilo;
#line 362 "MB04DS.f"
	    i__5 = *n - i__;
#line 362 "MB04DS.f"
	    c__ = dasum_(&i__2, &a[*ilo + i__ * a_dim1], &c__1) + dasum_(&
		    i__3, &a[i__ + 1 + i__ * a_dim1], &c__1) + dasum_(&i__4, &
		    qg[i__ + *ilo * qg_dim1], ldqg) + dasum_(&i__5, &qg[i__ + 
		    1 + i__ * qg_dim1], &c__1);

/*              Compute inf-norms of row and column I. */

#line 369 "MB04DS.f"
	    i__2 = *n - *ilo + 1;
#line 369 "MB04DS.f"
	    ic = idamax_(&i__2, &a[i__ + *ilo * a_dim1], lda);
#line 370 "MB04DS.f"
	    maxr = (d__1 = a[i__ + (ic + *ilo - 1) * a_dim1], abs(d__1));
#line 371 "MB04DS.f"
	    if (i__ > 1) {
#line 372 "MB04DS.f"
		i__2 = i__ - 1;
#line 372 "MB04DS.f"
		ic = idamax_(&i__2, &qg[(i__ + 1) * qg_dim1 + 1], &c__1);
/* Computing MAX */
#line 373 "MB04DS.f"
		d__2 = maxr, d__3 = (d__1 = qg[ic + (i__ + 1) * qg_dim1], abs(
			d__1));
#line 373 "MB04DS.f"
		maxr = max(d__2,d__3);
#line 374 "MB04DS.f"
	    }
#line 375 "MB04DS.f"
	    if (*n > i__) {
#line 376 "MB04DS.f"
		i__2 = *n - i__;
#line 376 "MB04DS.f"
		ic = idamax_(&i__2, &qg[i__ + (i__ + 2) * qg_dim1], ldqg);
/* Computing MAX */
#line 377 "MB04DS.f"
		d__2 = maxr, d__3 = (d__1 = qg[i__ + (ic + i__ + 1) * qg_dim1]
			, abs(d__1));
#line 377 "MB04DS.f"
		maxr = max(d__2,d__3);
#line 378 "MB04DS.f"
	    }
#line 379 "MB04DS.f"
	    ic = idamax_(n, &a[i__ * a_dim1 + 1], &c__1);
#line 380 "MB04DS.f"
	    maxc = (d__1 = a[ic + i__ * a_dim1], abs(d__1));
#line 381 "MB04DS.f"
	    if (i__ > *ilo) {
#line 382 "MB04DS.f"
		i__2 = i__ - *ilo;
#line 382 "MB04DS.f"
		ic = idamax_(&i__2, &qg[i__ + *ilo * qg_dim1], ldqg);
/* Computing MAX */
#line 383 "MB04DS.f"
		d__2 = maxc, d__3 = (d__1 = qg[i__ + (ic + *ilo - 1) * 
			qg_dim1], abs(d__1));
#line 383 "MB04DS.f"
		maxc = max(d__2,d__3);
#line 384 "MB04DS.f"
	    }
#line 385 "MB04DS.f"
	    if (*n > i__) {
#line 386 "MB04DS.f"
		i__2 = *n - i__;
#line 386 "MB04DS.f"
		ic = idamax_(&i__2, &qg[i__ + 1 + i__ * qg_dim1], &c__1);
/* Computing MAX */
#line 387 "MB04DS.f"
		d__2 = maxc, d__3 = (d__1 = qg[ic + i__ + i__ * qg_dim1], abs(
			d__1));
#line 387 "MB04DS.f"
		maxc = max(d__2,d__3);
#line 388 "MB04DS.f"
	    }

#line 390 "MB04DS.f"
	    if (c__ == 0. || r__ == 0.) {
#line 390 "MB04DS.f"
		goto L190;
#line 390 "MB04DS.f"
	    }
#line 392 "MB04DS.f"
	    g = r__ / sclfac;
#line 393 "MB04DS.f"
	    f = 1.;
#line 394 "MB04DS.f"
	    s = c__ + r__;
#line 395 "MB04DS.f"
L150:
/* Computing MAX */
#line 396 "MB04DS.f"
	    d__1 = max(f,c__);
/* Computing MIN */
#line 396 "MB04DS.f"
	    d__2 = min(r__,g);
#line 396 "MB04DS.f"
	    if (c__ >= g || max(d__1,maxc) >= sfmax2 || min(d__2,maxr) <= 
		    sfmin2) {
#line 396 "MB04DS.f"
		goto L160;
#line 396 "MB04DS.f"
	    }
#line 399 "MB04DS.f"
	    f *= sclfac;
#line 400 "MB04DS.f"
	    g /= sclfac;
#line 401 "MB04DS.f"
	    c__ *= sclfac;
#line 402 "MB04DS.f"
	    r__ /= sclfac;
#line 403 "MB04DS.f"
	    maxc *= sclfac;
#line 404 "MB04DS.f"
	    maxr /= sclfac;
#line 405 "MB04DS.f"
	    goto L150;

#line 407 "MB04DS.f"
L160:
#line 408 "MB04DS.f"
	    g = c__ / sclfac;
#line 409 "MB04DS.f"
L170:
/* Computing MIN */
#line 410 "MB04DS.f"
	    d__1 = min(f,c__), d__1 = min(d__1,g);
#line 410 "MB04DS.f"
	    if (g < r__ || max(r__,maxr) >= sfmax2 || min(d__1,maxc) <= 
		    sfmin2) {
#line 410 "MB04DS.f"
		goto L180;
#line 410 "MB04DS.f"
	    }
#line 413 "MB04DS.f"
	    f /= sclfac;
#line 414 "MB04DS.f"
	    g /= sclfac;
#line 415 "MB04DS.f"
	    c__ /= sclfac;
#line 416 "MB04DS.f"
	    r__ *= sclfac;
#line 417 "MB04DS.f"
	    maxc /= sclfac;
#line 418 "MB04DS.f"
	    maxr *= sclfac;
#line 419 "MB04DS.f"
	    goto L170;

#line 421 "MB04DS.f"
L180:

/*              Now balance if necessary. */

#line 425 "MB04DS.f"
	    if (c__ + r__ >= s * .95) {
#line 425 "MB04DS.f"
		goto L190;
#line 425 "MB04DS.f"
	    }
#line 427 "MB04DS.f"
	    if (f < 1. && scale[i__] < 1.) {
#line 428 "MB04DS.f"
		if (f * scale[i__] <= sfmin1) {
#line 428 "MB04DS.f"
		    goto L190;
#line 428 "MB04DS.f"
		}
#line 430 "MB04DS.f"
	    }
#line 431 "MB04DS.f"
	    if (f > 1. && scale[i__] > 1.) {
#line 432 "MB04DS.f"
		if (scale[i__] >= sfmax1 / f) {
#line 432 "MB04DS.f"
		    goto L190;
#line 432 "MB04DS.f"
		}
#line 434 "MB04DS.f"
	    }
#line 435 "MB04DS.f"
	    conv = FALSE_;
#line 436 "MB04DS.f"
	    scale[i__] *= f;
#line 437 "MB04DS.f"
	    i__2 = i__ - *ilo;
#line 437 "MB04DS.f"
	    drscl_(&i__2, &f, &a[i__ + *ilo * a_dim1], lda);
#line 438 "MB04DS.f"
	    i__2 = *n - i__;
#line 438 "MB04DS.f"
	    drscl_(&i__2, &f, &a[i__ + (i__ + 1) * a_dim1], lda);
#line 439 "MB04DS.f"
	    i__2 = i__ - 1;
#line 439 "MB04DS.f"
	    dscal_(&i__2, &f, &a[i__ * a_dim1 + 1], &c__1);
#line 440 "MB04DS.f"
	    i__2 = *n - i__;
#line 440 "MB04DS.f"
	    dscal_(&i__2, &f, &a[i__ + 1 + i__ * a_dim1], &c__1);
#line 441 "MB04DS.f"
	    i__2 = i__ - 1;
#line 441 "MB04DS.f"
	    drscl_(&i__2, &f, &qg[(i__ + 1) * qg_dim1 + 1], &c__1);
#line 442 "MB04DS.f"
	    i__2 = *n - i__;
#line 442 "MB04DS.f"
	    drscl_(&i__2, &f, &qg[i__ + (i__ + 2) * qg_dim1], ldqg);
#line 443 "MB04DS.f"
	    i__2 = i__ - *ilo;
#line 443 "MB04DS.f"
	    dscal_(&i__2, &f, &qg[i__ + *ilo * qg_dim1], ldqg);
#line 444 "MB04DS.f"
	    i__2 = *n - i__;
#line 444 "MB04DS.f"
	    dscal_(&i__2, &f, &qg[i__ + 1 + i__ * qg_dim1], &c__1);
#line 445 "MB04DS.f"
L190:
#line 445 "MB04DS.f"
	    ;
#line 445 "MB04DS.f"
	}
#line 446 "MB04DS.f"
	if (! conv) {
#line 446 "MB04DS.f"
	    goto L140;
#line 446 "MB04DS.f"
	}
#line 447 "MB04DS.f"
    }
#line 448 "MB04DS.f"
    return 0;
/* *** Last line of MB04DS *** */
} /* mb04ds_ */

