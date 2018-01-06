#line 1 "MB04MD.f"
/* MB04MD.f -- translated by f2c (version 20100827).
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

#line 1 "MB04MD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb04md_(integer *n, doublereal *maxred, doublereal *a, 
	integer *lda, doublereal *scale, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__, f, g;
    static integer i__, j;
    static doublereal r__, s, ca, ra;
    static integer ica, ira;
    static doublereal sred;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal anorm, sfmin1, sfmin2, sfmax1, sfmax2;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical noconv;
    static doublereal maxnrm;


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

/*     To reduce the 1-norm of a general real matrix A by balancing. */
/*     This involves diagonal similarity transformations applied */
/*     iteratively to A to make the rows and columns as close in norm as */
/*     possible. */

/*     This routine can be used instead LAPACK Library routine DGEBAL, */
/*     when no reduction of the 1-norm of the matrix is possible with */
/*     DGEBAL, as for upper triangular matrices. LAPACK Library routine */
/*     DGEBAK, with parameters ILO = 1, IHI = N, and JOB = 'S', should */
/*     be used to apply the backward transformation. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     MAXRED  (input/output) DOUBLE PRECISION */
/*             On entry, the maximum allowed reduction in the 1-norm of */
/*             A (in an iteration) if zero rows or columns are */
/*             encountered. */
/*             If MAXRED > 0.0, MAXRED must be larger than one (to enable */
/*             the norm reduction). */
/*             If MAXRED <= 0.0, then the value 10.0 for MAXRED is */
/*             used. */
/*             On exit, if the 1-norm of the given matrix A is non-zero, */
/*             the ratio between the 1-norm of the given matrix and the */
/*             1-norm of the balanced matrix. Usually, this ratio will be */
/*             larger than one, but it can sometimes be one, or even less */
/*             than one (for instance, for some companion matrices). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the input matrix A. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the balanced matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     SCALE   (output) DOUBLE PRECISION array, dimension (N) */
/*             The scaling factors applied to A.  If D(j) is the scaling */
/*             factor applied to row and column j, then SCALE(j) = D(j), */
/*             for j = 1,...,N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit. */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Balancing consists of applying a diagonal similarity */
/*     transformation inv(D) * A * D to make the 1-norms of each row */
/*     of A and its corresponding column nearly equal. */

/*     Information about the diagonal matrix D is returned in the vector */
/*     SCALE. */

/*     REFERENCES */

/*     [1] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routine MB04AD by T.W.C. Williams, */
/*     Kingston Polytechnic, United Kingdom, October 1984. */
/*     This subroutine is based on LAPACK routine DGEBAL, and routine */
/*     BALABC (A. Varga, German Aerospace Research Establishment, DLR). */


/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Balancing, eigenvalue, matrix algebra, matrix operations, */
/*     similarity transformation. */

/*  ********************************************************************* */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the scalar input arguments. */

#line 153 "MB04MD.f"
    /* Parameter adjustments */
#line 153 "MB04MD.f"
    a_dim1 = *lda;
#line 153 "MB04MD.f"
    a_offset = 1 + a_dim1;
#line 153 "MB04MD.f"
    a -= a_offset;
#line 153 "MB04MD.f"
    --scale;
#line 153 "MB04MD.f"

#line 153 "MB04MD.f"
    /* Function Body */
#line 153 "MB04MD.f"
    *info = 0;

#line 155 "MB04MD.f"
    if (*n < 0) {
#line 156 "MB04MD.f"
	*info = -1;
#line 157 "MB04MD.f"
    } else if (*maxred > 0. && *maxred < 1.) {
#line 158 "MB04MD.f"
	*info = -2;
#line 159 "MB04MD.f"
    } else if (*lda < max(1,*n)) {
#line 160 "MB04MD.f"
	*info = -4;
#line 161 "MB04MD.f"
    }
#line 162 "MB04MD.f"
    if (*info != 0) {
#line 163 "MB04MD.f"
	i__1 = -(*info);
#line 163 "MB04MD.f"
	xerbla_("MB04MD", &i__1, (ftnlen)6);
#line 164 "MB04MD.f"
	return 0;
#line 165 "MB04MD.f"
    }

#line 167 "MB04MD.f"
    if (*n == 0) {
#line 167 "MB04MD.f"
	return 0;
#line 167 "MB04MD.f"
    }

#line 170 "MB04MD.f"
    i__1 = *n;
#line 170 "MB04MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 171 "MB04MD.f"
	scale[i__] = 1.;
#line 172 "MB04MD.f"
/* L10: */
#line 172 "MB04MD.f"
    }

/*     Compute the 1-norm of matrix A and exit if it is zero. */

#line 176 "MB04MD.f"
    anorm = dlange_("1-norm", n, n, &a[a_offset], lda, &scale[1], (ftnlen)6);
#line 177 "MB04MD.f"
    if (anorm == 0.) {
#line 177 "MB04MD.f"
	return 0;
#line 177 "MB04MD.f"
    }

/*     Set some machine parameters and the maximum reduction in the */
/*     1-norm of A if zero rows or columns are encountered. */

#line 183 "MB04MD.f"
    sfmin1 = dlamch_("S", (ftnlen)1) / dlamch_("P", (ftnlen)1);
#line 184 "MB04MD.f"
    sfmax1 = 1. / sfmin1;
#line 185 "MB04MD.f"
    sfmin2 = sfmin1 * 10.;
#line 186 "MB04MD.f"
    sfmax2 = 1. / sfmin2;

#line 188 "MB04MD.f"
    sred = *maxred;
#line 189 "MB04MD.f"
    if (sred <= 0.) {
#line 189 "MB04MD.f"
	sred = 10.;
#line 189 "MB04MD.f"
    }

/* Computing MAX */
#line 191 "MB04MD.f"
    d__1 = anorm / sred;
#line 191 "MB04MD.f"
    maxnrm = max(d__1,sfmin1);

/*     Balance the matrix. */

/*     Iterative loop for norm reduction. */

#line 197 "MB04MD.f"
L20:
#line 198 "MB04MD.f"
    noconv = FALSE_;

#line 200 "MB04MD.f"
    i__1 = *n;
#line 200 "MB04MD.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 201 "MB04MD.f"
	c__ = 0.;
#line 202 "MB04MD.f"
	r__ = 0.;

#line 204 "MB04MD.f"
	i__2 = *n;
#line 204 "MB04MD.f"
	for (j = 1; j <= i__2; ++j) {
#line 205 "MB04MD.f"
	    if (j == i__) {
#line 205 "MB04MD.f"
		goto L30;
#line 205 "MB04MD.f"
	    }
#line 207 "MB04MD.f"
	    c__ += (d__1 = a[j + i__ * a_dim1], abs(d__1));
#line 208 "MB04MD.f"
	    r__ += (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 209 "MB04MD.f"
L30:
#line 209 "MB04MD.f"
	    ;
#line 209 "MB04MD.f"
	}
#line 210 "MB04MD.f"
	ica = idamax_(n, &a[i__ * a_dim1 + 1], &c__1);
#line 211 "MB04MD.f"
	ca = (d__1 = a[ica + i__ * a_dim1], abs(d__1));
#line 212 "MB04MD.f"
	ira = idamax_(n, &a[i__ + a_dim1], lda);
#line 213 "MB04MD.f"
	ra = (d__1 = a[i__ + ira * a_dim1], abs(d__1));

/*        Special case of zero C and/or R. */

#line 217 "MB04MD.f"
	if (c__ == 0. && r__ == 0.) {
#line 217 "MB04MD.f"
	    goto L80;
#line 217 "MB04MD.f"
	}
#line 219 "MB04MD.f"
	if (c__ == 0.) {
#line 220 "MB04MD.f"
	    if (r__ <= maxnrm) {
#line 220 "MB04MD.f"
		goto L80;
#line 220 "MB04MD.f"
	    }
#line 222 "MB04MD.f"
	    c__ = maxnrm;
#line 223 "MB04MD.f"
	}
#line 224 "MB04MD.f"
	if (r__ == 0.) {
#line 225 "MB04MD.f"
	    if (c__ <= maxnrm) {
#line 225 "MB04MD.f"
		goto L80;
#line 225 "MB04MD.f"
	    }
#line 227 "MB04MD.f"
	    r__ = maxnrm;
#line 228 "MB04MD.f"
	}

/*        Guard against zero C or R due to underflow. */

#line 232 "MB04MD.f"
	g = r__ / 10.;
#line 233 "MB04MD.f"
	f = 1.;
#line 234 "MB04MD.f"
	s = c__ + r__;
#line 235 "MB04MD.f"
L40:
/* Computing MAX */
#line 236 "MB04MD.f"
	d__1 = max(f,c__);
/* Computing MIN */
#line 236 "MB04MD.f"
	d__2 = min(r__,g);
#line 236 "MB04MD.f"
	if (c__ >= g || max(d__1,ca) >= sfmax2 || min(d__2,ra) <= sfmin2) {
#line 236 "MB04MD.f"
	    goto L50;
#line 236 "MB04MD.f"
	}
#line 238 "MB04MD.f"
	f *= 10.;
#line 239 "MB04MD.f"
	c__ *= 10.;
#line 240 "MB04MD.f"
	ca *= 10.;
#line 241 "MB04MD.f"
	r__ /= 10.;
#line 242 "MB04MD.f"
	g /= 10.;
#line 243 "MB04MD.f"
	ra /= 10.;
#line 244 "MB04MD.f"
	goto L40;

#line 246 "MB04MD.f"
L50:
#line 247 "MB04MD.f"
	g = c__ / 10.;
#line 248 "MB04MD.f"
L60:
/* Computing MIN */
#line 249 "MB04MD.f"
	d__1 = min(f,c__), d__1 = min(d__1,g);
#line 249 "MB04MD.f"
	if (g < r__ || max(r__,ra) >= sfmax2 || min(d__1,ca) <= sfmin2) {
#line 249 "MB04MD.f"
	    goto L70;
#line 249 "MB04MD.f"
	}
#line 251 "MB04MD.f"
	f /= 10.;
#line 252 "MB04MD.f"
	c__ /= 10.;
#line 253 "MB04MD.f"
	g /= 10.;
#line 254 "MB04MD.f"
	ca /= 10.;
#line 255 "MB04MD.f"
	r__ *= 10.;
#line 256 "MB04MD.f"
	ra *= 10.;
#line 257 "MB04MD.f"
	goto L60;

/*        Now balance. */

#line 261 "MB04MD.f"
L70:
#line 262 "MB04MD.f"
	if (c__ + r__ >= s * .95) {
#line 262 "MB04MD.f"
	    goto L80;
#line 262 "MB04MD.f"
	}
#line 264 "MB04MD.f"
	if (f < 1. && scale[i__] < 1.) {
#line 265 "MB04MD.f"
	    if (f * scale[i__] <= sfmin1) {
#line 265 "MB04MD.f"
		goto L80;
#line 265 "MB04MD.f"
	    }
#line 267 "MB04MD.f"
	}
#line 268 "MB04MD.f"
	if (f > 1. && scale[i__] > 1.) {
#line 269 "MB04MD.f"
	    if (scale[i__] >= sfmax1 / f) {
#line 269 "MB04MD.f"
		goto L80;
#line 269 "MB04MD.f"
	    }
#line 271 "MB04MD.f"
	}
#line 272 "MB04MD.f"
	g = 1. / f;
#line 273 "MB04MD.f"
	scale[i__] *= f;
#line 274 "MB04MD.f"
	noconv = TRUE_;

#line 276 "MB04MD.f"
	dscal_(n, &g, &a[i__ + a_dim1], lda);
#line 277 "MB04MD.f"
	dscal_(n, &f, &a[i__ * a_dim1 + 1], &c__1);

#line 279 "MB04MD.f"
L80:
#line 279 "MB04MD.f"
	;
#line 279 "MB04MD.f"
    }

#line 281 "MB04MD.f"
    if (noconv) {
#line 281 "MB04MD.f"
	goto L20;
#line 281 "MB04MD.f"
    }

/*     Set the norm reduction parameter. */

#line 286 "MB04MD.f"
    *maxred = anorm / dlange_("1-norm", n, n, &a[a_offset], lda, &scale[1], (
	    ftnlen)6);

#line 288 "MB04MD.f"
    return 0;
/* *** End of MB04MD *** */
} /* mb04md_ */

