#line 1 "SB04NX.f"
/* SB04NX.f -- translated by f2c (version 20100827).
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

#line 1 "SB04NX.f"
/* Table of constant values */

static integer c__2 = 2;
static doublereal c_b6 = 0.;
static integer c__1 = 1;

/* Subroutine */ int sb04nx_(char *rc, char *ul, integer *m, doublereal *a, 
	integer *lda, doublereal *lambd1, doublereal *lambd2, doublereal *
	lambd3, doublereal *lambd4, doublereal *d__, doublereal *tol, integer 
	*iwork, doublereal *dwork, integer *lddwor, integer *info, ftnlen 
	rc_len, ftnlen ul_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, dwork_dim1, dwork_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal c__;
    static integer j;
    static doublereal r__, s;
    static integer j1, j2, m2, mj, ml;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static char trans[1];
    extern /* Subroutine */ int dtrsv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dlaset_(char *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), dlartg_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);


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

/*     To solve a system of equations in Hessenberg form with two */
/*     consecutive offdiagonals and two right-hand sides. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     RC      CHARACTER*1 */
/*             Indicates processing by columns or rows, as follows: */
/*             = 'R':  Row transformations are applied; */
/*             = 'C':  Column transformations are applied. */

/*     UL      CHARACTER*1 */
/*             Indicates whether AB is upper or lower Hessenberg matrix, */
/*             as follows: */
/*             = 'U':  AB is upper Hessenberg; */
/*             = 'L':  AB is lower Hessenberg. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrix A.  M >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             The leading M-by-M part of this array must contain a */
/*             matrix A in Hessenberg form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     LAMBD1, (input) DOUBLE PRECISION */
/*     LAMBD2, These variables must contain the 2-by-2 block to be added */
/*     LAMBD3, to the diagonal blocks of A. */
/*     LAMBD4 */

/*     D       (input/output) DOUBLE PRECISION array, dimension (2*M) */
/*             On entry, this array must contain the two right-hand */
/*             side vectors of the Hessenberg system, stored row-wise. */
/*             On exit, if INFO = 0, this array contains the two solution */
/*             vectors of the Hessenberg system, stored row-wise. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used to test for near singularity of */
/*             the triangular factor R of the Hessenberg matrix. A matrix */
/*             whose estimated condition number is less than 1/TOL is */
/*             considered to be nonsingular. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDDWOR,2*M+3) */
/*             The leading 2*M-by-2*M part of this array is used for */
/*             computing the triangular factor of the QR decomposition */
/*             of the Hessenberg matrix. The remaining 6*M elements are */
/*             used as workspace for the computation of the reciprocal */
/*             condition estimate. */

/*     LDDWOR  INTEGER */
/*             The leading dimension of array DWORK. */
/*             LDDWOR >= MAX(1,2*M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if the Hessenberg matrix is (numerically) singular. */
/*                   That is, its estimated reciprocal condition number */
/*                   is less than or equal to TOL. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Aug. 1997. */
/*     Supersedes Release 2.0 routine SB04BX by M. Vanbegin, and */
/*     P. Van Dooren, Philips Research Laboratory, Brussels, Belgium. */

/*     REVISIONS */

/*     - */

/*     Note that RC, UL, M and LDA must be such that the value of the */
/*     LOGICAL variable OK in the following statement is true. */

/*      OK = ( ( UL.EQ.'U' ) .OR. ( UL.EQ.'u' ) .OR. */
/*             ( UL.EQ.'L' ) .OR. ( UL.EQ.'l' ) ) */
/*           .AND. */
/*           ( ( RC.EQ.'R' ) .OR. ( RC.EQ.'r' ) .OR. */
/*             ( RC.EQ.'C' ) .OR. ( RC.EQ.'c' ) ) */
/*           .AND. */
/*           ( M.GE.0 ) */
/*           .AND. */
/*           ( LDA.GE.MAX( 1, M ) ) */
/*           .AND. */
/*           ( LDDWOR.GE.MAX( 1, 2*M ) ) */

/*     KEYWORDS */

/*     Hessenberg form, orthogonal transformation, real Schur form, */
/*     Sylvester equation. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 154 "SB04NX.f"
    /* Parameter adjustments */
#line 154 "SB04NX.f"
    a_dim1 = *lda;
#line 154 "SB04NX.f"
    a_offset = 1 + a_dim1;
#line 154 "SB04NX.f"
    a -= a_offset;
#line 154 "SB04NX.f"
    --d__;
#line 154 "SB04NX.f"
    --iwork;
#line 154 "SB04NX.f"
    dwork_dim1 = *lddwor;
#line 154 "SB04NX.f"
    dwork_offset = 1 + dwork_dim1;
#line 154 "SB04NX.f"
    dwork -= dwork_offset;
#line 154 "SB04NX.f"

#line 154 "SB04NX.f"
    /* Function Body */
#line 154 "SB04NX.f"
    *info = 0;

/*     For speed, no tests on the input scalar arguments are made. */
/*     Quick return if possible. */

#line 159 "SB04NX.f"
    if (*m == 0) {
#line 159 "SB04NX.f"
	return 0;
#line 159 "SB04NX.f"
    }

#line 162 "SB04NX.f"
    m2 = *m << 1;
#line 163 "SB04NX.f"
    if (lsame_(ul, "U", (ftnlen)1, (ftnlen)1)) {

#line 165 "SB04NX.f"
	i__1 = *m;
#line 165 "SB04NX.f"
	for (j = 1; j <= i__1; ++j) {
#line 166 "SB04NX.f"
	    j2 = j << 1;
/* Computing MIN */
#line 167 "SB04NX.f"
	    i__2 = *m, i__3 = j + 1;
#line 167 "SB04NX.f"
	    ml = min(i__2,i__3);
#line 168 "SB04NX.f"
	    dlaset_("Full", &m2, &c__2, &c_b6, &c_b6, &dwork[(j2 - 1) * 
		    dwork_dim1 + 1], lddwor, (ftnlen)4);
#line 170 "SB04NX.f"
	    dcopy_(&ml, &a[j * a_dim1 + 1], &c__1, &dwork[(j2 - 1) * 
		    dwork_dim1 + 1], &c__2);
#line 171 "SB04NX.f"
	    dcopy_(&ml, &a[j * a_dim1 + 1], &c__1, &dwork[j2 * dwork_dim1 + 2]
		    , &c__2);
#line 172 "SB04NX.f"
	    dwork[j2 - 1 + (j2 - 1) * dwork_dim1] += *lambd1;
#line 173 "SB04NX.f"
	    dwork[j2 + (j2 - 1) * dwork_dim1] = *lambd3;
#line 174 "SB04NX.f"
	    dwork[j2 - 1 + j2 * dwork_dim1] = *lambd2;
#line 175 "SB04NX.f"
	    dwork[j2 + j2 * dwork_dim1] += *lambd4;
#line 176 "SB04NX.f"
/* L20: */
#line 176 "SB04NX.f"
	}

#line 178 "SB04NX.f"
	if (lsame_(rc, "R", (ftnlen)1, (ftnlen)1)) {
#line 179 "SB04NX.f"
	    *(unsigned char *)trans = 'N';

/*           A is an upper Hessenberg matrix, row transformations. */

#line 183 "SB04NX.f"
	    i__1 = m2 - 1;
#line 183 "SB04NX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 184 "SB04NX.f"
		mj = m2 - j;
#line 185 "SB04NX.f"
		if (j < m2 - 1) {
#line 186 "SB04NX.f"
		    if (dwork[j + 2 + j * dwork_dim1] != 0.) {
#line 187 "SB04NX.f"
			dlartg_(&dwork[j + 1 + j * dwork_dim1], &dwork[j + 2 
				+ j * dwork_dim1], &c__, &s, &r__);
#line 188 "SB04NX.f"
			dwork[j + 1 + j * dwork_dim1] = r__;
#line 189 "SB04NX.f"
			dwork[j + 2 + j * dwork_dim1] = 0.;
#line 190 "SB04NX.f"
			drot_(&mj, &dwork[j + 1 + (j + 1) * dwork_dim1], 
				lddwor, &dwork[j + 2 + (j + 1) * dwork_dim1], 
				lddwor, &c__, &s);
#line 192 "SB04NX.f"
			drot_(&c__1, &d__[j + 1], &c__1, &d__[j + 2], &c__1, &
				c__, &s);
#line 193 "SB04NX.f"
		    }
#line 194 "SB04NX.f"
		}
#line 195 "SB04NX.f"
		if (dwork[j + 1 + j * dwork_dim1] != 0.) {
#line 196 "SB04NX.f"
		    dlartg_(&dwork[j + j * dwork_dim1], &dwork[j + 1 + j * 
			    dwork_dim1], &c__, &s, &r__);
#line 197 "SB04NX.f"
		    dwork[j + j * dwork_dim1] = r__;
#line 198 "SB04NX.f"
		    dwork[j + 1 + j * dwork_dim1] = 0.;
#line 199 "SB04NX.f"
		    drot_(&mj, &dwork[j + (j + 1) * dwork_dim1], lddwor, &
			    dwork[j + 1 + (j + 1) * dwork_dim1], lddwor, &c__,
			     &s);
#line 201 "SB04NX.f"
		    drot_(&c__1, &d__[j], &c__1, &d__[j + 1], &c__1, &c__, &s)
			    ;
#line 202 "SB04NX.f"
		}
#line 203 "SB04NX.f"
/* L40: */
#line 203 "SB04NX.f"
	    }

#line 205 "SB04NX.f"
	} else {
#line 206 "SB04NX.f"
	    *(unsigned char *)trans = 'T';

/*           A is an upper Hessenberg matrix, column transformations. */

#line 210 "SB04NX.f"
	    i__1 = m2 - 1;
#line 210 "SB04NX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 211 "SB04NX.f"
		mj = m2 - j;
#line 212 "SB04NX.f"
		if (j < m2 - 1) {
#line 213 "SB04NX.f"
		    if (dwork[mj + 1 + (mj - 1) * dwork_dim1] != 0.) {
#line 214 "SB04NX.f"
			dlartg_(&dwork[mj + 1 + mj * dwork_dim1], &dwork[mj + 
				1 + (mj - 1) * dwork_dim1], &c__, &s, &r__);
#line 216 "SB04NX.f"
			dwork[mj + 1 + mj * dwork_dim1] = r__;
#line 217 "SB04NX.f"
			dwork[mj + 1 + (mj - 1) * dwork_dim1] = 0.;
#line 218 "SB04NX.f"
			drot_(&mj, &dwork[mj * dwork_dim1 + 1], &c__1, &dwork[
				(mj - 1) * dwork_dim1 + 1], &c__1, &c__, &s);
#line 220 "SB04NX.f"
			drot_(&c__1, &d__[mj], &c__1, &d__[mj - 1], &c__1, &
				c__, &s);
#line 221 "SB04NX.f"
		    }
#line 222 "SB04NX.f"
		}
#line 223 "SB04NX.f"
		if (dwork[mj + 1 + mj * dwork_dim1] != 0.) {
#line 224 "SB04NX.f"
		    dlartg_(&dwork[mj + 1 + (mj + 1) * dwork_dim1], &dwork[mj 
			    + 1 + mj * dwork_dim1], &c__, &s, &r__);
#line 226 "SB04NX.f"
		    dwork[mj + 1 + (mj + 1) * dwork_dim1] = r__;
#line 227 "SB04NX.f"
		    dwork[mj + 1 + mj * dwork_dim1] = 0.;
#line 228 "SB04NX.f"
		    drot_(&mj, &dwork[(mj + 1) * dwork_dim1 + 1], &c__1, &
			    dwork[mj * dwork_dim1 + 1], &c__1, &c__, &s);
#line 230 "SB04NX.f"
		    drot_(&c__1, &d__[mj + 1], &c__1, &d__[mj], &c__1, &c__, &
			    s);
#line 231 "SB04NX.f"
		}
#line 232 "SB04NX.f"
/* L60: */
#line 232 "SB04NX.f"
	    }

#line 234 "SB04NX.f"
	}
#line 235 "SB04NX.f"
    } else {

#line 237 "SB04NX.f"
	i__1 = *m;
#line 237 "SB04NX.f"
	for (j = 1; j <= i__1; ++j) {
#line 238 "SB04NX.f"
	    j2 = j << 1;
/* Computing MAX */
#line 239 "SB04NX.f"
	    i__2 = j - 1;
#line 239 "SB04NX.f"
	    j1 = max(i__2,1);
/* Computing MIN */
#line 240 "SB04NX.f"
	    i__2 = *m - j + 2;
#line 240 "SB04NX.f"
	    ml = min(i__2,*m);
#line 241 "SB04NX.f"
	    dlaset_("Full", &m2, &c__2, &c_b6, &c_b6, &dwork[(j2 - 1) * 
		    dwork_dim1 + 1], lddwor, (ftnlen)4);
#line 243 "SB04NX.f"
	    dcopy_(&ml, &a[j1 + j * a_dim1], &c__1, &dwork[(j1 << 1) - 1 + (
		    j2 - 1) * dwork_dim1], &c__2);
#line 244 "SB04NX.f"
	    dcopy_(&ml, &a[j1 + j * a_dim1], &c__1, &dwork[(j1 << 1) + j2 * 
		    dwork_dim1], &c__2);
#line 245 "SB04NX.f"
	    dwork[j2 - 1 + (j2 - 1) * dwork_dim1] += *lambd1;
#line 246 "SB04NX.f"
	    dwork[j2 + (j2 - 1) * dwork_dim1] = *lambd3;
#line 247 "SB04NX.f"
	    dwork[j2 - 1 + j2 * dwork_dim1] = *lambd2;
#line 248 "SB04NX.f"
	    dwork[j2 + j2 * dwork_dim1] += *lambd4;
#line 249 "SB04NX.f"
/* L80: */
#line 249 "SB04NX.f"
	}

#line 251 "SB04NX.f"
	if (lsame_(rc, "R", (ftnlen)1, (ftnlen)1)) {
#line 252 "SB04NX.f"
	    *(unsigned char *)trans = 'N';

/*           A is a lower Hessenberg matrix, row transformations. */

#line 256 "SB04NX.f"
	    i__1 = m2 - 1;
#line 256 "SB04NX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 257 "SB04NX.f"
		mj = m2 - j;
#line 258 "SB04NX.f"
		if (j < m2 - 1) {
#line 259 "SB04NX.f"
		    if (dwork[mj - 1 + (mj + 1) * dwork_dim1] != 0.) {
#line 260 "SB04NX.f"
			dlartg_(&dwork[mj + (mj + 1) * dwork_dim1], &dwork[mj 
				- 1 + (mj + 1) * dwork_dim1], &c__, &s, &r__);
#line 262 "SB04NX.f"
			dwork[mj + (mj + 1) * dwork_dim1] = r__;
#line 263 "SB04NX.f"
			dwork[mj - 1 + (mj + 1) * dwork_dim1] = 0.;
#line 264 "SB04NX.f"
			drot_(&mj, &dwork[mj + dwork_dim1], lddwor, &dwork[mj 
				- 1 + dwork_dim1], lddwor, &c__, &s);
#line 266 "SB04NX.f"
			drot_(&c__1, &d__[mj], &c__1, &d__[mj - 1], &c__1, &
				c__, &s);
#line 267 "SB04NX.f"
		    }
#line 268 "SB04NX.f"
		}
#line 269 "SB04NX.f"
		if (dwork[mj + (mj + 1) * dwork_dim1] != 0.) {
#line 270 "SB04NX.f"
		    dlartg_(&dwork[mj + 1 + (mj + 1) * dwork_dim1], &dwork[mj 
			    + (mj + 1) * dwork_dim1], &c__, &s, &r__);
#line 272 "SB04NX.f"
		    dwork[mj + 1 + (mj + 1) * dwork_dim1] = r__;
#line 273 "SB04NX.f"
		    dwork[mj + (mj + 1) * dwork_dim1] = 0.;
#line 274 "SB04NX.f"
		    drot_(&mj, &dwork[mj + 1 + dwork_dim1], lddwor, &dwork[mj 
			    + dwork_dim1], lddwor, &c__, &s);
#line 276 "SB04NX.f"
		    drot_(&c__1, &d__[mj + 1], &c__1, &d__[mj], &c__1, &c__, &
			    s);
#line 277 "SB04NX.f"
		}
#line 278 "SB04NX.f"
/* L100: */
#line 278 "SB04NX.f"
	    }

#line 280 "SB04NX.f"
	} else {
#line 281 "SB04NX.f"
	    *(unsigned char *)trans = 'T';

/*           A is a lower Hessenberg matrix, column transformations. */

#line 285 "SB04NX.f"
	    i__1 = m2 - 1;
#line 285 "SB04NX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 286 "SB04NX.f"
		mj = m2 - j;
#line 287 "SB04NX.f"
		if (j < m2 - 1) {
#line 288 "SB04NX.f"
		    if (dwork[j + (j + 2) * dwork_dim1] != 0.) {
#line 289 "SB04NX.f"
			dlartg_(&dwork[j + (j + 1) * dwork_dim1], &dwork[j + (
				j + 2) * dwork_dim1], &c__, &s, &r__);
#line 290 "SB04NX.f"
			dwork[j + (j + 1) * dwork_dim1] = r__;
#line 291 "SB04NX.f"
			dwork[j + (j + 2) * dwork_dim1] = 0.;
#line 292 "SB04NX.f"
			drot_(&mj, &dwork[j + 1 + (j + 1) * dwork_dim1], &
				c__1, &dwork[j + 1 + (j + 2) * dwork_dim1], &
				c__1, &c__, &s);
#line 294 "SB04NX.f"
			drot_(&c__1, &d__[j + 1], &c__1, &d__[j + 2], &c__1, &
				c__, &s);
#line 295 "SB04NX.f"
		    }
#line 296 "SB04NX.f"
		}
#line 297 "SB04NX.f"
		if (dwork[j + (j + 1) * dwork_dim1] != 0.) {
#line 298 "SB04NX.f"
		    dlartg_(&dwork[j + j * dwork_dim1], &dwork[j + (j + 1) * 
			    dwork_dim1], &c__, &s, &r__);
#line 299 "SB04NX.f"
		    dwork[j + j * dwork_dim1] = r__;
#line 300 "SB04NX.f"
		    dwork[j + (j + 1) * dwork_dim1] = 0.;
#line 301 "SB04NX.f"
		    drot_(&mj, &dwork[j + 1 + j * dwork_dim1], &c__1, &dwork[
			    j + 1 + (j + 1) * dwork_dim1], &c__1, &c__, &s);
#line 303 "SB04NX.f"
		    drot_(&c__1, &d__[j], &c__1, &d__[j + 1], &c__1, &c__, &s)
			    ;
#line 304 "SB04NX.f"
		}
#line 305 "SB04NX.f"
/* L120: */
#line 305 "SB04NX.f"
	    }

#line 307 "SB04NX.f"
	}
#line 308 "SB04NX.f"
    }

#line 310 "SB04NX.f"
    dtrcon_("1-norm", ul, "Non-unit", &m2, &dwork[dwork_offset], lddwor, &
	    rcond, &dwork[(m2 + 1) * dwork_dim1 + 1], &iwork[1], info, (
	    ftnlen)6, (ftnlen)1, (ftnlen)8);
#line 312 "SB04NX.f"
    if (rcond <= *tol) {
#line 313 "SB04NX.f"
	*info = 1;
#line 314 "SB04NX.f"
    } else {
#line 315 "SB04NX.f"
	dtrsv_(ul, trans, "Non-unit", &m2, &dwork[dwork_offset], lddwor, &d__[
		1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 316 "SB04NX.f"
    }

#line 318 "SB04NX.f"
    return 0;
/* *** Last line of SB04NX *** */
} /* sb04nx_ */

