#line 1 "SB04RX.f"
/* SB04RX.f -- translated by f2c (version 20100827).
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

#line 1 "SB04RX.f"
/* Table of constant values */

static integer c__2 = 2;
static doublereal c_b6 = 0.;
static integer c__1 = 1;

/* Subroutine */ int sb04rx_(char *rc, char *ul, integer *m, doublereal *a, 
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
	    doublereal *, integer *, doublereal *, doublereal *), dscal_(
	    integer *, doublereal *, doublereal *, integer *);
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

/*     To solve a system of equations in quasi-Hessenberg form */
/*     (Hessenberg form plus two consecutive offdiagonals) with two */
/*     right-hand sides. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     RC      CHARACTER*1 */
/*             Indicates processing by columns or rows, as follows: */
/*             = 'R':  Row transformations are applied; */
/*             = 'C':  Column transformations are applied. */

/*     UL      CHARACTER*1 */
/*             Indicates whether A is upper or lower Hessenberg matrix, */
/*             as follows: */
/*             = 'U':  A is upper Hessenberg; */
/*             = 'L':  A is lower Hessenberg. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The order of the matrix A.  M >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,M) */
/*             The leading M-by-M part of this array must contain a */
/*             matrix A in Hessenberg form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     LAMBD1, (input) DOUBLE PRECISION */
/*     LAMBD2, These variables must contain the 2-by-2 block to be */
/*     LAMBD3, multiplied to the elements of A. */
/*     LAMBD4 */

/*     D       (input/output) DOUBLE PRECISION array, dimension (2*M) */
/*             On entry, this array must contain the two right-hand */
/*             side vectors of the quasi-Hessenberg system, stored */
/*             row-wise. */
/*             On exit, if INFO = 0, this array contains the two solution */
/*             vectors of the quasi-Hessenberg system, stored row-wise. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used to test for near singularity of */
/*             the triangular factor R of the quasi-Hessenberg matrix. */
/*             A matrix whose estimated condition number is less */
/*             than 1/TOL is considered to be nonsingular. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (2*M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDDWOR,2*M+3) */
/*             The leading 2*M-by-2*M part of this array is used for */
/*             computing the triangular factor of the QR decomposition */
/*             of the quasi-Hessenberg matrix. The remaining 6*M elements */
/*             are used as workspace for the computation of the */
/*             reciprocal condition estimate. */

/*     LDDWOR  INTEGER */
/*             The leading dimension of array DWORK. */
/*             LDDWOR >= MAX(1,2*M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if the quasi-Hessenberg matrix is (numerically) */
/*                   singular. That is, its estimated reciprocal */
/*                   condition number is less than or equal to TOL. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTORS */

/*     D. Sima, University of Bucharest, May 2000. */

/*     REVISIONS */

/*     - */

/*     Note that RC, UL, M, LDA, and LDDWOR must be such that the value */
/*     of the LOGICAL variable OK in the following statement is true. */

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

/*     These conditions are not checked by the routine. */

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

#line 157 "SB04RX.f"
    /* Parameter adjustments */
#line 157 "SB04RX.f"
    a_dim1 = *lda;
#line 157 "SB04RX.f"
    a_offset = 1 + a_dim1;
#line 157 "SB04RX.f"
    a -= a_offset;
#line 157 "SB04RX.f"
    --d__;
#line 157 "SB04RX.f"
    --iwork;
#line 157 "SB04RX.f"
    dwork_dim1 = *lddwor;
#line 157 "SB04RX.f"
    dwork_offset = 1 + dwork_dim1;
#line 157 "SB04RX.f"
    dwork -= dwork_offset;
#line 157 "SB04RX.f"

#line 157 "SB04RX.f"
    /* Function Body */
#line 157 "SB04RX.f"
    *info = 0;

/*     For speed, no tests on the input scalar arguments are made. */
/*     Quick return if possible. */

#line 162 "SB04RX.f"
    if (*m == 0) {
#line 162 "SB04RX.f"
	return 0;
#line 162 "SB04RX.f"
    }

#line 165 "SB04RX.f"
    m2 = *m << 1;
#line 166 "SB04RX.f"
    if (lsame_(ul, "U", (ftnlen)1, (ftnlen)1)) {

#line 168 "SB04RX.f"
	i__1 = *m;
#line 168 "SB04RX.f"
	for (j = 1; j <= i__1; ++j) {
#line 169 "SB04RX.f"
	    j2 = j << 1;
/* Computing MIN */
#line 170 "SB04RX.f"
	    i__2 = *m, i__3 = j + 1;
#line 170 "SB04RX.f"
	    ml = min(i__2,i__3);
#line 171 "SB04RX.f"
	    dlaset_("Full", &m2, &c__2, &c_b6, &c_b6, &dwork[(j2 - 1) * 
		    dwork_dim1 + 1], lddwor, (ftnlen)4);
#line 173 "SB04RX.f"
	    dcopy_(&ml, &a[j * a_dim1 + 1], &c__1, &dwork[(j2 - 1) * 
		    dwork_dim1 + 1], &c__2);
#line 174 "SB04RX.f"
	    dscal_(&ml, lambd1, &dwork[(j2 - 1) * dwork_dim1 + 1], &c__2);
#line 175 "SB04RX.f"
	    dcopy_(&ml, &a[j * a_dim1 + 1], &c__1, &dwork[(j2 - 1) * 
		    dwork_dim1 + 2], &c__2);
#line 176 "SB04RX.f"
	    dscal_(&ml, lambd3, &dwork[(j2 - 1) * dwork_dim1 + 2], &c__2);
#line 177 "SB04RX.f"
	    dcopy_(&ml, &a[j * a_dim1 + 1], &c__1, &dwork[j2 * dwork_dim1 + 1]
		    , &c__2);
#line 178 "SB04RX.f"
	    dscal_(&ml, lambd2, &dwork[j2 * dwork_dim1 + 1], &c__2);
#line 179 "SB04RX.f"
	    dcopy_(&ml, &a[j * a_dim1 + 1], &c__1, &dwork[j2 * dwork_dim1 + 2]
		    , &c__2);
#line 180 "SB04RX.f"
	    dscal_(&ml, lambd4, &dwork[j2 * dwork_dim1 + 2], &c__2);

#line 182 "SB04RX.f"
	    dwork[j2 - 1 + (j2 - 1) * dwork_dim1] += 1.;
#line 183 "SB04RX.f"
	    dwork[j2 + j2 * dwork_dim1] += 1.;
#line 184 "SB04RX.f"
/* L20: */
#line 184 "SB04RX.f"
	}

#line 186 "SB04RX.f"
	if (lsame_(rc, "R", (ftnlen)1, (ftnlen)1)) {
#line 187 "SB04RX.f"
	    *(unsigned char *)trans = 'N';

/*           A is an upper Hessenberg matrix, row transformations. */

#line 191 "SB04RX.f"
	    i__1 = m2 - 1;
#line 191 "SB04RX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 192 "SB04RX.f"
		mj = m2 - j;
#line 193 "SB04RX.f"
		if (j % 2 == 1 && j < m2 - 2) {
#line 194 "SB04RX.f"
		    if (dwork[j + 3 + j * dwork_dim1] != 0.) {
#line 195 "SB04RX.f"
			dlartg_(&dwork[j + 2 + j * dwork_dim1], &dwork[j + 3 
				+ j * dwork_dim1], &c__, &s, &r__);
#line 196 "SB04RX.f"
			dwork[j + 2 + j * dwork_dim1] = r__;
#line 197 "SB04RX.f"
			dwork[j + 3 + j * dwork_dim1] = 0.;
#line 198 "SB04RX.f"
			drot_(&mj, &dwork[j + 2 + (j + 1) * dwork_dim1], 
				lddwor, &dwork[j + 3 + (j + 1) * dwork_dim1], 
				lddwor, &c__, &s);
#line 200 "SB04RX.f"
			drot_(&c__1, &d__[j + 2], &c__1, &d__[j + 3], &c__1, &
				c__, &s);
#line 201 "SB04RX.f"
		    }
#line 202 "SB04RX.f"
		}
#line 203 "SB04RX.f"
		if (j < m2 - 1) {
#line 204 "SB04RX.f"
		    if (dwork[j + 2 + j * dwork_dim1] != 0.) {
#line 205 "SB04RX.f"
			dlartg_(&dwork[j + 1 + j * dwork_dim1], &dwork[j + 2 
				+ j * dwork_dim1], &c__, &s, &r__);
#line 206 "SB04RX.f"
			dwork[j + 1 + j * dwork_dim1] = r__;
#line 207 "SB04RX.f"
			dwork[j + 2 + j * dwork_dim1] = 0.;
#line 208 "SB04RX.f"
			drot_(&mj, &dwork[j + 1 + (j + 1) * dwork_dim1], 
				lddwor, &dwork[j + 2 + (j + 1) * dwork_dim1], 
				lddwor, &c__, &s);
#line 210 "SB04RX.f"
			drot_(&c__1, &d__[j + 1], &c__1, &d__[j + 2], &c__1, &
				c__, &s);
#line 211 "SB04RX.f"
		    }
#line 212 "SB04RX.f"
		}
#line 213 "SB04RX.f"
		if (dwork[j + 1 + j * dwork_dim1] != 0.) {
#line 214 "SB04RX.f"
		    dlartg_(&dwork[j + j * dwork_dim1], &dwork[j + 1 + j * 
			    dwork_dim1], &c__, &s, &r__);
#line 215 "SB04RX.f"
		    dwork[j + j * dwork_dim1] = r__;
#line 216 "SB04RX.f"
		    dwork[j + 1 + j * dwork_dim1] = 0.;
#line 217 "SB04RX.f"
		    drot_(&mj, &dwork[j + (j + 1) * dwork_dim1], lddwor, &
			    dwork[j + 1 + (j + 1) * dwork_dim1], lddwor, &c__,
			     &s);
#line 219 "SB04RX.f"
		    drot_(&c__1, &d__[j], &c__1, &d__[j + 1], &c__1, &c__, &s)
			    ;
#line 220 "SB04RX.f"
		}
#line 221 "SB04RX.f"
/* L40: */
#line 221 "SB04RX.f"
	    }

#line 223 "SB04RX.f"
	} else {
#line 224 "SB04RX.f"
	    *(unsigned char *)trans = 'T';

/*           A is an upper Hessenberg matrix, column transformations. */

#line 228 "SB04RX.f"
	    i__1 = m2 - 1;
#line 228 "SB04RX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 229 "SB04RX.f"
		mj = m2 - j;
#line 230 "SB04RX.f"
		if (j % 2 == 1 && j < m2 - 2) {
#line 231 "SB04RX.f"
		    if (dwork[mj + 1 + (mj - 2) * dwork_dim1] != 0.) {
#line 232 "SB04RX.f"
			dlartg_(&dwork[mj + 1 + (mj - 1) * dwork_dim1], &
				dwork[mj + 1 + (mj - 2) * dwork_dim1], &c__, &
				s, &r__);
#line 234 "SB04RX.f"
			dwork[mj + 1 + (mj - 1) * dwork_dim1] = r__;
#line 235 "SB04RX.f"
			dwork[mj + 1 + (mj - 2) * dwork_dim1] = 0.;
#line 236 "SB04RX.f"
			drot_(&mj, &dwork[(mj - 1) * dwork_dim1 + 1], &c__1, &
				dwork[(mj - 2) * dwork_dim1 + 1], &c__1, &c__,
				 &s);
#line 238 "SB04RX.f"
			drot_(&c__1, &d__[mj - 1], &c__1, &d__[mj - 2], &c__1,
				 &c__, &s);
#line 239 "SB04RX.f"
		    }
#line 240 "SB04RX.f"
		}
#line 241 "SB04RX.f"
		if (j < m2 - 1) {
#line 242 "SB04RX.f"
		    if (dwork[mj + 1 + (mj - 1) * dwork_dim1] != 0.) {
#line 243 "SB04RX.f"
			dlartg_(&dwork[mj + 1 + mj * dwork_dim1], &dwork[mj + 
				1 + (mj - 1) * dwork_dim1], &c__, &s, &r__);
#line 245 "SB04RX.f"
			dwork[mj + 1 + mj * dwork_dim1] = r__;
#line 246 "SB04RX.f"
			dwork[mj + 1 + (mj - 1) * dwork_dim1] = 0.;
#line 247 "SB04RX.f"
			drot_(&mj, &dwork[mj * dwork_dim1 + 1], &c__1, &dwork[
				(mj - 1) * dwork_dim1 + 1], &c__1, &c__, &s);
#line 249 "SB04RX.f"
			drot_(&c__1, &d__[mj], &c__1, &d__[mj - 1], &c__1, &
				c__, &s);
#line 250 "SB04RX.f"
		    }
#line 251 "SB04RX.f"
		}
#line 252 "SB04RX.f"
		if (dwork[mj + 1 + mj * dwork_dim1] != 0.) {
#line 253 "SB04RX.f"
		    dlartg_(&dwork[mj + 1 + (mj + 1) * dwork_dim1], &dwork[mj 
			    + 1 + mj * dwork_dim1], &c__, &s, &r__);
#line 255 "SB04RX.f"
		    dwork[mj + 1 + (mj + 1) * dwork_dim1] = r__;
#line 256 "SB04RX.f"
		    dwork[mj + 1 + mj * dwork_dim1] = 0.;
#line 257 "SB04RX.f"
		    drot_(&mj, &dwork[(mj + 1) * dwork_dim1 + 1], &c__1, &
			    dwork[mj * dwork_dim1 + 1], &c__1, &c__, &s);
#line 259 "SB04RX.f"
		    drot_(&c__1, &d__[mj + 1], &c__1, &d__[mj], &c__1, &c__, &
			    s);
#line 260 "SB04RX.f"
		}
#line 261 "SB04RX.f"
/* L60: */
#line 261 "SB04RX.f"
	    }

#line 263 "SB04RX.f"
	}
#line 264 "SB04RX.f"
    } else {

#line 266 "SB04RX.f"
	i__1 = *m;
#line 266 "SB04RX.f"
	for (j = 1; j <= i__1; ++j) {
#line 267 "SB04RX.f"
	    j2 = j << 1;
/* Computing MAX */
#line 268 "SB04RX.f"
	    i__2 = j - 1;
#line 268 "SB04RX.f"
	    j1 = max(i__2,1);
/* Computing MIN */
#line 269 "SB04RX.f"
	    i__2 = *m - j + 2;
#line 269 "SB04RX.f"
	    ml = min(i__2,*m);
#line 270 "SB04RX.f"
	    dlaset_("Full", &m2, &c__2, &c_b6, &c_b6, &dwork[(j2 - 1) * 
		    dwork_dim1 + 1], lddwor, (ftnlen)4);
#line 272 "SB04RX.f"
	    dcopy_(&ml, &a[j1 + j * a_dim1], &c__1, &dwork[(j1 << 1) - 1 + (
		    j2 - 1) * dwork_dim1], &c__2);
#line 273 "SB04RX.f"
	    dscal_(&ml, lambd1, &dwork[(j1 << 1) - 1 + (j2 - 1) * dwork_dim1],
		     &c__2);
#line 274 "SB04RX.f"
	    dcopy_(&ml, &a[j1 + j * a_dim1], &c__1, &dwork[(j1 << 1) + (j2 - 
		    1) * dwork_dim1], &c__2);
#line 275 "SB04RX.f"
	    dscal_(&ml, lambd3, &dwork[(j1 << 1) + (j2 - 1) * dwork_dim1], &
		    c__2);
#line 276 "SB04RX.f"
	    dcopy_(&ml, &a[j1 + j * a_dim1], &c__1, &dwork[(j1 << 1) - 1 + j2 
		    * dwork_dim1], &c__2);
#line 277 "SB04RX.f"
	    dscal_(&ml, lambd2, &dwork[(j1 << 1) - 1 + j2 * dwork_dim1], &
		    c__2);
#line 278 "SB04RX.f"
	    dcopy_(&ml, &a[j1 + j * a_dim1], &c__1, &dwork[(j1 << 1) + j2 * 
		    dwork_dim1], &c__2);
#line 279 "SB04RX.f"
	    dscal_(&ml, lambd4, &dwork[(j1 << 1) + j2 * dwork_dim1], &c__2);

#line 281 "SB04RX.f"
	    dwork[j2 - 1 + (j2 - 1) * dwork_dim1] += 1.;
#line 282 "SB04RX.f"
	    dwork[j2 + j2 * dwork_dim1] += 1.;
#line 283 "SB04RX.f"
/* L80: */
#line 283 "SB04RX.f"
	}

#line 285 "SB04RX.f"
	if (lsame_(rc, "R", (ftnlen)1, (ftnlen)1)) {
#line 286 "SB04RX.f"
	    *(unsigned char *)trans = 'N';

/*           A is a lower Hessenberg matrix, row transformations. */

#line 290 "SB04RX.f"
	    i__1 = m2 - 1;
#line 290 "SB04RX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 291 "SB04RX.f"
		mj = m2 - j;
#line 292 "SB04RX.f"
		if (j % 2 == 1 && j < m2 - 2) {
#line 293 "SB04RX.f"
		    if (dwork[mj - 2 + (mj + 1) * dwork_dim1] != 0.) {
#line 294 "SB04RX.f"
			dlartg_(&dwork[mj - 1 + (mj + 1) * dwork_dim1], &
				dwork[mj - 2 + (mj + 1) * dwork_dim1], &c__, &
				s, &r__);
#line 296 "SB04RX.f"
			dwork[mj - 1 + (mj + 1) * dwork_dim1] = r__;
#line 297 "SB04RX.f"
			dwork[mj - 2 + (mj + 1) * dwork_dim1] = 0.;
#line 298 "SB04RX.f"
			drot_(&mj, &dwork[mj - 1 + dwork_dim1], lddwor, &
				dwork[mj - 2 + dwork_dim1], lddwor, &c__, &s);
#line 300 "SB04RX.f"
			drot_(&c__1, &d__[mj - 1], &c__1, &d__[mj - 2], &c__1,
				 &c__, &s);
#line 301 "SB04RX.f"
		    }
#line 302 "SB04RX.f"
		}
#line 303 "SB04RX.f"
		if (j < m2 - 1) {
#line 304 "SB04RX.f"
		    if (dwork[mj - 1 + (mj + 1) * dwork_dim1] != 0.) {
#line 305 "SB04RX.f"
			dlartg_(&dwork[mj + (mj + 1) * dwork_dim1], &dwork[mj 
				- 1 + (mj + 1) * dwork_dim1], &c__, &s, &r__);
#line 307 "SB04RX.f"
			dwork[mj + (mj + 1) * dwork_dim1] = r__;
#line 308 "SB04RX.f"
			dwork[mj - 1 + (mj + 1) * dwork_dim1] = 0.;
#line 309 "SB04RX.f"
			drot_(&mj, &dwork[mj + dwork_dim1], lddwor, &dwork[mj 
				- 1 + dwork_dim1], lddwor, &c__, &s);
#line 311 "SB04RX.f"
			drot_(&c__1, &d__[mj], &c__1, &d__[mj - 1], &c__1, &
				c__, &s);
#line 312 "SB04RX.f"
		    }
#line 313 "SB04RX.f"
		}
#line 314 "SB04RX.f"
		if (dwork[mj + (mj + 1) * dwork_dim1] != 0.) {
#line 315 "SB04RX.f"
		    dlartg_(&dwork[mj + 1 + (mj + 1) * dwork_dim1], &dwork[mj 
			    + (mj + 1) * dwork_dim1], &c__, &s, &r__);
#line 317 "SB04RX.f"
		    dwork[mj + 1 + (mj + 1) * dwork_dim1] = r__;
#line 318 "SB04RX.f"
		    dwork[mj + (mj + 1) * dwork_dim1] = 0.;
#line 319 "SB04RX.f"
		    drot_(&mj, &dwork[mj + 1 + dwork_dim1], lddwor, &dwork[mj 
			    + dwork_dim1], lddwor, &c__, &s);
#line 321 "SB04RX.f"
		    drot_(&c__1, &d__[mj + 1], &c__1, &d__[mj], &c__1, &c__, &
			    s);
#line 322 "SB04RX.f"
		}
#line 323 "SB04RX.f"
/* L100: */
#line 323 "SB04RX.f"
	    }

#line 325 "SB04RX.f"
	} else {
#line 326 "SB04RX.f"
	    *(unsigned char *)trans = 'T';

/*           A is a lower Hessenberg matrix, column transformations. */

#line 330 "SB04RX.f"
	    i__1 = m2 - 1;
#line 330 "SB04RX.f"
	    for (j = 1; j <= i__1; ++j) {
#line 331 "SB04RX.f"
		mj = m2 - j;
#line 332 "SB04RX.f"
		if (j % 2 == 1 && j < m2 - 2) {
#line 333 "SB04RX.f"
		    if (dwork[j + (j + 3) * dwork_dim1] != 0.) {
#line 334 "SB04RX.f"
			dlartg_(&dwork[j + (j + 2) * dwork_dim1], &dwork[j + (
				j + 3) * dwork_dim1], &c__, &s, &r__);
#line 335 "SB04RX.f"
			dwork[j + (j + 2) * dwork_dim1] = r__;
#line 336 "SB04RX.f"
			dwork[j + (j + 3) * dwork_dim1] = 0.;
#line 337 "SB04RX.f"
			drot_(&mj, &dwork[j + 1 + (j + 2) * dwork_dim1], &
				c__1, &dwork[j + 1 + (j + 3) * dwork_dim1], &
				c__1, &c__, &s);
#line 339 "SB04RX.f"
			drot_(&c__1, &d__[j + 2], &c__1, &d__[j + 3], &c__1, &
				c__, &s);
#line 340 "SB04RX.f"
		    }
#line 341 "SB04RX.f"
		}
#line 342 "SB04RX.f"
		if (j < m2 - 1) {
#line 343 "SB04RX.f"
		    if (dwork[j + (j + 2) * dwork_dim1] != 0.) {
#line 344 "SB04RX.f"
			dlartg_(&dwork[j + (j + 1) * dwork_dim1], &dwork[j + (
				j + 2) * dwork_dim1], &c__, &s, &r__);
#line 345 "SB04RX.f"
			dwork[j + (j + 1) * dwork_dim1] = r__;
#line 346 "SB04RX.f"
			dwork[j + (j + 2) * dwork_dim1] = 0.;
#line 347 "SB04RX.f"
			drot_(&mj, &dwork[j + 1 + (j + 1) * dwork_dim1], &
				c__1, &dwork[j + 1 + (j + 2) * dwork_dim1], &
				c__1, &c__, &s);
#line 349 "SB04RX.f"
			drot_(&c__1, &d__[j + 1], &c__1, &d__[j + 2], &c__1, &
				c__, &s);
#line 350 "SB04RX.f"
		    }
#line 351 "SB04RX.f"
		}
#line 352 "SB04RX.f"
		if (dwork[j + (j + 1) * dwork_dim1] != 0.) {
#line 353 "SB04RX.f"
		    dlartg_(&dwork[j + j * dwork_dim1], &dwork[j + (j + 1) * 
			    dwork_dim1], &c__, &s, &r__);
#line 354 "SB04RX.f"
		    dwork[j + j * dwork_dim1] = r__;
#line 355 "SB04RX.f"
		    dwork[j + (j + 1) * dwork_dim1] = 0.;
#line 356 "SB04RX.f"
		    drot_(&mj, &dwork[j + 1 + j * dwork_dim1], &c__1, &dwork[
			    j + 1 + (j + 1) * dwork_dim1], &c__1, &c__, &s);
#line 358 "SB04RX.f"
		    drot_(&c__1, &d__[j], &c__1, &d__[j + 1], &c__1, &c__, &s)
			    ;
#line 359 "SB04RX.f"
		}
#line 360 "SB04RX.f"
/* L120: */
#line 360 "SB04RX.f"
	    }

#line 362 "SB04RX.f"
	}
#line 363 "SB04RX.f"
    }

#line 365 "SB04RX.f"
    dtrcon_("1-norm", ul, "Non-unit", &m2, &dwork[dwork_offset], lddwor, &
	    rcond, &dwork[(m2 + 1) * dwork_dim1 + 1], &iwork[1], info, (
	    ftnlen)6, (ftnlen)1, (ftnlen)8);
#line 367 "SB04RX.f"
    if (rcond <= *tol) {
#line 368 "SB04RX.f"
	*info = 1;
#line 369 "SB04RX.f"
    } else {
#line 370 "SB04RX.f"
	dtrsv_(ul, trans, "Non-unit", &m2, &dwork[dwork_offset], lddwor, &d__[
		1], &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 371 "SB04RX.f"
    }

#line 373 "SB04RX.f"
    return 0;
/* *** Last line of SB04RX *** */
} /* sb04rx_ */

