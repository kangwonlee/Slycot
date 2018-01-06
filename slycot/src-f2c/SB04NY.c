#line 1 "SB04NY.f"
/* SB04NY.f -- translated by f2c (version 20100827).
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

#line 1 "SB04NY.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int sb04ny_(char *rc, char *ul, integer *m, doublereal *a, 
	integer *lda, doublereal *lambda, doublereal *d__, doublereal *tol, 
	integer *iwork, doublereal *dwork, integer *lddwor, integer *info, 
	ftnlen rc_len, ftnlen ul_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, dwork_dim1, dwork_offset, i__1, i__2, i__3;

    /* Local variables */
    static doublereal c__;
    static integer j;
    static doublereal r__, s;
    static integer j1, mj;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal rcond;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static char trans[1];
    extern /* Subroutine */ int dtrsv_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen, ftnlen, 
	    ftnlen), dlartg_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *), dtrcon_(char *, char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen);


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

/*     To solve a system of equations in Hessenberg form with one */
/*     offdiagonal and one right-hand side. */

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

/*     LAMBDA  (input) DOUBLE PRECISION */
/*             This variable must contain the value to be added to the */
/*             diagonal elements of A. */

/*     D       (input/output) DOUBLE PRECISION array, dimension (M) */
/*             On entry, this array must contain the right-hand side */
/*             vector of the Hessenberg system. */
/*             On exit, if INFO = 0, this array contains the solution */
/*             vector of the Hessenberg system. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used to test for near singularity of */
/*             the triangular factor R of the Hessenberg matrix. A matrix */
/*             whose estimated condition number is less than 1/TOL is */
/*             considered to be nonsingular. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDDWOR,M+3) */
/*             The leading M-by-M part of this array is used for */
/*             computing the triangular factor of the QR decomposition */
/*             of the Hessenberg matrix. The remaining 3*M elements are */
/*             used as workspace for the computation of the reciprocal */
/*             condition estimate. */

/*     LDDWOR  INTEGER */
/*             The leading dimension of array DWORK.  LDDWOR >= MAX(1,M). */

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
/*     Supersedes Release 2.0 routine SB04BY by M. Vanbegin, and */
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
/*           ( LDDWOR.GE.MAX( 1, M ) ) */

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

#line 152 "SB04NY.f"
    /* Parameter adjustments */
#line 152 "SB04NY.f"
    a_dim1 = *lda;
#line 152 "SB04NY.f"
    a_offset = 1 + a_dim1;
#line 152 "SB04NY.f"
    a -= a_offset;
#line 152 "SB04NY.f"
    --d__;
#line 152 "SB04NY.f"
    --iwork;
#line 152 "SB04NY.f"
    dwork_dim1 = *lddwor;
#line 152 "SB04NY.f"
    dwork_offset = 1 + dwork_dim1;
#line 152 "SB04NY.f"
    dwork -= dwork_offset;
#line 152 "SB04NY.f"

#line 152 "SB04NY.f"
    /* Function Body */
#line 152 "SB04NY.f"
    *info = 0;

/*     For speed, no tests on the input scalar arguments are made. */
/*     Quick return if possible. */

#line 157 "SB04NY.f"
    if (*m == 0) {
#line 157 "SB04NY.f"
	return 0;
#line 157 "SB04NY.f"
    }

#line 160 "SB04NY.f"
    if (lsame_(ul, "U", (ftnlen)1, (ftnlen)1)) {

#line 162 "SB04NY.f"
	i__1 = *m;
#line 162 "SB04NY.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 163 "SB04NY.f"
	    i__3 = j + 1;
#line 163 "SB04NY.f"
	    i__2 = min(i__3,*m);
#line 163 "SB04NY.f"
	    dcopy_(&i__2, &a[j * a_dim1 + 1], &c__1, &dwork[j * dwork_dim1 + 
		    1], &c__1);
#line 164 "SB04NY.f"
	    dwork[j + j * dwork_dim1] += *lambda;
#line 165 "SB04NY.f"
/* L20: */
#line 165 "SB04NY.f"
	}

#line 167 "SB04NY.f"
	if (lsame_(rc, "R", (ftnlen)1, (ftnlen)1)) {
#line 168 "SB04NY.f"
	    *(unsigned char *)trans = 'N';

/*           A is an upper Hessenberg matrix, row transformations. */

#line 172 "SB04NY.f"
	    i__1 = *m - 1;
#line 172 "SB04NY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 173 "SB04NY.f"
		mj = *m - j;
#line 174 "SB04NY.f"
		if (dwork[j + 1 + j * dwork_dim1] != 0.) {
#line 175 "SB04NY.f"
		    dlartg_(&dwork[j + j * dwork_dim1], &dwork[j + 1 + j * 
			    dwork_dim1], &c__, &s, &r__);
#line 176 "SB04NY.f"
		    dwork[j + j * dwork_dim1] = r__;
#line 177 "SB04NY.f"
		    dwork[j + 1 + j * dwork_dim1] = 0.;
#line 178 "SB04NY.f"
		    drot_(&mj, &dwork[j + (j + 1) * dwork_dim1], lddwor, &
			    dwork[j + 1 + (j + 1) * dwork_dim1], lddwor, &c__,
			     &s);
#line 180 "SB04NY.f"
		    drot_(&c__1, &d__[j], &c__1, &d__[j + 1], &c__1, &c__, &s)
			    ;
#line 181 "SB04NY.f"
		}
#line 182 "SB04NY.f"
/* L40: */
#line 182 "SB04NY.f"
	    }

#line 184 "SB04NY.f"
	} else {
#line 185 "SB04NY.f"
	    *(unsigned char *)trans = 'T';

/*           A is an upper Hessenberg matrix, column transformations. */

#line 189 "SB04NY.f"
	    i__1 = *m - 1;
#line 189 "SB04NY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 190 "SB04NY.f"
		mj = *m - j;
#line 191 "SB04NY.f"
		if (dwork[mj + 1 + mj * dwork_dim1] != 0.) {
#line 192 "SB04NY.f"
		    dlartg_(&dwork[mj + 1 + (mj + 1) * dwork_dim1], &dwork[mj 
			    + 1 + mj * dwork_dim1], &c__, &s, &r__);
#line 194 "SB04NY.f"
		    dwork[mj + 1 + (mj + 1) * dwork_dim1] = r__;
#line 195 "SB04NY.f"
		    dwork[mj + 1 + mj * dwork_dim1] = 0.;
#line 196 "SB04NY.f"
		    drot_(&mj, &dwork[(mj + 1) * dwork_dim1 + 1], &c__1, &
			    dwork[mj * dwork_dim1 + 1], &c__1, &c__, &s);
#line 198 "SB04NY.f"
		    drot_(&c__1, &d__[mj + 1], &c__1, &d__[mj], &c__1, &c__, &
			    s);
#line 199 "SB04NY.f"
		}
#line 200 "SB04NY.f"
/* L60: */
#line 200 "SB04NY.f"
	    }

#line 202 "SB04NY.f"
	}
#line 203 "SB04NY.f"
    } else {

#line 205 "SB04NY.f"
	i__1 = *m;
#line 205 "SB04NY.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 206 "SB04NY.f"
	    i__2 = j - 1;
#line 206 "SB04NY.f"
	    j1 = max(i__2,1);
#line 207 "SB04NY.f"
	    i__2 = *m - j1 + 1;
#line 207 "SB04NY.f"
	    dcopy_(&i__2, &a[j1 + j * a_dim1], &c__1, &dwork[j1 + j * 
		    dwork_dim1], &c__1);
#line 208 "SB04NY.f"
	    dwork[j + j * dwork_dim1] += *lambda;
#line 209 "SB04NY.f"
/* L80: */
#line 209 "SB04NY.f"
	}

#line 211 "SB04NY.f"
	if (lsame_(rc, "R", (ftnlen)1, (ftnlen)1)) {
#line 212 "SB04NY.f"
	    *(unsigned char *)trans = 'N';

/*           A is a lower Hessenberg matrix, row transformations. */

#line 216 "SB04NY.f"
	    i__1 = *m - 1;
#line 216 "SB04NY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 217 "SB04NY.f"
		mj = *m - j;
#line 218 "SB04NY.f"
		if (dwork[mj + (mj + 1) * dwork_dim1] != 0.) {
#line 219 "SB04NY.f"
		    dlartg_(&dwork[mj + 1 + (mj + 1) * dwork_dim1], &dwork[mj 
			    + (mj + 1) * dwork_dim1], &c__, &s, &r__);
#line 221 "SB04NY.f"
		    dwork[mj + 1 + (mj + 1) * dwork_dim1] = r__;
#line 222 "SB04NY.f"
		    dwork[mj + (mj + 1) * dwork_dim1] = 0.;
#line 223 "SB04NY.f"
		    drot_(&mj, &dwork[mj + 1 + dwork_dim1], lddwor, &dwork[mj 
			    + dwork_dim1], lddwor, &c__, &s);
#line 225 "SB04NY.f"
		    drot_(&c__1, &d__[mj + 1], &c__1, &d__[mj], &c__1, &c__, &
			    s);
#line 226 "SB04NY.f"
		}
#line 227 "SB04NY.f"
/* L100: */
#line 227 "SB04NY.f"
	    }

#line 229 "SB04NY.f"
	} else {
#line 230 "SB04NY.f"
	    *(unsigned char *)trans = 'T';

/*           A is a lower Hessenberg matrix, column transformations. */

#line 234 "SB04NY.f"
	    i__1 = *m - 1;
#line 234 "SB04NY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 235 "SB04NY.f"
		mj = *m - j;
#line 236 "SB04NY.f"
		if (dwork[j + (j + 1) * dwork_dim1] != 0.) {
#line 237 "SB04NY.f"
		    dlartg_(&dwork[j + j * dwork_dim1], &dwork[j + (j + 1) * 
			    dwork_dim1], &c__, &s, &r__);
#line 238 "SB04NY.f"
		    dwork[j + j * dwork_dim1] = r__;
#line 239 "SB04NY.f"
		    dwork[j + (j + 1) * dwork_dim1] = 0.;
#line 240 "SB04NY.f"
		    drot_(&mj, &dwork[j + 1 + j * dwork_dim1], &c__1, &dwork[
			    j + 1 + (j + 1) * dwork_dim1], &c__1, &c__, &s);

#line 243 "SB04NY.f"
		    drot_(&c__1, &d__[j], &c__1, &d__[j + 1], &c__1, &c__, &s)
			    ;
#line 244 "SB04NY.f"
		}
#line 245 "SB04NY.f"
/* L120: */
#line 245 "SB04NY.f"
	    }

#line 247 "SB04NY.f"
	}
#line 248 "SB04NY.f"
    }

#line 250 "SB04NY.f"
    dtrcon_("1-norm", ul, "Non-unit", m, &dwork[dwork_offset], lddwor, &rcond,
	     &dwork[(*m + 1) * dwork_dim1 + 1], &iwork[1], info, (ftnlen)6, (
	    ftnlen)1, (ftnlen)8);
#line 252 "SB04NY.f"
    if (rcond <= *tol) {
#line 253 "SB04NY.f"
	*info = 1;
#line 254 "SB04NY.f"
    } else {
#line 255 "SB04NY.f"
	dtrsv_(ul, trans, "Non-unit", m, &dwork[dwork_offset], lddwor, &d__[1]
		, &c__1, (ftnlen)1, (ftnlen)1, (ftnlen)8);
#line 256 "SB04NY.f"
    }

#line 258 "SB04NY.f"
    return 0;
/* *** Last line of SB04NY *** */
} /* sb04ny_ */

