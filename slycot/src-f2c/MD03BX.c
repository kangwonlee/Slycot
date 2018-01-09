#line 1 "MD03BX.f"
/* MD03BX.f -- translated by f2c (version 20100827).
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

#line 1 "MD03BX.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int md03bx_(integer *m, integer *n, doublereal *fnorm, 
	doublereal *j, integer *ldj, doublereal *e, doublereal *jnorms, 
	doublereal *gnorm, integer *ipvt, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, l;
    static doublereal sum;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer itau;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dgeqp3_(integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), dlacpy_(char *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen), dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer wrkopt;


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

/*     To compute the QR factorization with column pivoting of an */
/*     m-by-n matrix J (m >= n), that is, J*P = Q*R, where Q is a matrix */
/*     with orthogonal columns, P a permutation matrix, and R an upper */
/*     trapezoidal matrix with diagonal elements of nonincreasing */
/*     magnitude, and to apply the transformation Q' on the error */
/*     vector e (in-situ). The 1-norm of the scaled gradient is also */
/*     returned. The matrix J could be the Jacobian of a nonlinear least */
/*     squares problem. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the Jacobian matrix J.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the Jacobian matrix J. */
/*             M >= N >= 0. */

/*     FNORM   (input) DOUBLE PRECISION */
/*             The Euclidean norm of the vector e.  FNORM >= 0. */

/*     J       (input/output) DOUBLE PRECISION array, dimension (LDJ, N) */
/*             On entry, the leading M-by-N part of this array must */
/*             contain the Jacobian matrix J. */
/*             On exit, the leading N-by-N upper triangular part of this */
/*             array contains the upper triangular factor R of the */
/*             Jacobian matrix. Note that for efficiency of the later */
/*             calculations, the matrix R is delivered with the leading */
/*             dimension MAX(1,N), possibly much smaller than the value */
/*             of LDJ on entry. */

/*     LDJ     (input/output) INTEGER */
/*             The leading dimension of array J. */
/*             On entry, LDJ >= MAX(1,M). */
/*             On exit,  LDJ >= MAX(1,N). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (M) */
/*             On entry, this array must contain the error vector e. */
/*             On exit, this array contains the updated vector Q'*e. */

/*     JNORMS  (output) DOUBLE PRECISION array, dimension (N) */
/*             This array contains the Euclidean norms of the columns of */
/*             the Jacobian matrix, considered in the initial order. */

/*     GNORM   (output) DOUBLE PRECISION */
/*             If FNORM > 0, the 1-norm of the scaled vector */
/*             J'*Q'*e/FNORM, with each element i further divided by */
/*             JNORMS(i) (if JNORMS(i) is nonzero). */
/*             If FNORM = 0, the returned value of GNORM is 0. */

/*     IPVT    (output) INTEGER array, dimension (N) */
/*             This array defines the permutation matrix P such that */
/*             J*P = Q*R. Column j of P is column IPVT(j) of the identity */
/*             matrix. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 1,      if N = 0 or M = 1; */
/*             LDWORK >= 4*N+1,  if N > 1. */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The algorithm uses QR factorization with column pivoting of the */
/*     matrix J, J*P = Q*R, and applies the orthogonal matrix Q' to the */
/*     vector e. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Elementary matrix operations, Jacobian matrix, matrix algebra, */
/*     matrix operations. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 143 "MD03BX.f"
    /* Parameter adjustments */
#line 143 "MD03BX.f"
    --dwork;
#line 143 "MD03BX.f"
    --ipvt;
#line 143 "MD03BX.f"
    --jnorms;
#line 143 "MD03BX.f"
    --e;
#line 143 "MD03BX.f"
    --j;
#line 143 "MD03BX.f"

#line 143 "MD03BX.f"
    /* Function Body */
#line 143 "MD03BX.f"
    *info = 0;
#line 144 "MD03BX.f"
    if (*m < 0) {
#line 145 "MD03BX.f"
	*info = -1;
#line 146 "MD03BX.f"
    } else if (*n < 0 || *m < *n) {
#line 147 "MD03BX.f"
	*info = -2;
#line 148 "MD03BX.f"
    } else if (*fnorm < 0.) {
#line 149 "MD03BX.f"
	*info = -3;
#line 150 "MD03BX.f"
    } else if (*ldj < max(1,*m)) {
#line 151 "MD03BX.f"
	*info = -5;
#line 152 "MD03BX.f"
    } else {
#line 153 "MD03BX.f"
	if (*n == 0 || *m == 1) {
#line 154 "MD03BX.f"
	    jwork = 1;
#line 155 "MD03BX.f"
	} else {
#line 156 "MD03BX.f"
	    jwork = (*n << 2) + 1;
#line 157 "MD03BX.f"
	}
#line 158 "MD03BX.f"
	if (*ldwork < jwork) {
#line 158 "MD03BX.f"
	    *info = -11;
#line 158 "MD03BX.f"
	}
#line 160 "MD03BX.f"
    }

#line 162 "MD03BX.f"
    if (*info != 0) {

/*        Error return. */

#line 166 "MD03BX.f"
	i__1 = -(*info);
#line 166 "MD03BX.f"
	xerbla_("MD03BX", &i__1, (ftnlen)6);
#line 167 "MD03BX.f"
	return 0;
#line 168 "MD03BX.f"
    }

/*     Quick return if possible. */

#line 172 "MD03BX.f"
    *gnorm = 0.;
#line 173 "MD03BX.f"
    if (*n == 0) {
#line 174 "MD03BX.f"
	*ldj = 1;
#line 175 "MD03BX.f"
	dwork[1] = 1.;
#line 176 "MD03BX.f"
	return 0;
#line 177 "MD03BX.f"
    } else if (*m == 1) {
#line 178 "MD03BX.f"
	jnorms[1] = abs(j[1]);
#line 179 "MD03BX.f"
	if (*fnorm * j[1] != 0.) {
#line 179 "MD03BX.f"
	    *gnorm = (d__1 = e[1] / *fnorm, abs(d__1));
#line 179 "MD03BX.f"
	}
#line 181 "MD03BX.f"
	*ldj = 1;
#line 182 "MD03BX.f"
	ipvt[1] = 1;
#line 183 "MD03BX.f"
	dwork[1] = 1.;
#line 184 "MD03BX.f"
	return 0;
#line 185 "MD03BX.f"
    }

/*     Initialize the column pivoting indices. */

#line 189 "MD03BX.f"
    i__1 = *n;
#line 189 "MD03BX.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 190 "MD03BX.f"
	ipvt[i__] = 0;
#line 191 "MD03BX.f"
/* L10: */
#line 191 "MD03BX.f"
    }

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 199 "MD03BX.f"
    itau = 1;
#line 200 "MD03BX.f"
    jwork = itau + *n;
#line 201 "MD03BX.f"
    wrkopt = 1;

/*     Compute the QR factorization with pivoting of J, and apply Q' to */
/*     the vector e. */

/*     Workspace: need:    4*N + 1; */
/*                prefer:  3*N + ( N+1 )*NB. */

#line 209 "MD03BX.f"
    i__1 = *ldwork - jwork + 1;
#line 209 "MD03BX.f"
    dgeqp3_(m, n, &j[1], ldj, &ipvt[1], &dwork[itau], &dwork[jwork], &i__1, 
	    info);
/* Computing MAX */
#line 211 "MD03BX.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 211 "MD03BX.f"
    wrkopt = max(i__1,i__2);

/*     Workspace: need:    N + 1; */
/*                prefer:  N + NB. */

#line 216 "MD03BX.f"
    i__1 = *ldwork - jwork + 1;
#line 216 "MD03BX.f"
    dormqr_("Left", "Transpose", m, &c__1, n, &j[1], ldj, &dwork[itau], &e[1],
	     m, &dwork[jwork], &i__1, info, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 218 "MD03BX.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 218 "MD03BX.f"
    wrkopt = max(i__1,i__2);

#line 220 "MD03BX.f"
    if (*ldj > *n) {

/*        Reshape the array J to have the leading dimension N. */
/*        This destroys the details of the orthogonal matrix Q. */

#line 225 "MD03BX.f"
	dlacpy_("Upper", n, n, &j[1], ldj, &j[1], n, (ftnlen)5);
#line 226 "MD03BX.f"
	*ldj = *n;
#line 227 "MD03BX.f"
    }

/*     Compute the norm of the scaled gradient and original column norms. */

#line 231 "MD03BX.f"
    if (*fnorm != 0.) {

#line 233 "MD03BX.f"
	i__1 = *n;
#line 233 "MD03BX.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 234 "MD03BX.f"
	    l = ipvt[i__];
#line 235 "MD03BX.f"
	    jnorms[l] = dnrm2_(&i__, &j[(i__ - 1) * *ldj + 1], &c__1);
#line 236 "MD03BX.f"
	    if (jnorms[l] != 0.) {
#line 237 "MD03BX.f"
		sum = ddot_(&i__, &j[(i__ - 1) * *ldj + 1], &c__1, &e[1], &
			c__1) / *fnorm;
/* Computing MAX */
#line 238 "MD03BX.f"
		d__2 = *gnorm, d__3 = (d__1 = sum / jnorms[l], abs(d__1));
#line 238 "MD03BX.f"
		*gnorm = max(d__2,d__3);
#line 239 "MD03BX.f"
	    }
#line 240 "MD03BX.f"
/* L20: */
#line 240 "MD03BX.f"
	}

#line 242 "MD03BX.f"
    } else {

#line 244 "MD03BX.f"
	i__1 = *n;
#line 244 "MD03BX.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 245 "MD03BX.f"
	    l = ipvt[i__];
#line 246 "MD03BX.f"
	    jnorms[l] = dnrm2_(&i__, &j[(i__ - 1) * *ldj + 1], &c__1);
#line 247 "MD03BX.f"
/* L30: */
#line 247 "MD03BX.f"
	}

#line 249 "MD03BX.f"
    }

#line 251 "MD03BX.f"
    dwork[1] = (doublereal) wrkopt;
#line 252 "MD03BX.f"
    return 0;

/* *** Last line of MD03BX *** */
} /* md03bx_ */

