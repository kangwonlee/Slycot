#line 1 "MB04XY.f"
/* MB04XY.f -- translated by f2c (version 20100827).
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

#line 1 "MB04XY.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb04xy_(char *jobu, char *jobv, integer *m, integer *n, 
	doublereal *x, integer *ldx, doublereal *taup, doublereal *tauq, 
	doublereal *u, integer *ldu, doublereal *v, integer *ldv, logical *
	inul, integer *info, ftnlen jobu_len, ftnlen jobv_len)
{
    /* System generated locals */
    integer u_dim1, u_offset, v_dim1, v_offset, x_dim1, x_offset, i__1, i__2;

    /* Local variables */
    static integer i__, l, p, im, ioff, ncol;
    extern /* Subroutine */ int dlarf_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal dwork[1], first;
    static logical wantu, wantv, ljobua, ljobva;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static logical ljobus, ljobvs;


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

/*     To apply the Householder transformations Pj stored in factored */
/*     form into the columns of the array X, to the desired columns of */
/*     the matrix U by premultiplication, and/or the Householder */
/*     transformations Qj stored in factored form into the rows of the */
/*     array X, to the desired columns of the matrix V by */
/*     premultiplication. The Householder transformations Pj and Qj */
/*     are stored as produced by LAPACK Library routine DGEBRD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBU    CHARACTER*1 */
/*             Specifies whether to transform the columns in U as */
/*             follows: */
/*             = 'N':  Do not transform the columns in U; */
/*             = 'A':  Transform the columns in U (U has M columns); */
/*             = 'S':  Transform the columns in U (U has min(M,N) */
/*                     columns). */

/*     JOBV    CHARACTER*1 */
/*             Specifies whether to transform the columns in V as */
/*             follows: */
/*             = 'N':  Do not transform the columns in V; */
/*             = 'A':  Transform the columns in V (V has N columns); */
/*             = 'S':  Transform the columns in V (V has min(M,N) */
/*                     columns). */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix X.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix X.  N >= 0. */

/*     X       (input) DOUBLE PRECISION array, dimension (LDX,N) */
/*             The leading M-by-N part contains in the columns of its */
/*             lower triangle the Householder transformations Pj, and */
/*             in the rows of its upper triangle the Householder */
/*             transformations Qj in factored form. */
/*             X is modified by the routine but restored on exit. */

/*     LDX     INTEGER */
/*             The leading dimension of the array X.   LDX >= MAX(1,M). */

/*     TAUP    (input) DOUBLE PRECISION array, dimension (MIN(M,N)) */
/*             The scalar factors of the Householder transformations Pj. */

/*     TAUQ    (input) DOUBLE PRECISION array, dimension (MIN(M,N)) */
/*             The scalar factors of the Householder transformations Qj. */

/*     U       (input/output) DOUBLE PRECISION array, dimension (LDU,*) */
/*             On entry, U contains the M-by-M (if JOBU = 'A') or */
/*             M-by-min(M,N) (if JOBU = 'S') matrix U. */
/*             On exit, the Householder transformations Pj have been */
/*             applied to each column i of U corresponding to a parameter */
/*             INUL(i) = .TRUE. */
/*             NOTE that U is not referenced if JOBU = 'N'. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= MAX(1,M), if JOBU = 'A' or JOBU = 'S'; */
/*             LDU >= 1,        if JOBU = 'N'. */

/*     V       (input/output) DOUBLE PRECISION array, dimension (LDV,*) */
/*             On entry, V contains the N-by-N (if JOBV = 'A') or */
/*             N-by-min(M,N) (if JOBV = 'S') matrix V. */
/*             On exit, the Householder transformations Qj have been */
/*             applied to each column i of V corresponding to a parameter */
/*             INUL(i) = .TRUE. */
/*             NOTE that V is not referenced if JOBV = 'N'. */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= MAX(1,M), if JOBV = 'A' or JOBV = 'S'; */
/*             LDV >= 1,        if JOBV = 'N'. */

/*     INUL    (input) LOGICAL array, dimension (MAX(M,N)) */
/*             INUL(i) = .TRUE. if the i-th column of U and/or V is to be */
/*             transformed, and INUL(i) = .FALSE., otherwise. */
/*             (1 <= i <= MAX(M,N)). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The Householder transformations Pj or Qj are applied to the */
/*     columns of U or V indexed by I for which INUL(I) = .TRUE.. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1997. */
/*     Supersedes Release 2.0 routine MB04PZ by S. Van Huffel, Katholieke */
/*     University Leuven, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Bidiagonalization, orthogonal transformation, singular subspace, */
/*     singular value decomposition. */

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

#line 165 "MB04XY.f"
    /* Parameter adjustments */
#line 165 "MB04XY.f"
    x_dim1 = *ldx;
#line 165 "MB04XY.f"
    x_offset = 1 + x_dim1;
#line 165 "MB04XY.f"
    x -= x_offset;
#line 165 "MB04XY.f"
    --taup;
#line 165 "MB04XY.f"
    --tauq;
#line 165 "MB04XY.f"
    u_dim1 = *ldu;
#line 165 "MB04XY.f"
    u_offset = 1 + u_dim1;
#line 165 "MB04XY.f"
    u -= u_offset;
#line 165 "MB04XY.f"
    v_dim1 = *ldv;
#line 165 "MB04XY.f"
    v_offset = 1 + v_dim1;
#line 165 "MB04XY.f"
    v -= v_offset;
#line 165 "MB04XY.f"
    --inul;
#line 165 "MB04XY.f"

#line 165 "MB04XY.f"
    /* Function Body */
#line 165 "MB04XY.f"
    *info = 0;
#line 166 "MB04XY.f"
    ljobua = lsame_(jobu, "A", (ftnlen)1, (ftnlen)1);
#line 167 "MB04XY.f"
    ljobus = lsame_(jobu, "S", (ftnlen)1, (ftnlen)1);
#line 168 "MB04XY.f"
    ljobva = lsame_(jobv, "A", (ftnlen)1, (ftnlen)1);
#line 169 "MB04XY.f"
    ljobvs = lsame_(jobv, "S", (ftnlen)1, (ftnlen)1);
#line 170 "MB04XY.f"
    wantu = ljobua || ljobus;
#line 171 "MB04XY.f"
    wantv = ljobva || ljobvs;

/*     Test the input scalar arguments. */

#line 175 "MB04XY.f"
    if (! wantu && ! lsame_(jobu, "N", (ftnlen)1, (ftnlen)1)) {
#line 176 "MB04XY.f"
	*info = -1;
#line 177 "MB04XY.f"
    } else if (! wantv && ! lsame_(jobv, "N", (ftnlen)1, (ftnlen)1)) {
#line 178 "MB04XY.f"
	*info = -2;
#line 179 "MB04XY.f"
    } else if (*m < 0) {
#line 180 "MB04XY.f"
	*info = -3;
#line 181 "MB04XY.f"
    } else if (*n < 0) {
#line 182 "MB04XY.f"
	*info = -4;
#line 183 "MB04XY.f"
    } else if (*ldx < max(1,*m)) {
#line 184 "MB04XY.f"
	*info = -6;
#line 185 "MB04XY.f"
    } else if (wantu && *ldu < max(1,*m) || ! wantu && *ldu < 1) {
#line 187 "MB04XY.f"
	*info = -10;
#line 188 "MB04XY.f"
    } else if (wantv && *ldv < max(1,*n) || ! wantv && *ldv < 1) {
#line 190 "MB04XY.f"
	*info = -12;
#line 191 "MB04XY.f"
    }

#line 193 "MB04XY.f"
    if (*info != 0) {

/*        Error return */

#line 197 "MB04XY.f"
	i__1 = -(*info);
#line 197 "MB04XY.f"
	xerbla_("MB04XY", &i__1, (ftnlen)6);
#line 198 "MB04XY.f"
	return 0;
#line 199 "MB04XY.f"
    }

/*     Quick return if possible. */

#line 203 "MB04XY.f"
    p = min(*m,*n);
#line 204 "MB04XY.f"
    if (p == 0) {
#line 204 "MB04XY.f"
	return 0;
#line 204 "MB04XY.f"
    }

#line 207 "MB04XY.f"
    if (*m < *n) {
#line 208 "MB04XY.f"
	ioff = 1;
#line 209 "MB04XY.f"
    } else {
#line 210 "MB04XY.f"
	ioff = 0;
#line 211 "MB04XY.f"
    }

/*     Apply the Householder transformations Pj onto the desired */
/*     columns of U. */

/* Computing MIN */
#line 216 "MB04XY.f"
    i__1 = *m - 1;
#line 216 "MB04XY.f"
    im = min(i__1,*n);
#line 217 "MB04XY.f"
    if (wantu && im > 0) {
#line 218 "MB04XY.f"
	if (ljobua) {
#line 219 "MB04XY.f"
	    ncol = *m;
#line 220 "MB04XY.f"
	} else {
#line 221 "MB04XY.f"
	    ncol = p;
#line 222 "MB04XY.f"
	}

#line 224 "MB04XY.f"
	i__1 = ncol;
#line 224 "MB04XY.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 225 "MB04XY.f"
	    if (inul[i__]) {

#line 227 "MB04XY.f"
		for (l = im; l >= 1; --l) {
#line 228 "MB04XY.f"
		    if (taup[l] != 0.) {
#line 229 "MB04XY.f"
			first = x[l + ioff + l * x_dim1];
#line 230 "MB04XY.f"
			x[l + ioff + l * x_dim1] = 1.;
#line 231 "MB04XY.f"
			i__2 = *m - l + 1 - ioff;
#line 231 "MB04XY.f"
			dlarf_("Left", &i__2, &c__1, &x[l + ioff + l * x_dim1]
				, &c__1, &taup[l], &u[l + ioff + i__ * u_dim1]
				, ldu, dwork, (ftnlen)4);
#line 233 "MB04XY.f"
			x[l + ioff + l * x_dim1] = first;
#line 234 "MB04XY.f"
		    }
#line 235 "MB04XY.f"
/* L20: */
#line 235 "MB04XY.f"
		}

#line 237 "MB04XY.f"
	    }
#line 238 "MB04XY.f"
/* L40: */
#line 238 "MB04XY.f"
	}

#line 240 "MB04XY.f"
    }

/*     Apply the Householder transformations Qj onto the desired columns */
/*     of V. */

/* Computing MIN */
#line 245 "MB04XY.f"
    i__1 = *n - 1;
#line 245 "MB04XY.f"
    im = min(i__1,*m);
#line 246 "MB04XY.f"
    if (wantv && im > 0) {
#line 247 "MB04XY.f"
	if (ljobva) {
#line 248 "MB04XY.f"
	    ncol = *n;
#line 249 "MB04XY.f"
	} else {
#line 250 "MB04XY.f"
	    ncol = p;
#line 251 "MB04XY.f"
	}

#line 253 "MB04XY.f"
	i__1 = ncol;
#line 253 "MB04XY.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 254 "MB04XY.f"
	    if (inul[i__]) {

#line 256 "MB04XY.f"
		for (l = im; l >= 1; --l) {
#line 257 "MB04XY.f"
		    if (tauq[l] != 0.) {
#line 258 "MB04XY.f"
			first = x[l + (l + 1 - ioff) * x_dim1];
#line 259 "MB04XY.f"
			x[l + (l + 1 - ioff) * x_dim1] = 1.;
#line 260 "MB04XY.f"
			i__2 = *n - l + ioff;
#line 260 "MB04XY.f"
			dlarf_("Left", &i__2, &c__1, &x[l + (l + 1 - ioff) * 
				x_dim1], ldx, &tauq[l], &v[l + 1 - ioff + i__ 
				* v_dim1], ldv, dwork, (ftnlen)4);
#line 263 "MB04XY.f"
			x[l + (l + 1 - ioff) * x_dim1] = first;
#line 264 "MB04XY.f"
		    }
#line 265 "MB04XY.f"
/* L60: */
#line 265 "MB04XY.f"
		}

#line 267 "MB04XY.f"
	    }
#line 268 "MB04XY.f"
/* L80: */
#line 268 "MB04XY.f"
	}

#line 270 "MB04XY.f"
    }

#line 272 "MB04XY.f"
    return 0;
/* *** Last line of MB04XY *** */
} /* mb04xy_ */

