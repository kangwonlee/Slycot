#line 1 "MB02TD.f"
/* MB02TD.f -- translated by f2c (version 20100827).
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

#line 1 "MB02TD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb02td_(char *norm, integer *n, doublereal *hnorm, 
	doublereal *h__, integer *ldh, integer *ipiv, doublereal *rcond, 
	integer *iwork, doublereal *dwork, integer *info, ftnlen norm_len)
{
    /* System generated locals */
    integer h_dim1, h_offset, i__1;
    doublereal d__1;

    /* Local variables */
    static integer j;
    static doublereal t;
    static integer jp, ix, kase, kase1;
    static doublereal scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlacon_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *);
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dlatrs_(
	    char *, char *, char *, char *, integer *, doublereal *, integer *
	    , doublereal *, doublereal *, doublereal *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static logical onenrm;
    static doublereal hinvnm;
    static char normin[1];
    static doublereal smlnum;


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

/*     To estimate the reciprocal of the condition number of an upper */
/*     Hessenberg matrix H, in either the 1-norm or the infinity-norm, */
/*     using the LU factorization computed by MB02SD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     NORM    CHARACTER*1 */
/*             Specifies whether the 1-norm condition number or the */
/*             infinity-norm condition number is required: */
/*             = '1' or 'O':  1-norm; */
/*             = 'I':         Infinity-norm. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix H.  N >= 0. */

/*     HNORM   (input) DOUBLE PRECISION */
/*             If NORM = '1' or 'O', the 1-norm of the original matrix H. */
/*             If NORM = 'I', the infinity-norm of the original matrix H. */

/*     H       (input) DOUBLE PRECISION array, dimension (LDH,N) */
/*             The factors L and U from the factorization H = P*L*U */
/*             as computed by MB02SD. */

/*     LDH     INTEGER */
/*             The leading dimension of the array H.  LDH >= max(1,N). */

/*     IPIV    (input) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= i <= N, row i of the matrix */
/*             was interchanged with row IPIV(i). */

/*     RCOND   (output) DOUBLE PRECISION */
/*             The reciprocal of the condition number of the matrix H, */
/*             computed as RCOND = 1/(norm(H) * norm(inv(H))). */

/*     Workspace */

/*     IWORK   DOUBLE PRECISION array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (3*N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     An estimate is obtained for norm(inv(H)), and the reciprocal of */
/*     the condition number is computed as */
/*        RCOND = 1 / ( norm(H) * norm(inv(H)) ). */

/*     REFERENCES */

/*     - */

/*     NUMERICAL ASPECTS */
/*                                2 */
/*     The algorithm requires 0( N ) operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Hessenberg form, matrix algebra. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */

/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input parameters. */

#line 137 "MB02TD.f"
    /* Parameter adjustments */
#line 137 "MB02TD.f"
    h_dim1 = *ldh;
#line 137 "MB02TD.f"
    h_offset = 1 + h_dim1;
#line 137 "MB02TD.f"
    h__ -= h_offset;
#line 137 "MB02TD.f"
    --ipiv;
#line 137 "MB02TD.f"
    --iwork;
#line 137 "MB02TD.f"
    --dwork;
#line 137 "MB02TD.f"

#line 137 "MB02TD.f"
    /* Function Body */
#line 137 "MB02TD.f"
    *info = 0;
#line 138 "MB02TD.f"
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O", (ftnlen)1, (
	    ftnlen)1);
#line 139 "MB02TD.f"
    if (! onenrm && ! lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {
#line 140 "MB02TD.f"
	*info = -1;
#line 141 "MB02TD.f"
    } else if (*n < 0) {
#line 142 "MB02TD.f"
	*info = -2;
#line 143 "MB02TD.f"
    } else if (*hnorm < 0.) {
#line 144 "MB02TD.f"
	*info = -3;
#line 145 "MB02TD.f"
    } else if (*ldh < max(1,*n)) {
#line 146 "MB02TD.f"
	*info = -5;
#line 147 "MB02TD.f"
    }
#line 148 "MB02TD.f"
    if (*info != 0) {
#line 149 "MB02TD.f"
	i__1 = -(*info);
#line 149 "MB02TD.f"
	xerbla_("MB02TD", &i__1, (ftnlen)6);
#line 150 "MB02TD.f"
	return 0;
#line 151 "MB02TD.f"
    }

/*     Quick return if possible. */

#line 155 "MB02TD.f"
    *rcond = 0.;
#line 156 "MB02TD.f"
    if (*n == 0) {
#line 157 "MB02TD.f"
	*rcond = 1.;
#line 158 "MB02TD.f"
	return 0;
#line 159 "MB02TD.f"
    } else if (*hnorm == 0.) {
#line 160 "MB02TD.f"
	return 0;
#line 161 "MB02TD.f"
    }

#line 163 "MB02TD.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12);

/*     Estimate the norm of inv(H). */

#line 167 "MB02TD.f"
    hinvnm = 0.;
#line 168 "MB02TD.f"
    *(unsigned char *)normin = 'N';
#line 169 "MB02TD.f"
    if (onenrm) {
#line 170 "MB02TD.f"
	kase1 = 1;
#line 171 "MB02TD.f"
    } else {
#line 172 "MB02TD.f"
	kase1 = 2;
#line 173 "MB02TD.f"
    }
#line 174 "MB02TD.f"
    kase = 0;
#line 175 "MB02TD.f"
L10:
#line 176 "MB02TD.f"
    dlacon_(n, &dwork[*n + 1], &dwork[1], &iwork[1], &hinvnm, &kase);
#line 177 "MB02TD.f"
    if (kase != 0) {
#line 178 "MB02TD.f"
	if (kase == kase1) {

/*           Multiply by inv(L). */

#line 182 "MB02TD.f"
	    i__1 = *n - 1;
#line 182 "MB02TD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 183 "MB02TD.f"
		jp = ipiv[j];
#line 184 "MB02TD.f"
		t = dwork[jp];
#line 185 "MB02TD.f"
		if (jp != j) {
#line 186 "MB02TD.f"
		    dwork[jp] = dwork[j];
#line 187 "MB02TD.f"
		    dwork[j] = t;
#line 188 "MB02TD.f"
		}
#line 189 "MB02TD.f"
		dwork[j + 1] -= t * h__[j + 1 + j * h_dim1];
#line 190 "MB02TD.f"
/* L20: */
#line 190 "MB02TD.f"
	    }

/*           Multiply by inv(U). */

#line 194 "MB02TD.f"
	    dlatrs_("Upper", "No transpose", "Non-unit", normin, n, &h__[
		    h_offset], ldh, &dwork[1], &scale, &dwork[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)12, (ftnlen)8, (ftnlen)1);
#line 196 "MB02TD.f"
	} else {

/*           Multiply by inv(U'). */

#line 200 "MB02TD.f"
	    dlatrs_("Upper", "Transpose", "Non-unit", normin, n, &h__[
		    h_offset], ldh, &dwork[1], &scale, &dwork[(*n << 1) + 1], 
		    info, (ftnlen)5, (ftnlen)9, (ftnlen)8, (ftnlen)1);

/*           Multiply by inv(L'). */

#line 205 "MB02TD.f"
	    for (j = *n - 1; j >= 1; --j) {
#line 206 "MB02TD.f"
		dwork[j] -= h__[j + 1 + j * h_dim1] * dwork[j + 1];
#line 207 "MB02TD.f"
		jp = ipiv[j];
#line 208 "MB02TD.f"
		if (jp != j) {
#line 209 "MB02TD.f"
		    t = dwork[jp];
#line 210 "MB02TD.f"
		    dwork[jp] = dwork[j];
#line 211 "MB02TD.f"
		    dwork[j] = t;
#line 212 "MB02TD.f"
		}
#line 213 "MB02TD.f"
/* L30: */
#line 213 "MB02TD.f"
	    }
#line 214 "MB02TD.f"
	}

/*        Divide X by 1/SCALE if doing so will not cause overflow. */

#line 218 "MB02TD.f"
	*(unsigned char *)normin = 'Y';
#line 219 "MB02TD.f"
	if (scale != 1.) {
#line 220 "MB02TD.f"
	    ix = idamax_(n, &dwork[1], &c__1);
#line 221 "MB02TD.f"
	    if (scale < (d__1 = dwork[ix], abs(d__1)) * smlnum || scale == 0.)
		     {
#line 221 "MB02TD.f"
		goto L40;
#line 221 "MB02TD.f"
	    }
#line 223 "MB02TD.f"
	    drscl_(n, &scale, &dwork[1], &c__1);
#line 224 "MB02TD.f"
	}
#line 225 "MB02TD.f"
	goto L10;
#line 226 "MB02TD.f"
    }

/*     Compute the estimate of the reciprocal condition number. */

#line 230 "MB02TD.f"
    if (hinvnm != 0.) {
#line 230 "MB02TD.f"
	*rcond = 1. / hinvnm / *hnorm;
#line 230 "MB02TD.f"
    }

#line 233 "MB02TD.f"
L40:
#line 234 "MB02TD.f"
    return 0;
/* *** Last line of MB02TD *** */
} /* mb02td_ */

