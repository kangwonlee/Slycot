#line 1 "SB01FY.f"
/* SB01FY.f -- translated by f2c (version 20100827).
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

#line 1 "SB01FY.f"
/* Table of constant values */

static integer c__1 = 1;
static logical c_false = FALSE_;
static integer c_n1 = -1;
static integer c__2 = 2;
static doublereal c_b15 = 0.;
static doublereal c_b16 = 1.;

/* Subroutine */ int sb01fy_(logical *discr, integer *n, integer *m, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	f, integer *ldf, doublereal *v, integer *ldv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, f_dim1, f_offset, v_dim1, 
	    v_offset, i__1;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__;
    static doublereal u[4]	/* was [2][2] */, r11, r12, cs, r22, at[4]	
	    /* was [2][2] */, sn, temp;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal scale;
    extern /* Subroutine */ int mb04ox_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), drotg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *), sb03oy_(logical *, logical *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal dummy[4]	/* was [2][2] */;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlapy3_(doublereal 
	    *, doublereal *, doublereal *);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlatzm_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen), dtrtri_(char *, char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen);


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

/*     To compute the inner denominator of a right-coprime factorization */
/*     of a system of order N, where N is either 1 or 2. Specifically, */
/*     given the N-by-N unstable system state matrix A and the N-by-M */
/*     system input matrix B, an M-by-N state-feedback matrix F and */
/*     an M-by-M matrix V are constructed, such that the system */
/*     (A + B*F, B*V, F, V) is inner. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DISCR   LOGICAL */
/*             Specifies the type of system as follows: */
/*             = .FALSE.:  continuous-time system; */
/*             = .TRUE. :  discrete-time system. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A and also the number of rows of */
/*             the matrix B and the number of columns of the matrix F. */
/*             N is either 1 or 2. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrices B and V, and also */
/*             the number of rows of the matrix F.  M >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             system state matrix A whose eigenvalues must have positive */
/*             real parts if DISCR = .FALSE. or moduli greater than unity */
/*             if DISCR = .TRUE.. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= N. */

/*     F       (output) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading M-by-N part of this array contains the state- */
/*             feedback matrix F which assigns one eigenvalue (if N = 1) */
/*             or two eigenvalues (if N = 2) of the matrix A + B*F in */
/*             symmetric positions with respect to the imaginary axis */
/*             (if DISCR = .FALSE.) or the unit circle (if */
/*             DISCR = .TRUE.). */

/*     LDF     INTEGER */
/*             The leading dimension of array F.  LDF >= MAX(1,M). */

/*     V       (output) DOUBLE PRECISION array, dimension (LDV,M) */
/*             The leading M-by-M upper triangular part of this array */
/*             contains the input/output matrix V of the resulting inner */
/*             system in upper triangular form. */
/*             If DISCR = .FALSE., the resulting V is an identity matrix. */

/*     LDV     INTEGER */
/*             The leading dimension of array V.  LDF >= MAX(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if uncontrollability of the pair (A,B) is detected; */
/*             = 2:  if A is stable or at the stability limit; */
/*             = 3:  if N = 2 and A has a pair of real eigenvalues. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine RCFID2. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     Feb. 1999, A. Varga, DLR Oberpfaffenhofen. */

/*     KEYWORDS */

/*     Coprime factorization, eigenvalue, eigenvalue assignment, */
/*     feedback control, pole placement, state-space model. */

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

/*     For efficiency reasons, the parameters are not checked. */

#line 139 "SB01FY.f"
    /* Parameter adjustments */
#line 139 "SB01FY.f"
    a_dim1 = *lda;
#line 139 "SB01FY.f"
    a_offset = 1 + a_dim1;
#line 139 "SB01FY.f"
    a -= a_offset;
#line 139 "SB01FY.f"
    b_dim1 = *ldb;
#line 139 "SB01FY.f"
    b_offset = 1 + b_dim1;
#line 139 "SB01FY.f"
    b -= b_offset;
#line 139 "SB01FY.f"
    f_dim1 = *ldf;
#line 139 "SB01FY.f"
    f_offset = 1 + f_dim1;
#line 139 "SB01FY.f"
    f -= f_offset;
#line 139 "SB01FY.f"
    v_dim1 = *ldv;
#line 139 "SB01FY.f"
    v_offset = 1 + v_dim1;
#line 139 "SB01FY.f"
    v -= v_offset;
#line 139 "SB01FY.f"

#line 139 "SB01FY.f"
    /* Function Body */
#line 139 "SB01FY.f"
    *info = 0;

/*     Compute an N-by-N upper triangular R such that R'*R = B*B' and */
/*     find an upper triangular matrix U in the equation */

/*     A'*U'*U + U'*U*A = R'*R if DISCR = .FALSE. or */
/*     A'*U'*U*A - U'*U = R'*R if DISCR = .TRUE. . */

#line 147 "SB01FY.f"
    ma02ad_("Full", n, m, &b[b_offset], ldb, &f[f_offset], ldf, (ftnlen)4);

#line 149 "SB01FY.f"
    if (*n == 1) {

/*        The N = 1 case. */

#line 153 "SB01FY.f"
	if (*m > 1) {
#line 153 "SB01FY.f"
	    dlarfg_(m, &f[f_dim1 + 1], &f[f_dim1 + 2], &c__1, &temp);
#line 153 "SB01FY.f"
	}
#line 155 "SB01FY.f"
	r11 = (d__1 = f[f_dim1 + 1], abs(d__1));

/*        Make sure A is unstable or divergent and find U. */

#line 159 "SB01FY.f"
	if (*discr) {
#line 160 "SB01FY.f"
	    temp = (d__1 = a[a_dim1 + 1], abs(d__1));
#line 161 "SB01FY.f"
	    if (temp <= 1.) {
#line 162 "SB01FY.f"
		*info = 2;
#line 163 "SB01FY.f"
		return 0;
#line 164 "SB01FY.f"
	    } else {
#line 165 "SB01FY.f"
		temp = r11 / sqrt((temp - 1.) * (temp + 1.));
#line 166 "SB01FY.f"
	    }
#line 167 "SB01FY.f"
	} else {
#line 168 "SB01FY.f"
	    if (a[a_dim1 + 1] <= 0.) {
#line 169 "SB01FY.f"
		*info = 2;
#line 170 "SB01FY.f"
		return 0;
#line 171 "SB01FY.f"
	    } else {
#line 172 "SB01FY.f"
		temp = r11 / sqrt((d__1 = a[a_dim1 + 1] * 2., abs(d__1)));
#line 173 "SB01FY.f"
	    }
#line 174 "SB01FY.f"
	}
#line 175 "SB01FY.f"
	u[0] = temp;
#line 176 "SB01FY.f"
	scale = 1.;
#line 177 "SB01FY.f"
    } else {

/*        The N = 2 case. */

#line 181 "SB01FY.f"
	if (*m > 1) {
#line 182 "SB01FY.f"
	    dlarfg_(m, &f[f_dim1 + 1], &f[f_dim1 + 2], &c__1, &temp);
#line 183 "SB01FY.f"
	    i__1 = *n - 1;
#line 183 "SB01FY.f"
	    dlatzm_("Left", m, &i__1, &f[f_dim1 + 2], &c__1, &temp, &f[(
		    f_dim1 << 1) + 1], &f[(f_dim1 << 1) + 2], ldf, &v[
		    v_offset], (ftnlen)4);
#line 185 "SB01FY.f"
	}
#line 186 "SB01FY.f"
	r11 = f[f_dim1 + 1];
#line 187 "SB01FY.f"
	r12 = f[(f_dim1 << 1) + 1];
#line 188 "SB01FY.f"
	if (*m > 2) {
#line 188 "SB01FY.f"
	    i__1 = *m - 1;
#line 188 "SB01FY.f"
	    dlarfg_(&i__1, &f[(f_dim1 << 1) + 2], &f[(f_dim1 << 1) + 3], &
		    c__1, &temp);
#line 188 "SB01FY.f"
	}
#line 190 "SB01FY.f"
	if (*m == 1) {
#line 191 "SB01FY.f"
	    r22 = 0.;
#line 192 "SB01FY.f"
	} else {
#line 193 "SB01FY.f"
	    r22 = f[(f_dim1 << 1) + 2];
#line 194 "SB01FY.f"
	}
#line 195 "SB01FY.f"
	at[0] = a[a_dim1 + 1];
#line 196 "SB01FY.f"
	at[2] = a[a_dim1 + 2];
#line 197 "SB01FY.f"
	at[1] = a[(a_dim1 << 1) + 1];
#line 198 "SB01FY.f"
	at[3] = a[(a_dim1 << 1) + 2];
#line 199 "SB01FY.f"
	u[0] = r11;
#line 200 "SB01FY.f"
	u[2] = r12;
#line 201 "SB01FY.f"
	u[3] = r22;
#line 202 "SB01FY.f"
	sb03oy_(discr, &c_false, &c_n1, at, &c__2, u, &c__2, dummy, &c__2, &
		scale, info);
#line 204 "SB01FY.f"
	if (*info != 0) {
#line 205 "SB01FY.f"
	    if (*info != 4) {
#line 206 "SB01FY.f"
		*info = 2;
#line 207 "SB01FY.f"
	    } else {
#line 208 "SB01FY.f"
		*info = 3;
#line 209 "SB01FY.f"
	    }
#line 210 "SB01FY.f"
	    return 0;
#line 211 "SB01FY.f"
	}
#line 212 "SB01FY.f"
    }

/*     Check the controllability of the pair (A,B). */

/*     Warning. Only an exact controllability check is performed. */
/*              If the pair (A,B) is nearly uncontrollable, then */
/*              the computed results may be inaccurate. */

#line 220 "SB01FY.f"
    i__1 = *n;
#line 220 "SB01FY.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 221 "SB01FY.f"
	if (u[i__ + (i__ << 1) - 3] == 0.) {
#line 222 "SB01FY.f"
	    *info = 1;
#line 223 "SB01FY.f"
	    return 0;
#line 224 "SB01FY.f"
	}
#line 225 "SB01FY.f"
/* L10: */
#line 225 "SB01FY.f"
    }

/*     Set V = I. */

#line 229 "SB01FY.f"
    dlaset_("Upper", m, m, &c_b15, &c_b16, &v[v_offset], ldv, (ftnlen)5);

#line 231 "SB01FY.f"
    if (*discr) {

/*        Compute an upper triangular matrix V such that */
/*                                 -1 */
/*        V*V' = (I+B'*inv(U'*U)*B)  . */

/*        First compute F = B'*inv(U) and the Cholesky factorization */
/*        of I + F*F'. */

#line 240 "SB01FY.f"
	i__1 = *m;
#line 240 "SB01FY.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 241 "SB01FY.f"
	    f[i__ + f_dim1] = b[i__ * b_dim1 + 1] / u[0] * scale;
#line 242 "SB01FY.f"
/* L20: */
#line 242 "SB01FY.f"
	}
#line 243 "SB01FY.f"
	if (*n == 2) {
#line 244 "SB01FY.f"
	    i__1 = *m;
#line 244 "SB01FY.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 245 "SB01FY.f"
		f[i__ + (f_dim1 << 1)] = (b[i__ * b_dim1 + 2] - f[i__ + 
			f_dim1] * u[2]) / u[3] * scale;
#line 246 "SB01FY.f"
/* L30: */
#line 246 "SB01FY.f"
	    }
#line 247 "SB01FY.f"
	    mb04ox_(m, &v[v_offset], ldv, &f[(f_dim1 << 1) + 1], &c__1);
#line 248 "SB01FY.f"
	}
#line 249 "SB01FY.f"
	mb04ox_(m, &v[v_offset], ldv, &f[f_dim1 + 1], &c__1);
#line 250 "SB01FY.f"
	dtrtri_("Upper", "NonUnit", m, &v[v_offset], ldv, info, (ftnlen)5, (
		ftnlen)7);
#line 251 "SB01FY.f"
    }

/*     Compute the feedback matrix F as: */

/*     1)   If DISCR = .FALSE. */

/*             F = -B'*inv(U'*U); */

/*     2)   If DISCR = .TRUE. */
/*                                -1 */
/*             F = -B'*(U'*U+B*B')  *A. */

#line 263 "SB01FY.f"
    if (*n == 1) {
#line 264 "SB01FY.f"
	if (*discr) {
#line 265 "SB01FY.f"
	    temp = -a[a_dim1 + 1];
#line 266 "SB01FY.f"
	    r11 = dlapy2_(u, &r11);
#line 267 "SB01FY.f"
	    i__1 = *m;
#line 267 "SB01FY.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 268 "SB01FY.f"
		f[i__ + f_dim1] = b[i__ * b_dim1 + 1] / r11 / r11 * temp;
#line 269 "SB01FY.f"
/* L40: */
#line 269 "SB01FY.f"
	    }
#line 270 "SB01FY.f"
	} else {
#line 271 "SB01FY.f"
	    r11 = u[0];
#line 272 "SB01FY.f"
	    i__1 = *m;
#line 272 "SB01FY.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 273 "SB01FY.f"
		f[i__ + f_dim1] = -(b[i__ * b_dim1 + 1] / r11 / r11);
#line 274 "SB01FY.f"
/* L50: */
#line 274 "SB01FY.f"
	    }
#line 275 "SB01FY.f"
	}
#line 276 "SB01FY.f"
    } else {

/*        Set R = U  if DISCR = .FALSE. or compute the Cholesky */
/*        factorization of R'*R = U'*U+B*B' if DISCR = .TRUE.. */

#line 281 "SB01FY.f"
	if (*discr) {
#line 282 "SB01FY.f"
	    temp = u[0];
#line 283 "SB01FY.f"
	    drotg_(&r11, &temp, &cs, &sn);
#line 284 "SB01FY.f"
	    temp = -sn * r12 + cs * u[2];
#line 285 "SB01FY.f"
	    r12 = cs * r12 + sn * u[2];
#line 286 "SB01FY.f"
	    r22 = dlapy3_(&r22, &temp, &u[3]);
#line 287 "SB01FY.f"
	} else {
#line 288 "SB01FY.f"
	    r11 = u[0];
#line 289 "SB01FY.f"
	    r12 = u[2];
#line 290 "SB01FY.f"
	    r22 = u[3];
#line 291 "SB01FY.f"
	}

/*        Compute F = -B'*inv(R'*R). */

#line 295 "SB01FY.f"
	i__1 = *m;
#line 295 "SB01FY.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 296 "SB01FY.f"
	    f[i__ + f_dim1] = -b[i__ * b_dim1 + 1] / r11;
#line 297 "SB01FY.f"
	    f[i__ + (f_dim1 << 1)] = -(b[i__ * b_dim1 + 2] + f[i__ + f_dim1] *
		     r12) / r22;
#line 298 "SB01FY.f"
	    f[i__ + (f_dim1 << 1)] /= r22;
#line 299 "SB01FY.f"
	    f[i__ + f_dim1] = (f[i__ + f_dim1] - f[i__ + (f_dim1 << 1)] * r12)
		     / r11;
#line 300 "SB01FY.f"
/* L60: */
#line 300 "SB01FY.f"
	}
#line 301 "SB01FY.f"
	if (*discr) {

/*           Compute F <-- F*A. */

#line 305 "SB01FY.f"
	    i__1 = *m;
#line 305 "SB01FY.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 306 "SB01FY.f"
		temp = f[i__ + f_dim1] * a[a_dim1 + 1] + f[i__ + (f_dim1 << 1)
			] * a[a_dim1 + 2];
#line 307 "SB01FY.f"
		f[i__ + (f_dim1 << 1)] = f[i__ + f_dim1] * a[(a_dim1 << 1) + 
			1] + f[i__ + (f_dim1 << 1)] * a[(a_dim1 << 1) + 2];
#line 308 "SB01FY.f"
		f[i__ + f_dim1] = temp;
#line 309 "SB01FY.f"
/* L70: */
#line 309 "SB01FY.f"
	    }
#line 310 "SB01FY.f"
	}
#line 311 "SB01FY.f"
    }

#line 313 "SB01FY.f"
    return 0;
/* *** Last line of SB01FY *** */
} /* sb01fy_ */

