#line 1 "SB01BY.f"
/* SB01BY.f -- translated by f2c (version 20100827).
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

#line 1 "SB01BY.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b4 = 0.;
static integer c__2 = 2;

/* Subroutine */ int sb01by_(integer *n, integer *m, doublereal *s, 
	doublereal *p, doublereal *a, doublereal *b, doublereal *f, 
	doublereal *tol, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, f_dim1, f_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__;
    static integer j;
    static doublereal r__, x, y, z__, b1, b2, c0, c1, c3, c4, b21, c11, c12, 
	    c21, c22, cs, s12, cu, cv, s21;
    static integer ir;
    static doublereal rn, sn, wi, su, sv, wr, dc0, dc2, dc3, wi1, wr1, sig, 
	    tau1, tau2, absr;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal diffr;
    extern doublereal dlamc3_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlasv2_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlatzm_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen);


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

/*     To solve an N-by-N pole placement problem for the simple cases */
/*     N = 1 or N = 2: given the N-by-N matrix A and N-by-M matrix B, */
/*     construct an M-by-N matrix F such that A + B*F has prescribed */
/*     eigenvalues. These eigenvalues are specified by their sum S and */
/*     product P (if N = 2). The resulting F has minimum Frobenius norm. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A and also the number of rows of */
/*             the matrix B and the number of columns of the matrix F. */
/*             N is either 1, if a single real eigenvalue is prescribed */
/*             or 2, if a complex conjugate pair or a set of two real */
/*             eigenvalues are prescribed. */

/*     M       (input) INTEGER */
/*             The number of columns of the matrix B and also the number */
/*             of rows of the matrix F.  M >= 1. */

/*     S       (input) DOUBLE PRECISION */
/*             The sum of the prescribed eigenvalues if N = 2 or the */
/*             value of prescribed eigenvalue if N = 1. */

/*     P       (input) DOUBLE PRECISION */
/*             The product of the prescribed eigenvalues if N = 2. */
/*             Not referenced if N = 1. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (N,N) */
/*             On entry, this array must contain the N-by-N state */
/*             dynamics matrix whose eigenvalues have to be moved to */
/*             prescribed locations. */
/*             On exit, this array contains no useful information. */

/*     B       (input/output) DOUBLE PRECISION array, dimension (N,M) */
/*             On entry, this array must contain the N-by-M input/state */
/*             matrix B. */
/*             On exit, this array contains no useful information. */

/*     F       (output) DOUBLE PRECISION array, dimension (M,N) */
/*             The state feedback matrix F which assigns one pole or two */
/*             poles of the closed-loop matrix A + B*F. */
/*             If N = 2 and the pair (A,B) is not controllable */
/*             (INFO = 1), then F(1,1) and F(1,2) contain the elements of */
/*             an orthogonal rotation which can be used to remove the */
/*             uncontrollable part of the pair (A,B). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The absolute tolerance level below which the elements of A */
/*             and B are considered zero (used for controllability test). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (M) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  if uncontrollability of the pair (A,B) is detected. */

/*     CONTRIBUTOR */

/*     A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine SB01BY. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     May  2003, A. Varga, German Aerospace Center. */

/*     KEYWORDS */

/*     Eigenvalue, eigenvalue assignment, feedback control, pole */
/*     placement, state-space model. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     For efficiency reasons, the parameters are not checked. */

#line 132 "SB01BY.f"
    /* Parameter adjustments */
#line 132 "SB01BY.f"
    b_dim1 = *n;
#line 132 "SB01BY.f"
    b_offset = 1 + b_dim1;
#line 132 "SB01BY.f"
    b -= b_offset;
#line 132 "SB01BY.f"
    a_dim1 = *n;
#line 132 "SB01BY.f"
    a_offset = 1 + a_dim1;
#line 132 "SB01BY.f"
    a -= a_offset;
#line 132 "SB01BY.f"
    f_dim1 = *m;
#line 132 "SB01BY.f"
    f_offset = 1 + f_dim1;
#line 132 "SB01BY.f"
    f -= f_offset;
#line 132 "SB01BY.f"
    --dwork;
#line 132 "SB01BY.f"

#line 132 "SB01BY.f"
    /* Function Body */
#line 132 "SB01BY.f"
    *info = 0;
#line 133 "SB01BY.f"
    if (*n == 1) {

/*        The case N = 1. */

#line 137 "SB01BY.f"
	if (*m > 1) {
#line 137 "SB01BY.f"
	    dlarfg_(m, &b[b_dim1 + 1], &b[(b_dim1 << 1) + 1], n, &tau1);
#line 137 "SB01BY.f"
	}
#line 139 "SB01BY.f"
	b1 = b[b_dim1 + 1];
#line 140 "SB01BY.f"
	if (abs(b1) <= *tol) {

/*           The pair (A,B) is uncontrollable. */

#line 144 "SB01BY.f"
	    *info = 1;
#line 145 "SB01BY.f"
	    return 0;
#line 146 "SB01BY.f"
	}

#line 148 "SB01BY.f"
	f[f_dim1 + 1] = (*s - a[a_dim1 + 1]) / b1;
#line 149 "SB01BY.f"
	if (*m > 1) {
#line 150 "SB01BY.f"
	    i__1 = *m - 1;
#line 150 "SB01BY.f"
	    dlaset_("Full", &i__1, &c__1, &c_b4, &c_b4, &f[f_dim1 + 2], m, (
		    ftnlen)4);
#line 151 "SB01BY.f"
	    dlatzm_("Left", m, n, &b[(b_dim1 << 1) + 1], n, &tau1, &f[f_dim1 
		    + 1], &f[f_dim1 + 2], m, &dwork[1], (ftnlen)4);
#line 153 "SB01BY.f"
	}
#line 154 "SB01BY.f"
	return 0;
#line 155 "SB01BY.f"
    }

/*     In the sequel N = 2. */

/*     Compute the singular value decomposition of B in the form */

/*                    ( V  0 )                ( B1 0  ) */
/*     B = U*( G1 0 )*(      )*H2*H1 ,   G1 = (       ), */
/*                    ( 0  I )                ( 0  B2 ) */

/*               ( CU   SU )          ( CV   SV ) */
/*     where U = (         )  and V = (         )  are orthogonal */
/*               (-SU   CU )          (-SV   CV ) */

/*     rotations and H1 and H2 are elementary Householder reflectors. */
/*     ABS(B1) and ABS(B2) are the singular values of matrix B, */
/*     with ABS(B1) >= ABS(B2). */

/*     Reduce first B to the lower bidiagonal form  ( B1  0  ... 0 ). */
/*                                                  ( B21 B2 ... 0 ) */
#line 175 "SB01BY.f"
    if (*m == 1) {

/*        Initialization for the case M = 1; no reduction required. */

#line 179 "SB01BY.f"
	b1 = b[b_dim1 + 1];
#line 180 "SB01BY.f"
	b21 = b[b_dim1 + 2];
#line 181 "SB01BY.f"
	b2 = 0.;
#line 182 "SB01BY.f"
    } else {

/*        Postmultiply B with elementary Householder reflectors H1 */
/*        and H2. */

#line 187 "SB01BY.f"
	dlarfg_(m, &b[b_dim1 + 1], &b[(b_dim1 << 1) + 1], n, &tau1);
#line 188 "SB01BY.f"
	i__1 = *n - 1;
#line 188 "SB01BY.f"
	dlatzm_("Right", &i__1, m, &b[(b_dim1 << 1) + 1], n, &tau1, &b[b_dim1 
		+ 2], &b[(b_dim1 << 1) + 2], n, &dwork[1], (ftnlen)5);
#line 190 "SB01BY.f"
	b1 = b[b_dim1 + 1];
#line 191 "SB01BY.f"
	b21 = b[b_dim1 + 2];
#line 192 "SB01BY.f"
	if (*m > 2) {
#line 192 "SB01BY.f"
	    i__1 = *m - 1;
#line 192 "SB01BY.f"
	    dlarfg_(&i__1, &b[(b_dim1 << 1) + 2], &b[b_dim1 * 3 + 2], n, &
		    tau2);
#line 192 "SB01BY.f"
	}
#line 194 "SB01BY.f"
	b2 = b[(b_dim1 << 1) + 2];
#line 195 "SB01BY.f"
    }

/*     Reduce B to a diagonal form by premultiplying and postmultiplying */
/*     it with orthogonal rotations U and V, respectively, and order the */
/*     diagonal elements to have decreasing magnitudes. */
/*     Note: B2 has been set to zero if M = 1. Thus in the following */
/*     computations the case M = 1 need not to be distinguished. */
/*     Note also that LAPACK routine DLASV2 assumes an upper triangular */
/*     matrix, so the results should be adapted. */

#line 205 "SB01BY.f"
    dlasv2_(&b1, &b21, &b2, &x, &y, &su, &cu, &sv, &cv);
#line 206 "SB01BY.f"
    su = -su;
#line 207 "SB01BY.f"
    b1 = y;
#line 208 "SB01BY.f"
    b2 = x;

/*     Compute  A1 = U'*A*U. */

#line 212 "SB01BY.f"
    drot_(&c__2, &a[a_dim1 + 2], &c__2, &a[a_dim1 + 1], &c__2, &cu, &su);
#line 213 "SB01BY.f"
    drot_(&c__2, &a[(a_dim1 << 1) + 1], &c__1, &a[a_dim1 + 1], &c__1, &cu, &
	    su);

/*     Compute the rank of B and check the controllability of the */
/*     pair (A,B). */

#line 218 "SB01BY.f"
    ir = 0;
#line 219 "SB01BY.f"
    if (abs(b2) > *tol) {
#line 219 "SB01BY.f"
	++ir;
#line 219 "SB01BY.f"
    }
#line 220 "SB01BY.f"
    if (abs(b1) > *tol) {
#line 220 "SB01BY.f"
	++ir;
#line 220 "SB01BY.f"
    }
#line 221 "SB01BY.f"
    if (ir == 0 || ir == 1 && (d__1 = a[a_dim1 + 2], abs(d__1)) <= *tol) {
#line 222 "SB01BY.f"
	f[f_dim1 + 1] = cu;
#line 223 "SB01BY.f"
	f[(f_dim1 << 1) + 1] = -su;

/*        The pair (A,B) is uncontrollable. */

#line 227 "SB01BY.f"
	*info = 1;
#line 228 "SB01BY.f"
	return 0;
#line 229 "SB01BY.f"
    }

/*     Compute F1 which assigns N poles for the reduced pair (A1,G1). */

#line 233 "SB01BY.f"
    x = dlamc3_(&b1, &b2);
#line 234 "SB01BY.f"
    if (x == b1) {

/*        Rank one G1. */

#line 238 "SB01BY.f"
	f[f_dim1 + 1] = (*s - (a[a_dim1 + 1] + a[(a_dim1 << 1) + 2])) / b1;
#line 239 "SB01BY.f"
	f[(f_dim1 << 1) + 1] = -(a[(a_dim1 << 1) + 2] * (a[(a_dim1 << 1) + 2] 
		- *s) + a[a_dim1 + 2] * a[(a_dim1 << 1) + 1] + *p) / a[a_dim1 
		+ 2] / b1;
#line 241 "SB01BY.f"
	if (*m > 1) {
#line 242 "SB01BY.f"
	    f[f_dim1 + 2] = 0.;
#line 243 "SB01BY.f"
	    f[(f_dim1 << 1) + 2] = 0.;
#line 244 "SB01BY.f"
	}
#line 245 "SB01BY.f"
    } else {

/*        Rank two G1. */

#line 249 "SB01BY.f"
	z__ = (*s - (a[a_dim1 + 1] + a[(a_dim1 << 1) + 2])) / (b1 * b1 + b2 * 
		b2);
#line 250 "SB01BY.f"
	f[f_dim1 + 1] = b1 * z__;
#line 251 "SB01BY.f"
	f[(f_dim1 << 1) + 2] = b2 * z__;

/*        Compute an approximation for the minimum norm parameter */
/*        selection. */

#line 256 "SB01BY.f"
	x = a[a_dim1 + 1] + b1 * f[f_dim1 + 1];
#line 257 "SB01BY.f"
	c__ = x * (*s - x) - *p;
#line 258 "SB01BY.f"
	if (c__ >= 0.) {
#line 259 "SB01BY.f"
	    sig = 1.;
#line 260 "SB01BY.f"
	} else {
#line 261 "SB01BY.f"
	    sig = -1.;
#line 262 "SB01BY.f"
	}
#line 263 "SB01BY.f"
	s12 = b1 / b2;
#line 264 "SB01BY.f"
	s21 = b2 / b1;
#line 265 "SB01BY.f"
	c11 = 0.;
#line 266 "SB01BY.f"
	c12 = 1.;
#line 267 "SB01BY.f"
	c21 = sig * s12 * c__;
#line 268 "SB01BY.f"
	c22 = a[(a_dim1 << 1) + 1] - sig * s12 * a[a_dim1 + 2];
#line 269 "SB01BY.f"
	dlanv2_(&c11, &c12, &c21, &c22, &wr, &wi, &wr1, &wi1, &cs, &sn);
#line 270 "SB01BY.f"
	if ((d__1 = wr - a[(a_dim1 << 1) + 1], abs(d__1)) > (d__2 = wr1 - a[(
		a_dim1 << 1) + 1], abs(d__2))) {
#line 271 "SB01BY.f"
	    r__ = wr1;
#line 272 "SB01BY.f"
	} else {
#line 273 "SB01BY.f"
	    r__ = wr;
#line 274 "SB01BY.f"
	}

/*        Perform Newton iteration to solve the equation for minimum. */

#line 278 "SB01BY.f"
	c0 = -c__ * c__;
#line 279 "SB01BY.f"
	c1 = c__ * a[a_dim1 + 2];
#line 280 "SB01BY.f"
	c4 = s21 * s21;
#line 281 "SB01BY.f"
	c3 = -c4 * a[(a_dim1 << 1) + 1];
#line 282 "SB01BY.f"
	dc0 = c1;
#line 283 "SB01BY.f"
	dc2 = c3 * 3.;
#line 284 "SB01BY.f"
	dc3 = c4 * 4.;

#line 286 "SB01BY.f"
	for (j = 1; j <= 10; ++j) {
#line 287 "SB01BY.f"
	    x = c0 + r__ * (c1 + r__ * r__ * (c3 + r__ * c4));
#line 288 "SB01BY.f"
	    y = dc0 + r__ * r__ * (dc2 + r__ * dc3);
#line 289 "SB01BY.f"
	    if (y == 0.) {
#line 289 "SB01BY.f"
		goto L20;
#line 289 "SB01BY.f"
	    }
#line 290 "SB01BY.f"
	    rn = r__ - x / y;
#line 291 "SB01BY.f"
	    absr = abs(r__);
#line 292 "SB01BY.f"
	    diffr = (d__1 = r__ - rn, abs(d__1));
#line 293 "SB01BY.f"
	    z__ = dlamc3_(&absr, &diffr);
#line 294 "SB01BY.f"
	    if (z__ == absr) {
#line 294 "SB01BY.f"
		goto L20;
#line 294 "SB01BY.f"
	    }
#line 296 "SB01BY.f"
	    r__ = rn;
#line 297 "SB01BY.f"
/* L10: */
#line 297 "SB01BY.f"
	}

#line 299 "SB01BY.f"
L20:
#line 300 "SB01BY.f"
	if (r__ == 0.) {
#line 300 "SB01BY.f"
	    r__ = dlamch_("Epsilon", (ftnlen)7);
#line 300 "SB01BY.f"
	}
#line 301 "SB01BY.f"
	f[(f_dim1 << 1) + 1] = (r__ - a[(a_dim1 << 1) + 1]) / b1;
#line 302 "SB01BY.f"
	f[f_dim1 + 2] = (c__ / r__ - a[a_dim1 + 2]) / b2;
#line 303 "SB01BY.f"
    }

/*     Back-transform F1. Compute first F1*U'. */

#line 307 "SB01BY.f"
    i__1 = min(*m,2);
#line 307 "SB01BY.f"
    drot_(&i__1, &f[f_dim1 + 1], &c__1, &f[(f_dim1 << 1) + 1], &c__1, &cu, &
	    su);
#line 308 "SB01BY.f"
    if (*m == 1) {
#line 308 "SB01BY.f"
	return 0;
#line 308 "SB01BY.f"
    }

/*     Compute V'*F1. */

#line 313 "SB01BY.f"
    drot_(&c__2, &f[f_dim1 + 2], m, &f[f_dim1 + 1], m, &cv, &sv);

/*               ( F1 ) */
/*     Form  F = (    ) . */
/*               ( 0  ) */

#line 319 "SB01BY.f"
    if (*m > *n) {
#line 319 "SB01BY.f"
	i__1 = *m - *n;
#line 319 "SB01BY.f"
	dlaset_("Full", &i__1, n, &c_b4, &c_b4, &f[*n + 1 + f_dim1], m, (
		ftnlen)4);
#line 319 "SB01BY.f"
    }

/*     Compute H1*H2*F. */

#line 324 "SB01BY.f"
    if (*m > 2) {
#line 324 "SB01BY.f"
	i__1 = *m - 1;
#line 324 "SB01BY.f"
	dlatzm_("Left", &i__1, n, &b[b_dim1 * 3 + 2], n, &tau2, &f[f_dim1 + 2]
		, &f[f_dim1 + 3], m, &dwork[1], (ftnlen)4);
#line 324 "SB01BY.f"
    }
#line 327 "SB01BY.f"
    dlatzm_("Left", m, n, &b[(b_dim1 << 1) + 1], n, &tau1, &f[f_dim1 + 1], &f[
	    f_dim1 + 2], m, &dwork[1], (ftnlen)4);

#line 330 "SB01BY.f"
    return 0;
/* *** Last line of SB01BY *** */
} /* sb01by_ */

