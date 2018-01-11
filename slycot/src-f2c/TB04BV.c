#line 1 "TB04BV.f"
/* TB04BV.f -- translated by f2c (version 20100827).
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

#line 1 "TB04BV.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int tb04bv_(char *order, integer *p, integer *m, integer *md,
	 integer *ign, integer *ldign, integer *igd, integer *ldigd, 
	doublereal *gn, doublereal *gd, doublereal *d__, integer *ldd, 
	doublereal *tol, integer *info, ftnlen order_len)
{
    /* System generated locals */
    integer d_dim1, d_offset, igd_dim1, igd_offset, ign_dim1, ign_offset, 
	    i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, ii, nd, kk, km, nn;
    static doublereal dij, eps;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen);
    static logical ascend;
    extern integer idamax_(integer *, doublereal *, integer *);
    static doublereal toldef;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);


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

/*     To separate the strictly proper part G0 from the constant part D */
/*     of an P-by-M proper transfer function matrix G. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ORDER   CHARACTER*1 */
/*             Specifies the order in which the polynomial coefficients */
/*             of the transfer function matrix are stored, as follows: */
/*             = 'I':  Increasing order of powers of the indeterminate; */
/*             = 'D':  Decreasing order of powers of the indeterminate. */

/*     Input/Output Parameters */

/*     P       (input) INTEGER */
/*             The number of the system outputs.  P >= 0. */

/*     M       (input) INTEGER */
/*             The number of the system inputs.  M >= 0. */

/*     MD      (input) INTEGER */
/*             The maximum degree of the polynomials in G, plus 1, i.e., */
/*             MD = MAX(IGD(I,J)) + 1. */
/*                  I,J */

/*     IGN     (input/output) INTEGER array, dimension (LDIGN,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the degrees of the numerator polynomials in G: */
/*             the (i,j) element of IGN must contain the degree of the */
/*             numerator polynomial of the polynomial ratio G(i,j). */
/*             On exit, the leading P-by-M part of this array contains */
/*             the degrees of the numerator polynomials in G0. */

/*     LDIGN   INTEGER */
/*             The leading dimension of array IGN.  LDIGN >= max(1,P). */

/*     IGD     (input) INTEGER array, dimension (LDIGD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             degrees of the denominator polynomials in G (and G0): */
/*             the (i,j) element of IGD contains the degree of the */
/*             denominator polynomial of the polynomial ratio G(i,j). */

/*     LDIGD   INTEGER */
/*             The leading dimension of array IGD.  LDIGD >= max(1,P). */

/*     GN      (input/output) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             On entry, this array must contain the coefficients of the */
/*             numerator polynomials, Num(i,j), of the transfer function */
/*             matrix G. The polynomials are stored in a column-wise */
/*             order, i.e., Num(1,1), Num(2,1), ..., Num(P,1), Num(1,2), */
/*             Num(2,2), ..., Num(P,2), ..., Num(1,M), Num(2,M), ..., */
/*             Num(P,M); MD memory locations are reserved for each */
/*             polynomial, hence, the (i,j) polynomial is stored starting */
/*             from the location ((j-1)*P+i-1)*MD+1. The coefficients */
/*             appear in increasing or decreasing order of the powers */
/*             of the indeterminate, according to ORDER. */
/*             On exit, this array contains the coefficients of the */
/*             numerator polynomials of the strictly proper part G0 of */
/*             the transfer function matrix G, stored similarly. */

/*     GD      (input) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             This array must contain the coefficients of the */
/*             denominator polynomials, Den(i,j), of the transfer */
/*             function matrix G. The polynomials are stored as for the */
/*             numerator polynomials. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array contains the */
/*             matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= max(1,P). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in determining the degrees of */
/*             the numerators Num0(i,j) of the strictly proper part of */
/*             the transfer function matrix G. If the user sets TOL > 0, */
/*             then the given value of TOL is used as an absolute */
/*             tolerance; the leading coefficients with absolute value */
/*             less than TOL are considered neglijible. If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by TOLDEF = IGN(i,j)*EPS*NORM( Num(i,j) ) is used */
/*             instead, where EPS is the machine precision (see LAPACK */
/*             Library routine DLAMCH), and NORM denotes the infinity */
/*             norm (the maximum coefficient in absolute value). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the transfer function matrix is not proper; */
/*             = 2:  if a denominator polynomial is null. */

/*     METHOD */

/*     The (i,j) entry of the real matrix D is zero, if the degree of */
/*     Num(i,j), IGN(i,j), is less than the degree of Den(i,j), IGD(i,j), */
/*     and it is given by the ratio of the leading coefficients of */
/*     Num(i,j) and Den(i,j), if IGN(i,j) is equal to IGD(i,j), */
/*     for i = 1 : P, and for j = 1 : M. */

/*     FURTHER COMMENTS */

/*     For maximum efficiency of index calculations, GN and GD are */
/*     implemented as one-dimensional arrays. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2002. */
/*     Based on the BIMASC Library routine TMPRP by A. Varga. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Feb. 2004. */

/*     KEYWORDS */

/*     State-space representation, transfer function. */

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

/*     Test the input scalar parameters. */

#line 178 "TB04BV.f"
    /* Parameter adjustments */
#line 178 "TB04BV.f"
    ign_dim1 = *ldign;
#line 178 "TB04BV.f"
    ign_offset = 1 + ign_dim1;
#line 178 "TB04BV.f"
    ign -= ign_offset;
#line 178 "TB04BV.f"
    igd_dim1 = *ldigd;
#line 178 "TB04BV.f"
    igd_offset = 1 + igd_dim1;
#line 178 "TB04BV.f"
    igd -= igd_offset;
#line 178 "TB04BV.f"
    --gn;
#line 178 "TB04BV.f"
    --gd;
#line 178 "TB04BV.f"
    d_dim1 = *ldd;
#line 178 "TB04BV.f"
    d_offset = 1 + d_dim1;
#line 178 "TB04BV.f"
    d__ -= d_offset;
#line 178 "TB04BV.f"

#line 178 "TB04BV.f"
    /* Function Body */
#line 178 "TB04BV.f"
    *info = 0;
#line 179 "TB04BV.f"
    ascend = lsame_(order, "I", (ftnlen)1, (ftnlen)1);
#line 180 "TB04BV.f"
    if (! ascend && ! lsame_(order, "D", (ftnlen)1, (ftnlen)1)) {
#line 181 "TB04BV.f"
	*info = -1;
#line 182 "TB04BV.f"
    } else if (*p < 0) {
#line 183 "TB04BV.f"
	*info = -2;
#line 184 "TB04BV.f"
    } else if (*m < 0) {
#line 185 "TB04BV.f"
	*info = -3;
#line 186 "TB04BV.f"
    } else if (*md < 1) {
#line 187 "TB04BV.f"
	*info = -4;
#line 188 "TB04BV.f"
    } else if (*ldign < max(1,*p)) {
#line 189 "TB04BV.f"
	*info = -6;
#line 190 "TB04BV.f"
    } else if (*ldigd < max(1,*p)) {
#line 191 "TB04BV.f"
	*info = -8;
#line 192 "TB04BV.f"
    } else if (*ldd < max(1,*p)) {
#line 193 "TB04BV.f"
	*info = -12;
#line 194 "TB04BV.f"
    }

#line 196 "TB04BV.f"
    if (*info != 0) {

/*        Error return. */

#line 200 "TB04BV.f"
	i__1 = -(*info);
#line 200 "TB04BV.f"
	xerbla_("TB04BV", &i__1, (ftnlen)6);
#line 201 "TB04BV.f"
	return 0;
#line 202 "TB04BV.f"
    }

/*     Quick return if possible. */

#line 206 "TB04BV.f"
    if (min(*p,*m) == 0) {
#line 206 "TB04BV.f"
	return 0;
#line 206 "TB04BV.f"
    }

/*     Prepare the computation of the default tolerance. */

#line 211 "TB04BV.f"
    toldef = *tol;
#line 212 "TB04BV.f"
    if (toldef <= 0.) {
#line 212 "TB04BV.f"
	eps = dlamch_("Epsilon", (ftnlen)7);
#line 212 "TB04BV.f"
    }

#line 215 "TB04BV.f"
    k = 1;

#line 217 "TB04BV.f"
    if (ascend) {

/*        Polynomial coefficients are stored in increasing order. */

#line 221 "TB04BV.f"
	i__1 = *m;
#line 221 "TB04BV.f"
	for (j = 1; j <= i__1; ++j) {

#line 223 "TB04BV.f"
	    i__2 = *p;
#line 223 "TB04BV.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 224 "TB04BV.f"
		nn = ign[i__ + j * ign_dim1];
#line 225 "TB04BV.f"
		nd = igd[i__ + j * igd_dim1];
#line 226 "TB04BV.f"
		if (nn > nd) {

/*                 Error return: the transfer function matrix is */
/*                               not proper. */

#line 231 "TB04BV.f"
		    *info = 1;
#line 232 "TB04BV.f"
		    return 0;
#line 233 "TB04BV.f"
		} else if (nn < nd || nd == 0 && gn[k] == 0.) {
#line 235 "TB04BV.f"
		    d__[i__ + j * d_dim1] = 0.;
#line 236 "TB04BV.f"
		} else {

/*                 Here NN = ND. */

#line 240 "TB04BV.f"
		    kk = k + nn;

#line 242 "TB04BV.f"
		    if (gd[kk] == 0.) {

/*                    Error return: the denominator is null. */

#line 246 "TB04BV.f"
			*info = 2;
#line 247 "TB04BV.f"
			return 0;
#line 248 "TB04BV.f"
		    }

#line 250 "TB04BV.f"
		    dij = gn[kk] / gd[kk];
#line 251 "TB04BV.f"
		    d__[i__ + j * d_dim1] = dij;
#line 252 "TB04BV.f"
		    gn[kk] = 0.;
#line 253 "TB04BV.f"
		    if (nn > 0) {
#line 254 "TB04BV.f"
			d__1 = -dij;
#line 254 "TB04BV.f"
			daxpy_(&nn, &d__1, &gd[k], &c__1, &gn[k], &c__1);
#line 255 "TB04BV.f"
			if (*tol <= 0.) {
#line 255 "TB04BV.f"
			    toldef = (doublereal) nn * eps * (d__1 = gn[
				    idamax_(&nn, &gn[k], &c__1)], abs(d__1));
#line 255 "TB04BV.f"
			}
#line 258 "TB04BV.f"
			km = nn;
#line 259 "TB04BV.f"
			i__3 = km;
#line 259 "TB04BV.f"
			for (ii = 1; ii <= i__3; ++ii) {
#line 260 "TB04BV.f"
			    --kk;
#line 261 "TB04BV.f"
			    --nn;
#line 262 "TB04BV.f"
			    if ((d__1 = gn[kk], abs(d__1)) > toldef) {
#line 262 "TB04BV.f"
				goto L20;
#line 262 "TB04BV.f"
			    }
#line 264 "TB04BV.f"
/* L10: */
#line 264 "TB04BV.f"
			}

#line 266 "TB04BV.f"
L20:

#line 268 "TB04BV.f"
			ign[i__ + j * ign_dim1] = nn;
#line 269 "TB04BV.f"
		    }
#line 270 "TB04BV.f"
		}
#line 271 "TB04BV.f"
		k += *md;
#line 272 "TB04BV.f"
/* L30: */
#line 272 "TB04BV.f"
	    }

#line 274 "TB04BV.f"
/* L40: */
#line 274 "TB04BV.f"
	}

#line 276 "TB04BV.f"
    } else {

/*        Polynomial coefficients are stored in decreasing order. */

#line 280 "TB04BV.f"
	i__1 = *m;
#line 280 "TB04BV.f"
	for (j = 1; j <= i__1; ++j) {

#line 282 "TB04BV.f"
	    i__2 = *p;
#line 282 "TB04BV.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 283 "TB04BV.f"
		nn = ign[i__ + j * ign_dim1];
#line 284 "TB04BV.f"
		nd = igd[i__ + j * igd_dim1];
#line 285 "TB04BV.f"
		if (nn > nd) {

/*                 Error return: the transfer function matrix is */
/*                               not proper. */

#line 290 "TB04BV.f"
		    *info = 1;
#line 291 "TB04BV.f"
		    return 0;
#line 292 "TB04BV.f"
		} else if (nn < nd || nd == 0 && gn[k] == 0.) {
#line 294 "TB04BV.f"
		    d__[i__ + j * d_dim1] = 0.;
#line 295 "TB04BV.f"
		} else {

/*                 Here NN = ND. */

#line 299 "TB04BV.f"
		    kk = k;

#line 301 "TB04BV.f"
		    if (gd[kk] == 0.) {

/*                    Error return: the denominator is null. */

#line 305 "TB04BV.f"
			*info = 2;
#line 306 "TB04BV.f"
			return 0;
#line 307 "TB04BV.f"
		    }

#line 309 "TB04BV.f"
		    dij = gn[kk] / gd[kk];
#line 310 "TB04BV.f"
		    d__[i__ + j * d_dim1] = dij;
#line 311 "TB04BV.f"
		    gn[kk] = 0.;
#line 312 "TB04BV.f"
		    if (nn > 0) {
#line 313 "TB04BV.f"
			d__1 = -dij;
#line 313 "TB04BV.f"
			daxpy_(&nn, &d__1, &gd[k + 1], &c__1, &gn[k + 1], &
				c__1);
#line 314 "TB04BV.f"
			if (*tol <= 0.) {
#line 314 "TB04BV.f"
			    toldef = (doublereal) nn * eps * (d__1 = gn[
				    idamax_(&nn, &gn[k + 1], &c__1)], abs(
				    d__1));
#line 314 "TB04BV.f"
			}
#line 317 "TB04BV.f"
			km = nn;
#line 318 "TB04BV.f"
			i__3 = km;
#line 318 "TB04BV.f"
			for (ii = 1; ii <= i__3; ++ii) {
#line 319 "TB04BV.f"
			    ++kk;
#line 320 "TB04BV.f"
			    --nn;
#line 321 "TB04BV.f"
			    if ((d__1 = gn[kk], abs(d__1)) > toldef) {
#line 321 "TB04BV.f"
				goto L60;
#line 321 "TB04BV.f"
			    }
#line 323 "TB04BV.f"
/* L50: */
#line 323 "TB04BV.f"
			}

#line 325 "TB04BV.f"
L60:

#line 327 "TB04BV.f"
			ign[i__ + j * ign_dim1] = nn;
#line 328 "TB04BV.f"
			i__3 = nn;
#line 328 "TB04BV.f"
			for (ii = 0; ii <= i__3; ++ii) {
#line 329 "TB04BV.f"
			    gn[k + ii] = gn[kk + ii];
#line 330 "TB04BV.f"
/* L70: */
#line 330 "TB04BV.f"
			}

#line 332 "TB04BV.f"
		    }
#line 333 "TB04BV.f"
		}
#line 334 "TB04BV.f"
		k += *md;
#line 335 "TB04BV.f"
/* L80: */
#line 335 "TB04BV.f"
	    }

#line 337 "TB04BV.f"
/* L90: */
#line 337 "TB04BV.f"
	}

#line 339 "TB04BV.f"
    }

#line 341 "TB04BV.f"
    return 0;
/* *** Last line of TB04BV *** */
} /* tb04bv_ */

