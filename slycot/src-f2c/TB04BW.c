#line 1 "TB04BW.f"
/* TB04BW.f -- translated by f2c (version 20100827).
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

#line 1 "TB04BW.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int tb04bw_(char *order, integer *p, integer *m, integer *md,
	 integer *ign, integer *ldign, integer *igd, integer *ldigd, 
	doublereal *gn, doublereal *gd, doublereal *d__, integer *ldd, 
	integer *info, ftnlen order_len)
{
    /* System generated locals */
    integer d_dim1, d_offset, igd_dim1, igd_offset, ign_dim1, ign_offset, 
	    i__1, i__2, i__3;

    /* Local variables */
    static integer i__, j, k, ii, nd, kk, km, nn;
    static doublereal dij;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static logical ascend;
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

/*     To compute the sum of an P-by-M rational matrix G and a real */
/*     P-by-M matrix D. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ORDER   CHARACTER*1 */
/*             Specifies the order in which the polynomial coefficients */
/*             of the rational matrix are stored, as follows: */
/*             = 'I':  Increasing order of powers of the indeterminate; */
/*             = 'D':  Decreasing order of powers of the indeterminate. */

/*     Input/Output Parameters */

/*     P       (input) INTEGER */
/*             The number of the system outputs.  P >= 0. */

/*     M       (input) INTEGER */
/*             The number of the system inputs.  M >= 0. */

/*     MD      (input) INTEGER */
/*             The maximum degree of the polynomials in G, plus 1, i.e., */
/*             MD = MAX(IGN(I,J),IGD(I,J)) + 1. */
/*                  I,J */

/*     IGN     (input/output) INTEGER array, dimension (LDIGN,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the degrees of the numerator polynomials in G: */
/*             the (i,j) element of IGN must contain the degree of the */
/*             numerator polynomial of the polynomial ratio G(i,j). */
/*             On exit, the leading P-by-M part of this array contains */
/*             the degrees of the numerator polynomials in G + D. */

/*     LDIGN   INTEGER */
/*             The leading dimension of array IGN.  LDIGN >= max(1,P). */

/*     IGD     (input) INTEGER array, dimension (LDIGD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             degrees of the denominator polynomials in G (and G + D): */
/*             the (i,j) element of IGD contains the degree of the */
/*             denominator polynomial of the polynomial ratio G(i,j). */

/*     LDIGD   INTEGER */
/*             The leading dimension of array IGD.  LDIGD >= max(1,P). */

/*     GN      (input/output) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             On entry, this array must contain the coefficients of the */
/*             numerator polynomials, Num(i,j), of the rational matrix G. */
/*             The polynomials are stored in a column-wise order, i.e., */
/*             Num(1,1), Num(2,1), ..., Num(P,1), Num(1,2), Num(2,2), */
/*             ..., Num(P,2), ..., Num(1,M), Num(2,M), ..., Num(P,M); */
/*             MD memory locations are reserved for each polynomial, */
/*             hence, the (i,j) polynomial is stored starting from the */
/*             location ((j-1)*P+i-1)*MD+1. The coefficients appear in */
/*             increasing or decreasing order of the powers of the */
/*             indeterminate, according to ORDER. */
/*             On exit, this array contains the coefficients of the */
/*             numerator polynomials of the rational matrix G + D, */
/*             stored similarly. */

/*     GD      (input) DOUBLE PRECISION array, dimension (P*M*MD) */
/*             This array must contain the coefficients of the */
/*             denominator polynomials, Den(i,j), of the rational */
/*             matrix G. The polynomials are stored as for the */
/*             numerator polynomials. */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             matrix D. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= max(1,P). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The (i,j) entry of the real matrix D is added to the (i,j) entry */
/*     of the matrix G, g(i,j), which is a ratio of two polynomials, */
/*     for i = 1 : P, and for j = 1 : M. If g(i,j) = 0, it is assumed */
/*     that its denominator is 1. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable. */

/*     FURTHER COMMENTS */

/*     Often, the rational matrix G is found from a state-space */
/*     representation (A,B,C), and D corresponds to the direct */
/*     feedthrough matrix of the system. The sum G + D gives the */
/*     transfer function matrix of the system (A,B,C,D). */
/*     For maximum efficiency of index calculations, GN and GD are */
/*     implemented as one-dimensional arrays. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, May 2002. */
/*     Based on the BIMASC Library routine TMCADD by A. Varga. */

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

#line 165 "TB04BW.f"
    /* Parameter adjustments */
#line 165 "TB04BW.f"
    ign_dim1 = *ldign;
#line 165 "TB04BW.f"
    ign_offset = 1 + ign_dim1;
#line 165 "TB04BW.f"
    ign -= ign_offset;
#line 165 "TB04BW.f"
    igd_dim1 = *ldigd;
#line 165 "TB04BW.f"
    igd_offset = 1 + igd_dim1;
#line 165 "TB04BW.f"
    igd -= igd_offset;
#line 165 "TB04BW.f"
    --gn;
#line 165 "TB04BW.f"
    --gd;
#line 165 "TB04BW.f"
    d_dim1 = *ldd;
#line 165 "TB04BW.f"
    d_offset = 1 + d_dim1;
#line 165 "TB04BW.f"
    d__ -= d_offset;
#line 165 "TB04BW.f"

#line 165 "TB04BW.f"
    /* Function Body */
#line 165 "TB04BW.f"
    *info = 0;
#line 166 "TB04BW.f"
    ascend = lsame_(order, "I", (ftnlen)1, (ftnlen)1);
#line 167 "TB04BW.f"
    if (! ascend && ! lsame_(order, "D", (ftnlen)1, (ftnlen)1)) {
#line 168 "TB04BW.f"
	*info = -1;
#line 169 "TB04BW.f"
    } else if (*p < 0) {
#line 170 "TB04BW.f"
	*info = -2;
#line 171 "TB04BW.f"
    } else if (*m < 0) {
#line 172 "TB04BW.f"
	*info = -3;
#line 173 "TB04BW.f"
    } else if (*md < 1) {
#line 174 "TB04BW.f"
	*info = -4;
#line 175 "TB04BW.f"
    } else if (*ldign < max(1,*p)) {
#line 176 "TB04BW.f"
	*info = -6;
#line 177 "TB04BW.f"
    } else if (*ldigd < max(1,*p)) {
#line 178 "TB04BW.f"
	*info = -8;
#line 179 "TB04BW.f"
    } else if (*ldd < max(1,*p)) {
#line 180 "TB04BW.f"
	*info = -12;
#line 181 "TB04BW.f"
    }

#line 183 "TB04BW.f"
    if (*info != 0) {

/*        Error return. */

#line 187 "TB04BW.f"
	i__1 = -(*info);
#line 187 "TB04BW.f"
	xerbla_("TB04BW", &i__1, (ftnlen)6);
#line 188 "TB04BW.f"
	return 0;
#line 189 "TB04BW.f"
    }

/*     Quick return if possible. */

#line 193 "TB04BW.f"
    if (min(*p,*m) == 0) {
#line 193 "TB04BW.f"
	return 0;
#line 193 "TB04BW.f"
    }

#line 196 "TB04BW.f"
    k = 1;

#line 198 "TB04BW.f"
    if (ascend) {

/*        Polynomial coefficients are stored in increasing order. */

#line 202 "TB04BW.f"
	i__1 = *m;
#line 202 "TB04BW.f"
	for (j = 1; j <= i__1; ++j) {

#line 204 "TB04BW.f"
	    i__2 = *p;
#line 204 "TB04BW.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 205 "TB04BW.f"
		dij = d__[i__ + j * d_dim1];
#line 206 "TB04BW.f"
		if (dij != 0.) {
#line 207 "TB04BW.f"
		    nn = ign[i__ + j * ign_dim1];
#line 208 "TB04BW.f"
		    nd = igd[i__ + j * igd_dim1];
#line 209 "TB04BW.f"
		    if (nn == 0 && nd == 0) {
#line 210 "TB04BW.f"
			if (gn[k] == 0.) {
#line 211 "TB04BW.f"
			    gn[k] = dij;
#line 212 "TB04BW.f"
			} else {
#line 213 "TB04BW.f"
			    gn[k] += dij * gd[k];
#line 214 "TB04BW.f"
			}
#line 215 "TB04BW.f"
		    } else {
#line 216 "TB04BW.f"
			km = min(nn,nd) + 1;
#line 217 "TB04BW.f"
			daxpy_(&km, &dij, &gd[k], &c__1, &gn[k], &c__1);
#line 218 "TB04BW.f"
			if (nn < nd) {

#line 220 "TB04BW.f"
			    i__3 = k + nd;
#line 220 "TB04BW.f"
			    for (ii = k + km; ii <= i__3; ++ii) {
#line 221 "TB04BW.f"
				gn[ii] = dij * gd[ii];
#line 222 "TB04BW.f"
/* L10: */
#line 222 "TB04BW.f"
			    }

#line 224 "TB04BW.f"
			    ign[i__ + j * ign_dim1] = nd;
#line 225 "TB04BW.f"
			}
#line 226 "TB04BW.f"
		    }
#line 227 "TB04BW.f"
		}
#line 228 "TB04BW.f"
		k += *md;
#line 229 "TB04BW.f"
/* L20: */
#line 229 "TB04BW.f"
	    }

#line 231 "TB04BW.f"
/* L30: */
#line 231 "TB04BW.f"
	}

#line 233 "TB04BW.f"
    } else {

/*        Polynomial coefficients are stored in decreasing order. */

#line 237 "TB04BW.f"
	i__1 = *m;
#line 237 "TB04BW.f"
	for (j = 1; j <= i__1; ++j) {

#line 239 "TB04BW.f"
	    i__2 = *p;
#line 239 "TB04BW.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 240 "TB04BW.f"
		dij = d__[i__ + j * d_dim1];
#line 241 "TB04BW.f"
		if (dij != 0.) {
#line 242 "TB04BW.f"
		    nn = ign[i__ + j * ign_dim1];
#line 243 "TB04BW.f"
		    nd = igd[i__ + j * igd_dim1];
#line 244 "TB04BW.f"
		    if (nn == 0 && nd == 0) {
#line 245 "TB04BW.f"
			if (gn[k] == 0.) {
#line 246 "TB04BW.f"
			    gn[k] = dij;
#line 247 "TB04BW.f"
			} else {
#line 248 "TB04BW.f"
			    gn[k] += dij * gd[k];
#line 249 "TB04BW.f"
			}
#line 250 "TB04BW.f"
		    } else {
#line 251 "TB04BW.f"
			km = min(nn,nd) + 1;
#line 252 "TB04BW.f"
			if (nn < nd) {
#line 253 "TB04BW.f"
			    kk = k + nd - nn;

#line 255 "TB04BW.f"
			    i__3 = k;
#line 255 "TB04BW.f"
			    for (ii = k + nn; ii >= i__3; --ii) {
#line 256 "TB04BW.f"
				gn[ii + nd - nn] = gn[ii];
#line 257 "TB04BW.f"
/* L35: */
#line 257 "TB04BW.f"
			    }

#line 259 "TB04BW.f"
			    i__3 = kk - 1;
#line 259 "TB04BW.f"
			    for (ii = k; ii <= i__3; ++ii) {
#line 260 "TB04BW.f"
				gn[ii] = dij * gd[ii];
#line 261 "TB04BW.f"
/* L40: */
#line 261 "TB04BW.f"
			    }

#line 263 "TB04BW.f"
			    ign[i__ + j * ign_dim1] = nd;
#line 264 "TB04BW.f"
			    daxpy_(&km, &dij, &gd[kk], &c__1, &gn[kk], &c__1);
#line 265 "TB04BW.f"
			} else {
#line 266 "TB04BW.f"
			    kk = k + nn - nd;
#line 267 "TB04BW.f"
			    daxpy_(&km, &dij, &gd[k], &c__1, &gn[kk], &c__1);
#line 268 "TB04BW.f"
			}
#line 269 "TB04BW.f"
		    }
#line 270 "TB04BW.f"
		}
#line 271 "TB04BW.f"
		k += *md;
#line 272 "TB04BW.f"
/* L50: */
#line 272 "TB04BW.f"
	    }

#line 274 "TB04BW.f"
/* L60: */
#line 274 "TB04BW.f"
	}

#line 276 "TB04BW.f"
    }

#line 278 "TB04BW.f"
    return 0;
/* *** Last line of TB04BW *** */
} /* tb04bw_ */

