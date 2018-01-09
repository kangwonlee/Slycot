#line 1 "MB02UV.f"
/* MB02UV.f -- translated by f2c (version 20100827).
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

#line 1 "MB02UV.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b9 = -1.;

/* Subroutine */ int mb02uv_(integer *n, doublereal *a, integer *lda, integer 
	*ipiv, integer *jpiv, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, ip, jp;
    static doublereal eps;
    static integer ipv, jpv;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal smin, xmax;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dswap_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal bignum, smlnum;


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

/*     To compute an LU factorization, using complete pivoting, of the */
/*     N-by-N matrix A. The factorization has the form A = P * L * U * Q, */
/*     where P and Q are permutation matrices, L is lower triangular with */
/*     unit diagonal elements and U is upper triangular. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA, N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A to be factored. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the factors L and U from the factorization A = P*L*U*Q; */
/*             the unit diagonal elements of L are not stored. If U(k, k) */
/*             appears to be less than SMIN, U(k, k) is given the value */
/*             of SMIN, giving a nonsingular perturbed system. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1, N). */

/*     IPIV    (output) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= i <= N, row i of the */
/*             matrix has been interchanged with row IPIV(i). */

/*     JPIV    (output) INTEGER array, dimension (N) */
/*             The pivot indices; for 1 <= j <= N, column j of the */
/*             matrix has been interchanged with column JPIV(j). */

/*     Error indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = k:  U(k, k) is likely to produce owerflow if one tries */
/*                   to solve for x in Ax = b. So U is perturbed to get */
/*                   a nonsingular system. This is a warning. */

/*     FURTHER COMMENTS */

/*     In the interests of speed, this routine does not check the input */
/*     for errors. It should only be used to factorize matrices A of */
/*     very small order. */

/*     CONTRIBUTOR */

/*     Bo Kagstrom and Peter Poromaa, Univ. of Umea, Sweden, Nov. 1993. */

/*     REVISIONS */

/*     April 1998 (T. Penzl). */
/*     Sep. 1998 (V. Sima). */
/*     March 1999 (V. Sima). */
/*     March 2004 (V. Sima). */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Set constants to control owerflow. */
#line 104 "MB02UV.f"
    /* Parameter adjustments */
#line 104 "MB02UV.f"
    a_dim1 = *lda;
#line 104 "MB02UV.f"
    a_offset = 1 + a_dim1;
#line 104 "MB02UV.f"
    a -= a_offset;
#line 104 "MB02UV.f"
    --ipiv;
#line 104 "MB02UV.f"
    --jpiv;
#line 104 "MB02UV.f"

#line 104 "MB02UV.f"
    /* Function Body */
#line 104 "MB02UV.f"
    *info = 0;
#line 105 "MB02UV.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 106 "MB02UV.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / eps;
#line 107 "MB02UV.f"
    bignum = 1. / smlnum;
#line 108 "MB02UV.f"
    dlabad_(&smlnum, &bignum);

/*     Find max element in matrix A. */

#line 112 "MB02UV.f"
    ipv = 1;
#line 113 "MB02UV.f"
    jpv = 1;
#line 114 "MB02UV.f"
    xmax = 0.;
#line 115 "MB02UV.f"
    i__1 = *n;
#line 115 "MB02UV.f"
    for (jp = 1; jp <= i__1; ++jp) {
#line 116 "MB02UV.f"
	i__2 = *n;
#line 116 "MB02UV.f"
	for (ip = 1; ip <= i__2; ++ip) {
#line 117 "MB02UV.f"
	    if ((d__1 = a[ip + jp * a_dim1], abs(d__1)) > xmax) {
#line 118 "MB02UV.f"
		xmax = (d__1 = a[ip + jp * a_dim1], abs(d__1));
#line 119 "MB02UV.f"
		ipv = ip;
#line 120 "MB02UV.f"
		jpv = jp;
#line 121 "MB02UV.f"
	    }
#line 122 "MB02UV.f"
/* L20: */
#line 122 "MB02UV.f"
	}
#line 123 "MB02UV.f"
/* L40: */
#line 123 "MB02UV.f"
    }
/* Computing MAX */
#line 124 "MB02UV.f"
    d__1 = eps * xmax;
#line 124 "MB02UV.f"
    smin = max(d__1,smlnum);

/*     Swap rows. */

#line 128 "MB02UV.f"
    if (ipv != 1) {
#line 128 "MB02UV.f"
	dswap_(n, &a[ipv + a_dim1], lda, &a[a_dim1 + 1], lda);
#line 128 "MB02UV.f"
    }
#line 129 "MB02UV.f"
    ipiv[1] = ipv;

/*     Swap columns. */

#line 133 "MB02UV.f"
    if (jpv != 1) {
#line 133 "MB02UV.f"
	dswap_(n, &a[jpv * a_dim1 + 1], &c__1, &a[a_dim1 + 1], &c__1);
#line 133 "MB02UV.f"
    }
#line 134 "MB02UV.f"
    jpiv[1] = jpv;

/*     Check for singularity. */

#line 138 "MB02UV.f"
    if ((d__1 = a[a_dim1 + 1], abs(d__1)) < smin) {
#line 139 "MB02UV.f"
	*info = 1;
#line 140 "MB02UV.f"
	a[a_dim1 + 1] = smin;
#line 141 "MB02UV.f"
    }
#line 142 "MB02UV.f"
    if (*n > 1) {
#line 143 "MB02UV.f"
	i__1 = *n - 1;
#line 143 "MB02UV.f"
	d__1 = 1. / a[a_dim1 + 1];
#line 143 "MB02UV.f"
	dscal_(&i__1, &d__1, &a[a_dim1 + 2], &c__1);
#line 144 "MB02UV.f"
	i__1 = *n - 1;
#line 144 "MB02UV.f"
	i__2 = *n - 1;
#line 144 "MB02UV.f"
	dger_(&i__1, &i__2, &c_b9, &a[a_dim1 + 2], &c__1, &a[(a_dim1 << 1) + 
		1], lda, &a[(a_dim1 << 1) + 2], lda);
#line 146 "MB02UV.f"
    }

/*     Factorize the rest of A with complete pivoting. */
/*     Set pivots less than SMIN to SMIN. */

#line 151 "MB02UV.f"
    i__1 = *n - 1;
#line 151 "MB02UV.f"
    for (i__ = 2; i__ <= i__1; ++i__) {

/*        Find max element in remaining matrix. */

#line 155 "MB02UV.f"
	ipv = i__;
#line 156 "MB02UV.f"
	jpv = i__;
#line 157 "MB02UV.f"
	xmax = 0.;
#line 158 "MB02UV.f"
	i__2 = *n;
#line 158 "MB02UV.f"
	for (jp = i__; jp <= i__2; ++jp) {
#line 159 "MB02UV.f"
	    i__3 = *n;
#line 159 "MB02UV.f"
	    for (ip = i__; ip <= i__3; ++ip) {
#line 160 "MB02UV.f"
		if ((d__1 = a[ip + jp * a_dim1], abs(d__1)) > xmax) {
#line 161 "MB02UV.f"
		    xmax = (d__1 = a[ip + jp * a_dim1], abs(d__1));
#line 162 "MB02UV.f"
		    ipv = ip;
#line 163 "MB02UV.f"
		    jpv = jp;
#line 164 "MB02UV.f"
		}
#line 165 "MB02UV.f"
/* L60: */
#line 165 "MB02UV.f"
	    }
#line 166 "MB02UV.f"
/* L80: */
#line 166 "MB02UV.f"
	}

/*        Swap rows. */

#line 170 "MB02UV.f"
	if (ipv != i__) {
#line 170 "MB02UV.f"
	    dswap_(n, &a[ipv + a_dim1], lda, &a[i__ + a_dim1], lda);
#line 170 "MB02UV.f"
	}
#line 171 "MB02UV.f"
	ipiv[i__] = ipv;

/*        Swap columns. */

#line 175 "MB02UV.f"
	if (jpv != i__) {
#line 175 "MB02UV.f"
	    dswap_(n, &a[jpv * a_dim1 + 1], &c__1, &a[i__ * a_dim1 + 1], &
		    c__1);
#line 175 "MB02UV.f"
	}
#line 176 "MB02UV.f"
	jpiv[i__] = jpv;

/*        Check for almost singularity. */

#line 180 "MB02UV.f"
	if ((d__1 = a[i__ + i__ * a_dim1], abs(d__1)) < smin) {
#line 181 "MB02UV.f"
	    *info = i__;
#line 182 "MB02UV.f"
	    a[i__ + i__ * a_dim1] = smin;
#line 183 "MB02UV.f"
	}
#line 184 "MB02UV.f"
	i__2 = *n - i__;
#line 184 "MB02UV.f"
	d__1 = 1. / a[i__ + i__ * a_dim1];
#line 184 "MB02UV.f"
	dscal_(&i__2, &d__1, &a[i__ + 1 + i__ * a_dim1], &c__1);
#line 185 "MB02UV.f"
	i__2 = *n - i__;
#line 185 "MB02UV.f"
	i__3 = *n - i__;
#line 185 "MB02UV.f"
	dger_(&i__2, &i__3, &c_b9, &a[i__ + 1 + i__ * a_dim1], &c__1, &a[i__ 
		+ (i__ + 1) * a_dim1], lda, &a[i__ + 1 + (i__ + 1) * a_dim1], 
		lda);
#line 187 "MB02UV.f"
/* L100: */
#line 187 "MB02UV.f"
    }
#line 188 "MB02UV.f"
    if ((d__1 = a[*n + *n * a_dim1], abs(d__1)) < smin) {
#line 189 "MB02UV.f"
	*info = *n;
#line 190 "MB02UV.f"
	a[*n + *n * a_dim1] = smin;
#line 191 "MB02UV.f"
    }

#line 193 "MB02UV.f"
    return 0;
/* *** Last line of MB02UV *** */
} /* mb02uv_ */

