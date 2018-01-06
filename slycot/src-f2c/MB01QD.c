#line 1 "MB01QD.f"
/* MB01QD.f -- translated by f2c (version 20100827).
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

#line 1 "MB01QD.f"
/* Subroutine */ int mb01qd_(char *type__, integer *m, integer *n, integer *
	kl, integer *ku, doublereal *cfrom, doublereal *cto, integer *nbl, 
	integer *nrows, doublereal *a, integer *lda, integer *info, ftnlen 
	type_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, j, k, k1, k2, k3, k4;
    static doublereal mul, cto1;
    static logical done;
    static integer ifin, jfin;
    static doublereal ctoc;
    static integer jini;
    static logical noblc;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer itype;
    static doublereal cfrom1;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cfromc, bignum, smlnum;


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

/*     To multiply the M by N real matrix A by the real scalar CTO/CFROM. */
/*     This is done without over/underflow as long as the final result */
/*     CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that */
/*     A may be full, (block) upper triangular, (block) lower triangular, */
/*     (block) upper Hessenberg, or banded. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYPE    CHARACTER*1 */
/*             TYPE indices the storage type of the input matrix. */
/*             = 'G':  A is a full matrix. */
/*             = 'L':  A is a (block) lower triangular matrix. */
/*             = 'U':  A is a (block) upper triangular matrix. */
/*             = 'H':  A is a (block) upper Hessenberg matrix. */
/*             = 'B':  A is a symmetric band matrix with lower bandwidth */
/*                     KL and upper bandwidth KU and with the only the */
/*                     lower half stored. */
/*             = 'Q':  A is a symmetric band matrix with lower bandwidth */
/*                     KL and upper bandwidth KU and with the only the */
/*                     upper half stored. */
/*             = 'Z':  A is a band matrix with lower bandwidth KL and */
/*                     upper bandwidth KU. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of the matrix A.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of columns of the matrix A.  N >= 0. */

/*     KL      (input) INTEGER */
/*             The lower bandwidth of A.  Referenced only if TYPE = 'B', */
/*             'Q' or 'Z'. */

/*     KU      (input) INTEGER */
/*             The upper bandwidth of A.  Referenced only if TYPE = 'B', */
/*             'Q' or 'Z'. */

/*     CFROM   (input) DOUBLE PRECISION */
/*     CTO     (input) DOUBLE PRECISION */
/*             The matrix A is multiplied by CTO/CFROM. A(I,J) is */
/*             computed without over/underflow if the final result */
/*             CTO*A(I,J)/CFROM can be represented without over/ */
/*             underflow.  CFROM must be nonzero. */

/*     NBL     (input) INTEGER */
/*             The number of diagonal blocks of the matrix A, if it has a */
/*             block structure.  To specify that matrix A has no block */
/*             structure, set NBL = 0.  NBL >= 0. */

/*     NROWS   (input) INTEGER array, dimension max(1,NBL) */
/*             NROWS(i) contains the number of rows and columns of the */
/*             i-th diagonal block of matrix A.  The sum of the values */
/*             NROWS(i),  for  i = 1: NBL,  should be equal to min(M,N). */
/*             The array  NROWS  is not referenced if NBL = 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The matrix to be multiplied by CTO/CFROM.  See TYPE for */
/*             the storage type. */

/*     LDA     (input) INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,M). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             Not used in this implementation. */

/*     METHOD */

/*     Matrix A is multiplied by the real scalar CTO/CFROM, taking into */
/*     account the specified storage mode of the matrix. */
/*     MB01QD is a version of the LAPACK routine DLASCL, modified for */
/*     dealing with block triangular, or block Hessenberg matrices. */
/*     For efficiency, no tests of the input scalar parameters are */
/*     performed. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Nov. 1996. */

/*    ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 139 "MB01QD.f"
    /* Parameter adjustments */
#line 139 "MB01QD.f"
    --nrows;
#line 139 "MB01QD.f"
    a_dim1 = *lda;
#line 139 "MB01QD.f"
    a_offset = 1 + a_dim1;
#line 139 "MB01QD.f"
    a -= a_offset;
#line 139 "MB01QD.f"

#line 139 "MB01QD.f"
    /* Function Body */
#line 139 "MB01QD.f"
    if (lsame_(type__, "G", (ftnlen)1, (ftnlen)1)) {
#line 140 "MB01QD.f"
	itype = 0;
#line 141 "MB01QD.f"
    } else if (lsame_(type__, "L", (ftnlen)1, (ftnlen)1)) {
#line 142 "MB01QD.f"
	itype = 1;
#line 143 "MB01QD.f"
    } else if (lsame_(type__, "U", (ftnlen)1, (ftnlen)1)) {
#line 144 "MB01QD.f"
	itype = 2;
#line 145 "MB01QD.f"
    } else if (lsame_(type__, "H", (ftnlen)1, (ftnlen)1)) {
#line 146 "MB01QD.f"
	itype = 3;
#line 147 "MB01QD.f"
    } else if (lsame_(type__, "B", (ftnlen)1, (ftnlen)1)) {
#line 148 "MB01QD.f"
	itype = 4;
#line 149 "MB01QD.f"
    } else if (lsame_(type__, "Q", (ftnlen)1, (ftnlen)1)) {
#line 150 "MB01QD.f"
	itype = 5;
#line 151 "MB01QD.f"
    } else {
#line 152 "MB01QD.f"
	itype = 6;
#line 153 "MB01QD.f"
    }

/*     Quick return if possible. */

#line 157 "MB01QD.f"
    if (min(*m,*n) == 0) {
#line 157 "MB01QD.f"
	return 0;
#line 157 "MB01QD.f"
    }

/*     Get machine parameters. */

#line 162 "MB01QD.f"
    smlnum = dlamch_("S", (ftnlen)1);
#line 163 "MB01QD.f"
    bignum = 1. / smlnum;

#line 165 "MB01QD.f"
    cfromc = *cfrom;
#line 166 "MB01QD.f"
    ctoc = *cto;

#line 168 "MB01QD.f"
L10:
#line 169 "MB01QD.f"
    cfrom1 = cfromc * smlnum;
#line 170 "MB01QD.f"
    cto1 = ctoc / bignum;
#line 171 "MB01QD.f"
    if (abs(cfrom1) > abs(ctoc) && ctoc != 0.) {
#line 172 "MB01QD.f"
	mul = smlnum;
#line 173 "MB01QD.f"
	done = FALSE_;
#line 174 "MB01QD.f"
	cfromc = cfrom1;
#line 175 "MB01QD.f"
    } else if (abs(cto1) > abs(cfromc)) {
#line 176 "MB01QD.f"
	mul = bignum;
#line 177 "MB01QD.f"
	done = FALSE_;
#line 178 "MB01QD.f"
	ctoc = cto1;
#line 179 "MB01QD.f"
    } else {
#line 180 "MB01QD.f"
	mul = ctoc / cfromc;
#line 181 "MB01QD.f"
	done = TRUE_;
#line 182 "MB01QD.f"
    }

#line 184 "MB01QD.f"
    noblc = *nbl == 0;

#line 186 "MB01QD.f"
    if (itype == 0) {

/*        Full matrix */

#line 190 "MB01QD.f"
	i__1 = *n;
#line 190 "MB01QD.f"
	for (j = 1; j <= i__1; ++j) {
#line 191 "MB01QD.f"
	    i__2 = *m;
#line 191 "MB01QD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 192 "MB01QD.f"
		a[i__ + j * a_dim1] *= mul;
#line 193 "MB01QD.f"
/* L20: */
#line 193 "MB01QD.f"
	    }
#line 194 "MB01QD.f"
/* L30: */
#line 194 "MB01QD.f"
	}

#line 196 "MB01QD.f"
    } else if (itype == 1) {

#line 198 "MB01QD.f"
	if (noblc) {

/*           Lower triangular matrix */

#line 202 "MB01QD.f"
	    i__1 = *n;
#line 202 "MB01QD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 203 "MB01QD.f"
		i__2 = *m;
#line 203 "MB01QD.f"
		for (i__ = j; i__ <= i__2; ++i__) {
#line 204 "MB01QD.f"
		    a[i__ + j * a_dim1] *= mul;
#line 205 "MB01QD.f"
/* L40: */
#line 205 "MB01QD.f"
		}
#line 206 "MB01QD.f"
/* L50: */
#line 206 "MB01QD.f"
	    }

#line 208 "MB01QD.f"
	} else {

/*           Block lower triangular matrix */

#line 212 "MB01QD.f"
	    jfin = 0;
#line 213 "MB01QD.f"
	    i__1 = *nbl;
#line 213 "MB01QD.f"
	    for (k = 1; k <= i__1; ++k) {
#line 214 "MB01QD.f"
		jini = jfin + 1;
#line 215 "MB01QD.f"
		jfin += nrows[k];
#line 216 "MB01QD.f"
		i__2 = jfin;
#line 216 "MB01QD.f"
		for (j = jini; j <= i__2; ++j) {
#line 217 "MB01QD.f"
		    i__3 = *m;
#line 217 "MB01QD.f"
		    for (i__ = jini; i__ <= i__3; ++i__) {
#line 218 "MB01QD.f"
			a[i__ + j * a_dim1] *= mul;
#line 219 "MB01QD.f"
/* L60: */
#line 219 "MB01QD.f"
		    }
#line 220 "MB01QD.f"
/* L70: */
#line 220 "MB01QD.f"
		}
#line 221 "MB01QD.f"
/* L80: */
#line 221 "MB01QD.f"
	    }
#line 222 "MB01QD.f"
	}

#line 224 "MB01QD.f"
    } else if (itype == 2) {

#line 226 "MB01QD.f"
	if (noblc) {

/*           Upper triangular matrix */

#line 230 "MB01QD.f"
	    i__1 = *n;
#line 230 "MB01QD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 231 "MB01QD.f"
		i__2 = min(j,*m);
#line 231 "MB01QD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 232 "MB01QD.f"
		    a[i__ + j * a_dim1] *= mul;
#line 233 "MB01QD.f"
/* L90: */
#line 233 "MB01QD.f"
		}
#line 234 "MB01QD.f"
/* L100: */
#line 234 "MB01QD.f"
	    }

#line 236 "MB01QD.f"
	} else {

/*           Block upper triangular matrix */

#line 240 "MB01QD.f"
	    jfin = 0;
#line 241 "MB01QD.f"
	    i__1 = *nbl;
#line 241 "MB01QD.f"
	    for (k = 1; k <= i__1; ++k) {
#line 242 "MB01QD.f"
		jini = jfin + 1;
#line 243 "MB01QD.f"
		jfin += nrows[k];
#line 244 "MB01QD.f"
		if (k == *nbl) {
#line 244 "MB01QD.f"
		    jfin = *n;
#line 244 "MB01QD.f"
		}
#line 245 "MB01QD.f"
		i__2 = jfin;
#line 245 "MB01QD.f"
		for (j = jini; j <= i__2; ++j) {
#line 246 "MB01QD.f"
		    i__3 = min(jfin,*m);
#line 246 "MB01QD.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 247 "MB01QD.f"
			a[i__ + j * a_dim1] *= mul;
#line 248 "MB01QD.f"
/* L110: */
#line 248 "MB01QD.f"
		    }
#line 249 "MB01QD.f"
/* L120: */
#line 249 "MB01QD.f"
		}
#line 250 "MB01QD.f"
/* L130: */
#line 250 "MB01QD.f"
	    }
#line 251 "MB01QD.f"
	}

#line 253 "MB01QD.f"
    } else if (itype == 3) {

#line 255 "MB01QD.f"
	if (noblc) {

/*           Upper Hessenberg matrix */

#line 259 "MB01QD.f"
	    i__1 = *n;
#line 259 "MB01QD.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 260 "MB01QD.f"
		i__3 = j + 1;
#line 260 "MB01QD.f"
		i__2 = min(i__3,*m);
#line 260 "MB01QD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 261 "MB01QD.f"
		    a[i__ + j * a_dim1] *= mul;
#line 262 "MB01QD.f"
/* L140: */
#line 262 "MB01QD.f"
		}
#line 263 "MB01QD.f"
/* L150: */
#line 263 "MB01QD.f"
	    }

#line 265 "MB01QD.f"
	} else {

/*           Block upper Hessenberg matrix */

#line 269 "MB01QD.f"
	    jfin = 0;
#line 270 "MB01QD.f"
	    i__1 = *nbl;
#line 270 "MB01QD.f"
	    for (k = 1; k <= i__1; ++k) {
#line 271 "MB01QD.f"
		jini = jfin + 1;
#line 272 "MB01QD.f"
		jfin += nrows[k];

#line 274 "MB01QD.f"
		if (k == *nbl) {
#line 275 "MB01QD.f"
		    jfin = *n;
#line 276 "MB01QD.f"
		    ifin = *n;
#line 277 "MB01QD.f"
		} else {
#line 278 "MB01QD.f"
		    ifin = jfin + nrows[k + 1];
#line 279 "MB01QD.f"
		}

#line 281 "MB01QD.f"
		i__2 = jfin;
#line 281 "MB01QD.f"
		for (j = jini; j <= i__2; ++j) {
#line 282 "MB01QD.f"
		    i__3 = min(ifin,*m);
#line 282 "MB01QD.f"
		    for (i__ = 1; i__ <= i__3; ++i__) {
#line 283 "MB01QD.f"
			a[i__ + j * a_dim1] *= mul;
#line 284 "MB01QD.f"
/* L160: */
#line 284 "MB01QD.f"
		    }
#line 285 "MB01QD.f"
/* L170: */
#line 285 "MB01QD.f"
		}
#line 286 "MB01QD.f"
/* L180: */
#line 286 "MB01QD.f"
	    }
#line 287 "MB01QD.f"
	}

#line 289 "MB01QD.f"
    } else if (itype == 4) {

/*        Lower half of a symmetric band matrix */

#line 293 "MB01QD.f"
	k3 = *kl + 1;
#line 294 "MB01QD.f"
	k4 = *n + 1;
#line 295 "MB01QD.f"
	i__1 = *n;
#line 295 "MB01QD.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
#line 296 "MB01QD.f"
	    i__3 = k3, i__4 = k4 - j;
#line 296 "MB01QD.f"
	    i__2 = min(i__3,i__4);
#line 296 "MB01QD.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 297 "MB01QD.f"
		a[i__ + j * a_dim1] *= mul;
#line 298 "MB01QD.f"
/* L190: */
#line 298 "MB01QD.f"
	    }
#line 299 "MB01QD.f"
/* L200: */
#line 299 "MB01QD.f"
	}

#line 301 "MB01QD.f"
    } else if (itype == 5) {

/*        Upper half of a symmetric band matrix */

#line 305 "MB01QD.f"
	k1 = *ku + 2;
#line 306 "MB01QD.f"
	k3 = *ku + 1;
#line 307 "MB01QD.f"
	i__1 = *n;
#line 307 "MB01QD.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 308 "MB01QD.f"
	    i__2 = k1 - j;
#line 308 "MB01QD.f"
	    i__3 = k3;
#line 308 "MB01QD.f"
	    for (i__ = max(i__2,1); i__ <= i__3; ++i__) {
#line 309 "MB01QD.f"
		a[i__ + j * a_dim1] *= mul;
#line 310 "MB01QD.f"
/* L210: */
#line 310 "MB01QD.f"
	    }
#line 311 "MB01QD.f"
/* L220: */
#line 311 "MB01QD.f"
	}

#line 313 "MB01QD.f"
    } else if (itype == 6) {

/*        Band matrix */

#line 317 "MB01QD.f"
	k1 = *kl + *ku + 2;
#line 318 "MB01QD.f"
	k2 = *kl + 1;
#line 319 "MB01QD.f"
	k3 = (*kl << 1) + *ku + 1;
#line 320 "MB01QD.f"
	k4 = *kl + *ku + 1 + *m;
#line 321 "MB01QD.f"
	i__1 = *n;
#line 321 "MB01QD.f"
	for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 322 "MB01QD.f"
	    i__3 = k1 - j;
/* Computing MIN */
#line 322 "MB01QD.f"
	    i__4 = k3, i__5 = k4 - j;
#line 322 "MB01QD.f"
	    i__2 = min(i__4,i__5);
#line 322 "MB01QD.f"
	    for (i__ = max(i__3,k2); i__ <= i__2; ++i__) {
#line 323 "MB01QD.f"
		a[i__ + j * a_dim1] *= mul;
#line 324 "MB01QD.f"
/* L230: */
#line 324 "MB01QD.f"
	    }
#line 325 "MB01QD.f"
/* L240: */
#line 325 "MB01QD.f"
	}

#line 327 "MB01QD.f"
    }

#line 329 "MB01QD.f"
    if (! done) {
#line 329 "MB01QD.f"
	goto L10;
#line 329 "MB01QD.f"
    }

#line 332 "MB01QD.f"
    return 0;
/* *** Last line of MB01QD *** */
} /* mb01qd_ */

