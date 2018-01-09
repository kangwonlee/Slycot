#line 1 "NF01BV.f"
/* NF01BV.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BV.f"
/* Table of constant values */

static doublereal c_b7 = 0.;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b14 = 1.;

/* Subroutine */ int nf01bv_(char *stor, char *uplo, integer *n, integer *
	ipar, integer *lipar, doublereal *dpar, integer *ldpar, doublereal *j,
	 integer *ldj, doublereal *jtj, integer *ldjtj, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen stor_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer j_dim1, j_offset, i__1;

    /* Local variables */
    static doublereal c__;
    static integer i__, m, ii;
    static doublereal dum[1];
    static logical full;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *);
    static logical upper;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), dlaset_(char *, integer *, integer *,
	     doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);


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

/*     To compute the matrix J'*J + c*I, for the Jacobian J as received */
/*     from SLICOT Library routine NF01BY, for one output variable. */

/*     NOTE: this routine must have the same arguments as SLICOT Library */
/*     routine NF01BU. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     STOR    CHARACTER*1 */
/*             Specifies the storage scheme for the symmetric */
/*             matrix J'*J + c*I, as follows: */
/*             = 'F' :  full storage is used; */
/*             = 'P' :  packed storage is used. */

/*     UPLO    CHARACTER*1 */
/*             Specifies which part of the matrix J'*J + c*I is stored, */
/*             as follows: */
/*             = 'U' :  the upper triagular part is stored; */
/*             = 'L' :  the lower triagular part is stored. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The number of columns of the Jacobian matrix J.  N >= 0. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters describing the structure of the */
/*             matrix J, as follows: */
/*             IPAR(1) must contain the number of rows M of the Jacobian */
/*                     matrix J.  M >= 0. */
/*             IPAR is provided for compatibility with SLICOT Library */
/*             routine MD03AD. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 1. */

/*     DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*             The real parameters needed for solving the problem. */
/*             The entry DPAR(1) must contain the real scalar c. */

/*     LDPAR   (input) INTEGER */
/*             The length of the array DPAR.  LDPAR >= 1. */

/*     J       (input) DOUBLE PRECISION array, dimension (LDJ,N) */
/*             The leading M-by-N part of this array must contain the */
/*             Jacobian matrix J. */

/*     LDJ     INTEGER */
/*             The leading dimension of the array J.  LDJ >= MAX(1,M). */

/*     JTJ     (output) DOUBLE PRECISION array, */
/*                      dimension (LDJTJ,N),    if STOR = 'F', */
/*                      dimension (N*(N+1)/2),  if STOR = 'P'. */
/*             The leading N-by-N (if STOR = 'F'), or N*(N+1)/2 (if */
/*             STOR = 'P') part of this array contains the upper or */
/*             lower triangle of the matrix J'*J + c*I, depending on */
/*             UPLO = 'U', or UPLO = 'L', respectively, stored either as */
/*             a two-dimensional, or one-dimensional array, depending */
/*             on STOR. */

/*     LDJTJ   INTEGER */
/*             The leading dimension of the array JTJ. */
/*             LDJTJ >= MAX(1,N), if STOR = 'F'. */
/*             LDJTJ >= 1,        if STOR = 'P'. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             Currently, this array is not used. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The matrix product is computed columnn-wise, exploiting the */
/*     symmetry. BLAS 3 routine DSYRK is used if STOR = 'F', and BLAS 2 */
/*     routine DGEMV is used if STOR = 'P'. */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2001. */

/*     REVISIONS */

/*     V. Sima, March 2002. */

/*     KEYWORDS */

/*     Elementary matrix operations, matrix algebra, matrix operations, */
/*     Wiener system. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 152 "NF01BV.f"
    /* Parameter adjustments */
#line 152 "NF01BV.f"
    --ipar;
#line 152 "NF01BV.f"
    --dpar;
#line 152 "NF01BV.f"
    j_dim1 = *ldj;
#line 152 "NF01BV.f"
    j_offset = 1 + j_dim1;
#line 152 "NF01BV.f"
    j -= j_offset;
#line 152 "NF01BV.f"
    --jtj;
#line 152 "NF01BV.f"
    --dwork;
#line 152 "NF01BV.f"

#line 152 "NF01BV.f"
    /* Function Body */
#line 152 "NF01BV.f"
    *info = 0;
#line 153 "NF01BV.f"
    full = lsame_(stor, "F", (ftnlen)1, (ftnlen)1);
#line 154 "NF01BV.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

#line 156 "NF01BV.f"
    if (! (full || lsame_(stor, "P", (ftnlen)1, (ftnlen)1))) {
#line 157 "NF01BV.f"
	*info = -1;
#line 158 "NF01BV.f"
    } else if (! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 159 "NF01BV.f"
	*info = -2;
#line 160 "NF01BV.f"
    } else if (*n < 0) {
#line 161 "NF01BV.f"
	*info = -3;
#line 162 "NF01BV.f"
    } else if (*lipar < 1) {
#line 163 "NF01BV.f"
	*info = -5;
#line 164 "NF01BV.f"
    } else if (*ldpar < 1) {
#line 165 "NF01BV.f"
	*info = -7;
#line 166 "NF01BV.f"
    } else if (*ldjtj < 1 || full && *ldjtj < *n) {
#line 167 "NF01BV.f"
	*info = -11;
#line 168 "NF01BV.f"
    } else if (*ldwork < 0) {
#line 169 "NF01BV.f"
	*info = -13;
#line 170 "NF01BV.f"
    } else {
#line 171 "NF01BV.f"
	m = ipar[1];
#line 172 "NF01BV.f"
	if (m < 0) {
#line 173 "NF01BV.f"
	    *info = -4;
#line 174 "NF01BV.f"
	} else if (*ldj < max(1,m)) {
#line 175 "NF01BV.f"
	    *info = -9;
#line 176 "NF01BV.f"
	}
#line 177 "NF01BV.f"
    }

/*     Return if there are illegal arguments. */

#line 181 "NF01BV.f"
    if (*info != 0) {
#line 182 "NF01BV.f"
	i__1 = -(*info);
#line 182 "NF01BV.f"
	xerbla_("NF01BV", &i__1, (ftnlen)6);
#line 183 "NF01BV.f"
	return 0;
#line 184 "NF01BV.f"
    }

/*     Quick return if possible. */

#line 188 "NF01BV.f"
    c__ = dpar[1];
#line 189 "NF01BV.f"
    if (*n == 0) {
#line 190 "NF01BV.f"
	return 0;
#line 191 "NF01BV.f"
    } else if (m == 0) {
#line 192 "NF01BV.f"
	if (full) {
#line 193 "NF01BV.f"
	    dlaset_(uplo, n, n, &c_b7, &c__, &jtj[1], ldjtj, (ftnlen)1);
#line 194 "NF01BV.f"
	} else {
#line 195 "NF01BV.f"
	    dum[0] = 0.;
#line 196 "NF01BV.f"
	    i__1 = *n * (*n + 1) / 2;
#line 196 "NF01BV.f"
	    dcopy_(&i__1, dum, &c__0, &jtj[1], &c__1);
#line 197 "NF01BV.f"
	    if (upper) {
#line 198 "NF01BV.f"
		ii = 0;

#line 200 "NF01BV.f"
		i__1 = *n;
#line 200 "NF01BV.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 201 "NF01BV.f"
		    ii += i__;
#line 202 "NF01BV.f"
		    jtj[ii] = c__;
#line 203 "NF01BV.f"
/* L10: */
#line 203 "NF01BV.f"
		}

#line 205 "NF01BV.f"
	    } else {
#line 206 "NF01BV.f"
		ii = 1;

#line 208 "NF01BV.f"
		for (i__ = *n; i__ >= 1; --i__) {
#line 209 "NF01BV.f"
		    jtj[ii] = c__;
#line 210 "NF01BV.f"
		    ii += i__;
#line 211 "NF01BV.f"
/* L20: */
#line 211 "NF01BV.f"
		}

#line 213 "NF01BV.f"
	    }
#line 214 "NF01BV.f"
	}
#line 215 "NF01BV.f"
	return 0;
#line 216 "NF01BV.f"
    }

/*     Build a triangle of the matrix J'*J + c*I. */

#line 220 "NF01BV.f"
    if (full) {
#line 221 "NF01BV.f"
	dlaset_(uplo, n, n, &c_b7, &c__, &jtj[1], ldjtj, (ftnlen)1);
#line 222 "NF01BV.f"
	dsyrk_(uplo, "Transpose", n, &m, &c_b14, &j[j_offset], ldj, &c_b14, &
		jtj[1], ldjtj, (ftnlen)1, (ftnlen)9);
#line 224 "NF01BV.f"
    } else if (upper) {
#line 225 "NF01BV.f"
	ii = 0;

#line 227 "NF01BV.f"
	i__1 = *n;
#line 227 "NF01BV.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 228 "NF01BV.f"
	    dgemv_("Transpose", &m, &i__, &c_b14, &j[j_offset], ldj, &j[i__ * 
		    j_dim1 + 1], &c__1, &c_b7, &jtj[ii + 1], &c__1, (ftnlen)9)
		    ;
#line 230 "NF01BV.f"
	    ii += i__;
#line 231 "NF01BV.f"
	    jtj[ii] += c__;
#line 232 "NF01BV.f"
/* L30: */
#line 232 "NF01BV.f"
	}

#line 234 "NF01BV.f"
    } else {
#line 235 "NF01BV.f"
	ii = 1;

#line 237 "NF01BV.f"
	for (i__ = *n; i__ >= 1; --i__) {
#line 238 "NF01BV.f"
	    dgemv_("Transpose", &m, &i__, &c_b14, &j[(*n - i__ + 1) * j_dim1 
		    + 1], ldj, &j[(*n - i__ + 1) * j_dim1 + 1], &c__1, &c_b7, 
		    &jtj[ii], &c__1, (ftnlen)9);
#line 240 "NF01BV.f"
	    jtj[ii] += c__;
#line 241 "NF01BV.f"
	    ii += i__;
#line 242 "NF01BV.f"
/* L40: */
#line 242 "NF01BV.f"
	}

#line 244 "NF01BV.f"
    }

#line 246 "NF01BV.f"
    return 0;

/* *** Last line of NF01BV *** */
} /* nf01bv_ */

