#line 1 "MA02ID.f"
/* MA02ID.f -- translated by f2c (version 20100827).
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

#line 1 "MA02ID.f"
/* Table of constant values */

static integer c__1 = 1;

doublereal ma02id_(char *typ, char *norm, integer *n, doublereal *a, integer *
	lda, doublereal *qg, integer *ldqg, doublereal *dwork, ftnlen typ_len,
	 ftnlen norm_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, qg_dim1, qg_offset, i__1, i__2;
    doublereal ret_val, d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static logical lsh;
    static doublereal sum, dscl, temp, dsum, scale;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static doublereal value;
    extern doublereal dlapy2_(doublereal *, doublereal *), dlange_(char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    ftnlen);
    extern /* Subroutine */ int dlassq_(integer *, doublereal *, integer *, 
	    doublereal *, doublereal *);


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

/*     To compute the value of the one norm, or the Frobenius norm, or */
/*     the infinity norm, or the element of largest absolute value */
/*     of a real skew-Hamiltonian matrix */

/*                   [  A   G  ]          T         T */
/*             X  =  [       T ],   G = -G,   Q = -Q, */
/*                   [  Q   A  ] */

/*     or of a real Hamiltonian matrix */

/*                   [  A   G  ]          T         T */
/*             X  =  [       T ],   G =  G,   Q =  Q, */
/*                   [  Q  -A  ] */

/*     where A, G and Q are real n-by-n matrices. */

/*     Note that for this kind of matrices the infinity norm is equal */
/*     to the one norm. */

/*     FUNCTION VALUE */

/*     MA02ID  DOUBLE PRECISION */
/*             The computed norm. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYP     CHARACTER*1 */
/*             Specifies the type of the input matrix X: */
/*             = 'S':         X is skew-Hamiltonian; */
/*             = 'H':         X is Hamiltonian. */

/*     NORM    CHARACTER*1 */
/*             Specifies the value to be returned in MA02ID: */
/*             = '1' or 'O':  one norm of X; */
/*             = 'F' or 'E':  Frobenius norm of X; */
/*             = 'I':         infinity norm of X; */
/*             = 'M':         max(abs(X(i,j)). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     QG      (input) DOUBLE PRECISION array, dimension (LDQG,N+1) */
/*             On entry, the leading N-by-N+1 part of this array must */
/*             contain in columns 1:N the lower triangular part of the */
/*             matrix Q and in columns 2:N+1 the upper triangular part */
/*             of the matrix G. If TYP = 'S', the parts containing the */
/*             diagonal and the first supdiagonal of this array are not */
/*             referenced. */

/*     LDQG    INTEGER */
/*             The leading dimension of the array QG.  LDQG >= MAX(1,N). */

/*     Workspace */


/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             where LDWORK >= 2*N when NORM = '1', NORM = 'I' or */
/*             NORM = 'O'; otherwise, DWORK is not referenced. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DLANHA). */

/*     KEYWORDS */

/*     Elementary matrix operations, Hamiltonian matrix, skew-Hamiltonian */
/*     matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

#line 133 "MA02ID.f"
    /* Parameter adjustments */
#line 133 "MA02ID.f"
    a_dim1 = *lda;
#line 133 "MA02ID.f"
    a_offset = 1 + a_dim1;
#line 133 "MA02ID.f"
    a -= a_offset;
#line 133 "MA02ID.f"
    qg_dim1 = *ldqg;
#line 133 "MA02ID.f"
    qg_offset = 1 + qg_dim1;
#line 133 "MA02ID.f"
    qg -= qg_offset;
#line 133 "MA02ID.f"
    --dwork;
#line 133 "MA02ID.f"

#line 133 "MA02ID.f"
    /* Function Body */
#line 133 "MA02ID.f"
    lsh = lsame_(typ, "S", (ftnlen)1, (ftnlen)1);

#line 135 "MA02ID.f"
    if (*n == 0) {
#line 136 "MA02ID.f"
	value = 0.;

#line 138 "MA02ID.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1) && lsh) {

/*        Find max(abs(A(i,j))). */

#line 142 "MA02ID.f"
	value = dlange_("MaxElement", n, n, &a[a_offset], lda, &dwork[1], (
		ftnlen)10);
#line 143 "MA02ID.f"
	if (*n > 1) {
#line 144 "MA02ID.f"
	    i__1 = *n + 1;
#line 144 "MA02ID.f"
	    for (j = 1; j <= i__1; ++j) {
#line 145 "MA02ID.f"
		i__2 = j - 2;
#line 145 "MA02ID.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 146 "MA02ID.f"
		    d__2 = value, d__3 = (d__1 = qg[i__ + j * qg_dim1], abs(
			    d__1));
#line 146 "MA02ID.f"
		    value = max(d__2,d__3);
#line 147 "MA02ID.f"
/* L10: */
#line 147 "MA02ID.f"
		}
#line 148 "MA02ID.f"
		i__2 = *n;
#line 148 "MA02ID.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
/* Computing MAX */
#line 149 "MA02ID.f"
		    d__2 = value, d__3 = (d__1 = qg[i__ + j * qg_dim1], abs(
			    d__1));
#line 149 "MA02ID.f"
		    value = max(d__2,d__3);
#line 150 "MA02ID.f"
/* L20: */
#line 150 "MA02ID.f"
		}
#line 151 "MA02ID.f"
/* L30: */
#line 151 "MA02ID.f"
	    }
#line 152 "MA02ID.f"
	}

#line 154 "MA02ID.f"
    } else if (lsame_(norm, "M", (ftnlen)1, (ftnlen)1)) {

/*        Find max( abs( A(i,j) ), abs( QG(i,j) ) ). */

/* Computing MAX */
#line 158 "MA02ID.f"
	i__1 = *n + 1;
#line 158 "MA02ID.f"
	d__1 = dlange_("MaxElement", n, n, &a[a_offset], lda, &dwork[1], (
		ftnlen)10), d__2 = dlange_("MaxElement", n, &i__1, &qg[
		qg_offset], ldqg, &dwork[1], (ftnlen)10);
#line 158 "MA02ID.f"
	value = max(d__1,d__2);

#line 162 "MA02ID.f"
    } else if ((lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1' || lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) && lsh) {

/*        Find the column and row sums of A (in one pass). */

#line 167 "MA02ID.f"
	value = 0.;
#line 168 "MA02ID.f"
	i__1 = *n;
#line 168 "MA02ID.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 169 "MA02ID.f"
	    dwork[i__] = 0.;
#line 170 "MA02ID.f"
/* L40: */
#line 170 "MA02ID.f"
	}

#line 172 "MA02ID.f"
	i__1 = *n;
#line 172 "MA02ID.f"
	for (j = 1; j <= i__1; ++j) {
#line 173 "MA02ID.f"
	    sum = 0.;
#line 174 "MA02ID.f"
	    i__2 = *n;
#line 174 "MA02ID.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 175 "MA02ID.f"
		temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 176 "MA02ID.f"
		sum += temp;
#line 177 "MA02ID.f"
		dwork[i__] += temp;
#line 178 "MA02ID.f"
/* L50: */
#line 178 "MA02ID.f"
	    }
#line 179 "MA02ID.f"
	    dwork[*n + j] = sum;
#line 180 "MA02ID.f"
/* L60: */
#line 180 "MA02ID.f"
	}

/*        Compute the maximal absolute column sum. */

#line 184 "MA02ID.f"
	i__1 = *n + 1;
#line 184 "MA02ID.f"
	for (j = 1; j <= i__1; ++j) {
#line 185 "MA02ID.f"
	    i__2 = j - 2;
#line 185 "MA02ID.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 186 "MA02ID.f"
		temp = (d__1 = qg[i__ + j * qg_dim1], abs(d__1));
#line 187 "MA02ID.f"
		dwork[i__] += temp;
#line 188 "MA02ID.f"
		dwork[j - 1] += temp;
#line 189 "MA02ID.f"
/* L70: */
#line 189 "MA02ID.f"
	    }
#line 190 "MA02ID.f"
	    if (j < *n + 1) {
#line 191 "MA02ID.f"
		sum = dwork[*n + j];
#line 192 "MA02ID.f"
		i__2 = *n;
#line 192 "MA02ID.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 193 "MA02ID.f"
		    temp = (d__1 = qg[i__ + j * qg_dim1], abs(d__1));
#line 194 "MA02ID.f"
		    sum += temp;
#line 195 "MA02ID.f"
		    dwork[*n + i__] += temp;
#line 196 "MA02ID.f"
/* L80: */
#line 196 "MA02ID.f"
		}
#line 197 "MA02ID.f"
		value = max(value,sum);
#line 198 "MA02ID.f"
	    }
#line 199 "MA02ID.f"
/* L90: */
#line 199 "MA02ID.f"
	}
#line 200 "MA02ID.f"
	i__1 = *n;
#line 200 "MA02ID.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 201 "MA02ID.f"
	    d__1 = value, d__2 = dwork[i__];
#line 201 "MA02ID.f"
	    value = max(d__1,d__2);
#line 202 "MA02ID.f"
/* L100: */
#line 202 "MA02ID.f"
	}

#line 204 "MA02ID.f"
    } else if (lsame_(norm, "O", (ftnlen)1, (ftnlen)1) || *(unsigned char *)
	    norm == '1' || lsame_(norm, "I", (ftnlen)1, (ftnlen)1)) {

/*        Find the column and row sums of A (in one pass). */

#line 209 "MA02ID.f"
	value = 0.;
#line 210 "MA02ID.f"
	i__1 = *n;
#line 210 "MA02ID.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 211 "MA02ID.f"
	    dwork[i__] = 0.;
#line 212 "MA02ID.f"
/* L110: */
#line 212 "MA02ID.f"
	}

#line 214 "MA02ID.f"
	i__1 = *n;
#line 214 "MA02ID.f"
	for (j = 1; j <= i__1; ++j) {
#line 215 "MA02ID.f"
	    sum = 0.;
#line 216 "MA02ID.f"
	    i__2 = *n;
#line 216 "MA02ID.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 217 "MA02ID.f"
		temp = (d__1 = a[i__ + j * a_dim1], abs(d__1));
#line 218 "MA02ID.f"
		sum += temp;
#line 219 "MA02ID.f"
		dwork[i__] += temp;
#line 220 "MA02ID.f"
/* L120: */
#line 220 "MA02ID.f"
	    }
#line 221 "MA02ID.f"
	    dwork[*n + j] = sum;
#line 222 "MA02ID.f"
/* L130: */
#line 222 "MA02ID.f"
	}

/*        Compute the maximal absolute column sum. */

#line 226 "MA02ID.f"
	i__1 = *n + 1;
#line 226 "MA02ID.f"
	for (j = 1; j <= i__1; ++j) {
#line 227 "MA02ID.f"
	    i__2 = j - 2;
#line 227 "MA02ID.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 228 "MA02ID.f"
		temp = (d__1 = qg[i__ + j * qg_dim1], abs(d__1));
#line 229 "MA02ID.f"
		dwork[i__] += temp;
#line 230 "MA02ID.f"
		dwork[j - 1] += temp;
#line 231 "MA02ID.f"
/* L140: */
#line 231 "MA02ID.f"
	    }
#line 232 "MA02ID.f"
	    if (j > 1) {
#line 232 "MA02ID.f"
		dwork[j - 1] += (d__1 = qg[j - 1 + j * qg_dim1], abs(d__1));
#line 232 "MA02ID.f"
	    }
#line 234 "MA02ID.f"
	    if (j < *n + 1) {
#line 235 "MA02ID.f"
		sum = dwork[*n + j] + (d__1 = qg[j + j * qg_dim1], abs(d__1));
#line 236 "MA02ID.f"
		i__2 = *n;
#line 236 "MA02ID.f"
		for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 237 "MA02ID.f"
		    temp = (d__1 = qg[i__ + j * qg_dim1], abs(d__1));
#line 238 "MA02ID.f"
		    sum += temp;
#line 239 "MA02ID.f"
		    dwork[*n + i__] += temp;
#line 240 "MA02ID.f"
/* L150: */
#line 240 "MA02ID.f"
		}
#line 241 "MA02ID.f"
		value = max(value,sum);
#line 242 "MA02ID.f"
	    }
#line 243 "MA02ID.f"
/* L160: */
#line 243 "MA02ID.f"
	}
#line 244 "MA02ID.f"
	i__1 = *n;
#line 244 "MA02ID.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 245 "MA02ID.f"
	    d__1 = value, d__2 = dwork[i__];
#line 245 "MA02ID.f"
	    value = max(d__1,d__2);
#line 246 "MA02ID.f"
/* L170: */
#line 246 "MA02ID.f"
	}

#line 248 "MA02ID.f"
    } else if ((lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) && lsh) {

/*        Find normF(A). */

#line 253 "MA02ID.f"
	scale = 0.;
#line 254 "MA02ID.f"
	sum = 1.;
#line 255 "MA02ID.f"
	i__1 = *n;
#line 255 "MA02ID.f"
	for (j = 1; j <= i__1; ++j) {
#line 256 "MA02ID.f"
	    dlassq_(n, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 257 "MA02ID.f"
/* L180: */
#line 257 "MA02ID.f"
	}

/*        Add normF(G) and normF(Q). */

#line 261 "MA02ID.f"
	i__1 = *n + 1;
#line 261 "MA02ID.f"
	for (j = 1; j <= i__1; ++j) {
#line 262 "MA02ID.f"
	    if (j > 2) {
#line 262 "MA02ID.f"
		i__2 = j - 2;
#line 262 "MA02ID.f"
		dlassq_(&i__2, &qg[j * qg_dim1 + 1], &c__1, &scale, &sum);
#line 262 "MA02ID.f"
	    }
#line 264 "MA02ID.f"
	    if (j < *n) {
#line 264 "MA02ID.f"
		i__2 = *n - j;
#line 264 "MA02ID.f"
		dlassq_(&i__2, &qg[j + 1 + j * qg_dim1], &c__1, &scale, &sum);
#line 264 "MA02ID.f"
	    }
#line 266 "MA02ID.f"
/* L190: */
#line 266 "MA02ID.f"
	}
#line 267 "MA02ID.f"
	value = sqrt(2.) * scale * sqrt(sum);
#line 268 "MA02ID.f"
    } else if (lsame_(norm, "F", (ftnlen)1, (ftnlen)1) || lsame_(norm, "E", (
	    ftnlen)1, (ftnlen)1)) {
#line 269 "MA02ID.f"
	scale = 0.;
#line 270 "MA02ID.f"
	sum = 1.;
#line 271 "MA02ID.f"
	i__1 = *n;
#line 271 "MA02ID.f"
	for (j = 1; j <= i__1; ++j) {
#line 272 "MA02ID.f"
	    dlassq_(n, &a[j * a_dim1 + 1], &c__1, &scale, &sum);
#line 273 "MA02ID.f"
/* L200: */
#line 273 "MA02ID.f"
	}
#line 274 "MA02ID.f"
	dscl = 0.;
#line 275 "MA02ID.f"
	dsum = 1.;
#line 276 "MA02ID.f"
	i__1 = *n + 1;
#line 276 "MA02ID.f"
	for (j = 1; j <= i__1; ++j) {
#line 277 "MA02ID.f"
	    if (j > 1) {
#line 278 "MA02ID.f"
		i__2 = j - 2;
#line 278 "MA02ID.f"
		dlassq_(&i__2, &qg[j * qg_dim1 + 1], &c__1, &scale, &sum);
#line 279 "MA02ID.f"
		dlassq_(&c__1, &qg[j - 1 + j * qg_dim1], &c__1, &dscl, &dsum);
#line 280 "MA02ID.f"
	    }
#line 281 "MA02ID.f"
	    if (j < *n + 1) {
#line 282 "MA02ID.f"
		dlassq_(&c__1, &qg[j + j * qg_dim1], &c__1, &dscl, &dsum);
#line 283 "MA02ID.f"
		i__2 = *n - j;
#line 283 "MA02ID.f"
		dlassq_(&i__2, &qg[j + 1 + j * qg_dim1], &c__1, &scale, &sum);
#line 284 "MA02ID.f"
	    }
#line 285 "MA02ID.f"
/* L210: */
#line 285 "MA02ID.f"
	}
#line 286 "MA02ID.f"
	d__1 = sqrt(2.) * scale * sqrt(sum);
#line 286 "MA02ID.f"
	d__2 = dscl * sqrt(dsum);
#line 286 "MA02ID.f"
	value = dlapy2_(&d__1, &d__2);
#line 288 "MA02ID.f"
    }

#line 290 "MA02ID.f"
    ret_val = value;
#line 291 "MA02ID.f"
    return ret_val;
/* *** Last line of MA02ID *** */
} /* ma02id_ */

