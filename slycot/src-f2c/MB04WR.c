#line 1 "MB04WR.f"
/* MB04WR.f -- translated by f2c (version 20100827).
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

#line 1 "MB04WR.f"
/* Table of constant values */

static doublereal c_b9 = 0.;
static doublereal c_b10 = 1.;

/* Subroutine */ int mb04wr_(char *job, char *trans, integer *n, integer *ilo,
	 doublereal *q1, integer *ldq1, doublereal *q2, integer *ldq2, 
	doublereal *cs, doublereal *tau, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen job_len, ftnlen trans_len)
{
    /* System generated locals */
    integer q1_dim1, q1_offset, q2_dim1, q2_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, nh, ierr;
    extern /* Subroutine */ int mb04wd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical ltran, compu;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
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

/*     To generate orthogonal symplectic matrices U or V, defined as */
/*     products of symplectic reflectors and Givens rotators */

/*     U = diag( HU(1),HU(1) )  GU(1)  diag( FU(1),FU(1) ) */
/*         diag( HU(2),HU(2) )  GU(2)  diag( FU(2),FU(2) ) */
/*                              .... */
/*         diag( HU(n),HU(n) )  GU(n)  diag( FU(n),FU(n) ), */

/*     V = diag( HV(1),HV(1) )       GV(1)   diag( FV(1),FV(1) ) */
/*         diag( HV(2),HV(2) )       GV(2)   diag( FV(2),FV(2) ) */
/*                                   .... */
/*         diag( HV(n-1),HV(n-1) )  GV(n-1)  diag( FV(n-1),FV(n-1) ), */

/*     as returned by the SLICOT Library routines MB04TS or MB04TB. The */
/*     matrices U and V are returned in terms of their first N/2 rows: */

/*                 [  U1   U2 ]           [  V1   V2 ] */
/*             U = [          ],      V = [          ]. */
/*                 [ -U2   U1 ]           [ -V2   V1 ] */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     JOB     CHARACTER*1 */
/*             Specifies whether the matrix U or the matrix V is */
/*             required: */
/*             = 'U':  generate U; */
/*             = 'V':  generate V. */

/*     TRANS   CHARACTER*1 */
/*             If  JOB = 'U'  then TRANS must have the same value as */
/*             the argument TRANA in the previous call of MB04TS or */
/*             MB04TB. */
/*             If  JOB = 'V'  then TRANS must have the same value as */
/*             the argument TRANB in the previous call of MB04TS or */
/*             MB04TB. */

/*     N       (input) INTEGER */
/*             The order of the matrices Q1 and Q2. N >= 0. */

/*     ILO     (input) INTEGER */
/*             ILO must have the same value as in the previous call of */
/*             MB04TS or MB04TB. U and V are equal to the unit matrix */
/*             except in the submatrices */
/*             U([ilo:n n+ilo:2*n], [ilo:n n+ilo:2*n]) and */
/*             V([ilo+1:n n+ilo+1:2*n], [ilo+1:n n+ilo+1:2*n]), */
/*             respectively. */
/*             1 <= ILO <= N, if N > 0; ILO = 1, if N = 0. */

/*     Q1      (input/output) DOUBLE PRECISION array, dimension (LDQ1,N) */
/*             On entry, if  JOB = 'U'  and  TRANS = 'N'  then the */
/*             leading N-by-N part of this array must contain in its i-th */
/*             column the vector which defines the elementary reflector */
/*             FU(i). */
/*             If  JOB = 'U'  and  TRANS = 'T'  or  TRANS = 'C' then the */
/*             leading N-by-N part of this array must contain in its i-th */
/*             row the vector which defines the elementary reflector */
/*             FU(i). */
/*             If  JOB = 'V'  and  TRANS = 'N'  then the leading N-by-N */
/*             part of this array must contain in its i-th row the vector */
/*             which defines the elementary reflector FV(i). */
/*             If  JOB = 'V'  and  TRANS = 'T'  or  TRANS = 'C' then the */
/*             leading N-by-N part of this array must contain in its i-th */
/*             column the vector which defines the elementary reflector */
/*             FV(i). */
/*             On exit, if  JOB = 'U'  and  TRANS = 'N'  then the leading */
/*             N-by-N part of this array contains the matrix U1. */
/*             If  JOB = 'U'  and  TRANS = 'T'  or  TRANS = 'C' then the */
/*             leading N-by-N part of this array contains the matrix */
/*             U1**T. */
/*             If  JOB = 'V'  and  TRANS = 'N'  then the leading N-by-N */
/*             part of this array contains the matrix V1**T. */
/*             If  JOB = 'V'  and  TRANS = 'T'  or  TRANS = 'C' then the */
/*             leading N-by-N part of this array contains the matrix V1. */

/*     LDQ1    INTEGER */
/*             The leading dimension of the array Q1.  LDQ1 >= MAX(1,N). */

/*     Q2      (input/output) DOUBLE PRECISION array, dimension (LDQ2,N) */
/*             On entry, if  JOB = 'U'  then the leading N-by-N part of */
/*             this array must contain in its i-th column the vector */
/*             which defines the elementary reflector HU(i). */
/*             If  JOB = 'V'  then the leading N-by-N part of this array */
/*             must contain in its i-th row the vector which defines the */
/*             elementary reflector HV(i). */
/*             On exit, if  JOB = 'U'  then the leading N-by-N part of */
/*             this array contains the matrix U2. */
/*             If  JOB = 'V'  then the leading N-by-N part of this array */
/*             contains the matrix V2**T. */

/*     LDQ2    INTEGER */
/*             The leading dimension of the array Q2.  LDQ2 >= MAX(1,N). */

/*     CS      (input) DOUBLE PRECISION array, dimension (2N) */
/*             On entry, if  JOB = 'U'  then the first 2N elements of */
/*             this array must contain the cosines and sines of the */
/*             symplectic Givens rotators GU(i). */
/*             If  JOB = 'V'  then the first 2N-2 elements of this array */
/*             must contain the cosines and sines of the symplectic */
/*             Givens rotators GV(i). */

/*     TAU     (input) DOUBLE PRECISION array, dimension (N) */
/*             On entry, if  JOB = 'U'  then the first N elements of */
/*             this array must contain the scalar factors of the */
/*             elementary reflectors FU(i). */
/*             If  JOB = 'V'  then the first N-1 elements of this array */
/*             must contain the scalar factors of the elementary */
/*             reflectors FV(i). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -12,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,2*(N-ILO+1)). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     REFERENCES */

/*     [1] Benner, P., Mehrmann, V., and Xu, H. */
/*         A numerically stable, structure preserving method for */
/*         computing the eigenvalues of real Hamiltonian or symplectic */
/*         pencils. Numer. Math., Vol 78 (3), pp. 329-358, 1998. */

/*     [2] Kressner, D. */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine DOSGSU). */

/*     KEYWORDS */

/*     Elementary matrix operations, Hamiltonian matrix, orthogonal */
/*     symplectic matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 203 "MB04WR.f"
    /* Parameter adjustments */
#line 203 "MB04WR.f"
    q1_dim1 = *ldq1;
#line 203 "MB04WR.f"
    q1_offset = 1 + q1_dim1;
#line 203 "MB04WR.f"
    q1 -= q1_offset;
#line 203 "MB04WR.f"
    q2_dim1 = *ldq2;
#line 203 "MB04WR.f"
    q2_offset = 1 + q2_dim1;
#line 203 "MB04WR.f"
    q2 -= q2_offset;
#line 203 "MB04WR.f"
    --cs;
#line 203 "MB04WR.f"
    --tau;
#line 203 "MB04WR.f"
    --dwork;
#line 203 "MB04WR.f"

#line 203 "MB04WR.f"
    /* Function Body */
#line 203 "MB04WR.f"
    *info = 0;
#line 204 "MB04WR.f"
    ltran = lsame_(trans, "T", (ftnlen)1, (ftnlen)1) || lsame_(trans, "C", (
	    ftnlen)1, (ftnlen)1);
#line 205 "MB04WR.f"
    compu = lsame_(job, "U", (ftnlen)1, (ftnlen)1);
#line 206 "MB04WR.f"
    if (! compu && ! lsame_(job, "V", (ftnlen)1, (ftnlen)1)) {
#line 207 "MB04WR.f"
	*info = -1;
#line 208 "MB04WR.f"
    } else if (! ltran && ! lsame_(trans, "N", (ftnlen)1, (ftnlen)1)) {
#line 209 "MB04WR.f"
	*info = -2;
#line 210 "MB04WR.f"
    } else if (*n < 0) {
#line 211 "MB04WR.f"
	*info = -3;
#line 212 "MB04WR.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 213 "MB04WR.f"
	*info = -4;
#line 214 "MB04WR.f"
    } else if (*ldq1 < max(1,*n)) {
#line 215 "MB04WR.f"
	*info = -6;
#line 216 "MB04WR.f"
    } else if (*ldq2 < max(1,*n)) {
#line 217 "MB04WR.f"
	*info = -8;
#line 218 "MB04WR.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 218 "MB04WR.f"
	i__1 = 1, i__2 = *n - *ilo + 1 << 1;
#line 218 "MB04WR.f"
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
#line 219 "MB04WR.f"
	    i__1 = 1, i__2 = *n - *ilo + 1 << 1;
#line 219 "MB04WR.f"
	    dwork[1] = (doublereal) max(i__1,i__2);
#line 220 "MB04WR.f"
	    *info = -12;
#line 221 "MB04WR.f"
	}
#line 221 "MB04WR.f"
    }

/*     Return if there were illegal values. */

#line 225 "MB04WR.f"
    if (*info != 0) {
#line 226 "MB04WR.f"
	i__1 = -(*info);
#line 226 "MB04WR.f"
	xerbla_("MB04WR", &i__1, (ftnlen)6);
#line 227 "MB04WR.f"
	return 0;
#line 228 "MB04WR.f"
    }

/*     Quick return if possible. */

#line 232 "MB04WR.f"
    if (*n == 0) {
#line 233 "MB04WR.f"
	dwork[1] = 1.;
#line 234 "MB04WR.f"
	return 0;
#line 235 "MB04WR.f"
    }

#line 237 "MB04WR.f"
    if (compu) {
#line 238 "MB04WR.f"
	i__1 = *ilo - 1;
#line 238 "MB04WR.f"
	dlaset_("All", n, &i__1, &c_b9, &c_b10, &q1[q1_offset], ldq1, (ftnlen)
		3);
#line 239 "MB04WR.f"
	i__1 = *ilo - 1;
#line 239 "MB04WR.f"
	i__2 = *n - *ilo + 1;
#line 239 "MB04WR.f"
	dlaset_("All", &i__1, &i__2, &c_b9, &c_b9, &q1[*ilo * q1_dim1 + 1], 
		ldq1, (ftnlen)3);
#line 241 "MB04WR.f"
	i__1 = *ilo - 1;
#line 241 "MB04WR.f"
	dlaset_("All", n, &i__1, &c_b9, &c_b9, &q2[q2_offset], ldq2, (ftnlen)
		3);
#line 242 "MB04WR.f"
	i__1 = *ilo - 1;
#line 242 "MB04WR.f"
	i__2 = *n - *ilo + 1;
#line 242 "MB04WR.f"
	dlaset_("All", &i__1, &i__2, &c_b9, &c_b9, &q2[*ilo * q2_dim1 + 1], 
		ldq2, (ftnlen)3);
#line 244 "MB04WR.f"
	nh = *n - *ilo + 1;
#line 245 "MB04WR.f"
    }
#line 246 "MB04WR.f"
    if (compu && ! ltran) {

/*        Generate U1 and U2. */

#line 250 "MB04WR.f"
	if (nh > 0) {
#line 251 "MB04WR.f"
	    mb04wd_("No Transpose", "No Transpose", &nh, &nh, &nh, &q1[*ilo + 
		    *ilo * q1_dim1], ldq1, &q2[*ilo + *ilo * q2_dim1], ldq2, &
		    cs[*ilo], &tau[*ilo], &dwork[1], ldwork, &ierr, (ftnlen)
		    12, (ftnlen)12);
#line 254 "MB04WR.f"
	}
#line 255 "MB04WR.f"
    } else if (compu && ltran) {

/*        Generate U1**T and U2. */

#line 259 "MB04WR.f"
	if (nh > 0) {
#line 260 "MB04WR.f"
	    mb04wd_("Transpose", "No Transpose", &nh, &nh, &nh, &q1[*ilo + *
		    ilo * q1_dim1], ldq1, &q2[*ilo + *ilo * q2_dim1], ldq2, &
		    cs[*ilo], &tau[*ilo], &dwork[1], ldwork, &ierr, (ftnlen)9,
		     (ftnlen)12);
#line 263 "MB04WR.f"
	}
#line 264 "MB04WR.f"
    } else if (! compu && ! ltran) {

/*        Generate V1**T and V2**T. */

/*        Shift the vectors which define the elementary reflectors one */
/*        column to the bottom, and set the first ilo rows and */
/*        columns to those of the unit matrix. */

#line 272 "MB04WR.f"
	i__1 = *n;
#line 272 "MB04WR.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 273 "MB04WR.f"
	    i__2 = max(i__,*ilo) + 1;
#line 273 "MB04WR.f"
	    for (j = *n; j >= i__2; --j) {
#line 274 "MB04WR.f"
		q1[j + i__ * q1_dim1] = 0.;
#line 275 "MB04WR.f"
/* L10: */
#line 275 "MB04WR.f"
	    }
#line 276 "MB04WR.f"
	    i__2 = *ilo + 1;
#line 276 "MB04WR.f"
	    for (j = max(i__,*ilo); j >= i__2; --j) {
#line 277 "MB04WR.f"
		q1[j + i__ * q1_dim1] = q1[j - 1 + i__ * q1_dim1];
#line 278 "MB04WR.f"
/* L20: */
#line 278 "MB04WR.f"
	    }
#line 279 "MB04WR.f"
	    for (j = *ilo; j >= 1; --j) {
#line 280 "MB04WR.f"
		q1[j + i__ * q1_dim1] = 0.;
#line 281 "MB04WR.f"
/* L30: */
#line 281 "MB04WR.f"
	    }
#line 282 "MB04WR.f"
	    if (i__ <= *ilo) {
#line 282 "MB04WR.f"
		q1[i__ + i__ * q1_dim1] = 1.;
#line 282 "MB04WR.f"
	    }
#line 283 "MB04WR.f"
/* L40: */
#line 283 "MB04WR.f"
	}
#line 284 "MB04WR.f"
	i__1 = *n;
#line 284 "MB04WR.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 285 "MB04WR.f"
	    i__2 = max(i__,*ilo) + 1;
#line 285 "MB04WR.f"
	    for (j = *n; j >= i__2; --j) {
#line 286 "MB04WR.f"
		q2[j + i__ * q2_dim1] = 0.;
#line 287 "MB04WR.f"
/* L50: */
#line 287 "MB04WR.f"
	    }
#line 288 "MB04WR.f"
	    i__2 = *ilo + 1;
#line 288 "MB04WR.f"
	    for (j = max(i__,*ilo); j >= i__2; --j) {
#line 289 "MB04WR.f"
		q2[j + i__ * q2_dim1] = q2[j - 1 + i__ * q2_dim1];
#line 290 "MB04WR.f"
/* L60: */
#line 290 "MB04WR.f"
	    }
#line 291 "MB04WR.f"
	    for (j = *ilo; j >= 1; --j) {
#line 292 "MB04WR.f"
		q2[j + i__ * q2_dim1] = 0.;
#line 293 "MB04WR.f"
/* L70: */
#line 293 "MB04WR.f"
	    }
#line 294 "MB04WR.f"
/* L80: */
#line 294 "MB04WR.f"
	}

#line 296 "MB04WR.f"
	nh = *n - *ilo;
#line 297 "MB04WR.f"
	if (nh > 0) {
#line 298 "MB04WR.f"
	    mb04wd_("Transpose", "Transpose", &nh, &nh, &nh, &q1[*ilo + 1 + (*
		    ilo + 1) * q1_dim1], ldq1, &q2[*ilo + 1 + (*ilo + 1) * 
		    q2_dim1], ldq2, &cs[*ilo], &tau[*ilo], &dwork[1], ldwork, 
		    &ierr, (ftnlen)9, (ftnlen)9);
#line 301 "MB04WR.f"
	}
#line 302 "MB04WR.f"
    } else if (! compu && ltran) {

/*        Generate V1 and V2**T. */

/*        Shift the vectors which define the elementary reflectors one */
/*        column to the right/bottom, and set the first ilo rows and */
/*        columns to those of the unit matrix. */

#line 310 "MB04WR.f"
	i__1 = *ilo + 1;
#line 310 "MB04WR.f"
	for (j = *n; j >= i__1; --j) {
#line 311 "MB04WR.f"
	    i__2 = j - 1;
#line 311 "MB04WR.f"
	    for (i__ = 1; i__ <= i__2; ++i__) {
#line 312 "MB04WR.f"
		q1[i__ + j * q1_dim1] = 0.;
#line 313 "MB04WR.f"
/* L90: */
#line 313 "MB04WR.f"
	    }
#line 314 "MB04WR.f"
	    i__2 = *n;
#line 314 "MB04WR.f"
	    for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 315 "MB04WR.f"
		q1[i__ + j * q1_dim1] = q1[i__ + (j - 1) * q1_dim1];
#line 316 "MB04WR.f"
/* L100: */
#line 316 "MB04WR.f"
	    }
#line 317 "MB04WR.f"
/* L110: */
#line 317 "MB04WR.f"
	}
#line 318 "MB04WR.f"
	dlaset_("All", n, ilo, &c_b9, &c_b10, &q1[q1_offset], ldq1, (ftnlen)3)
		;
#line 319 "MB04WR.f"
	i__1 = *n;
#line 319 "MB04WR.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 320 "MB04WR.f"
	    i__2 = max(i__,*ilo) + 1;
#line 320 "MB04WR.f"
	    for (j = *n; j >= i__2; --j) {
#line 321 "MB04WR.f"
		q2[j + i__ * q2_dim1] = 0.;
#line 322 "MB04WR.f"
/* L120: */
#line 322 "MB04WR.f"
	    }
#line 323 "MB04WR.f"
	    i__2 = *ilo + 1;
#line 323 "MB04WR.f"
	    for (j = max(i__,*ilo); j >= i__2; --j) {
#line 324 "MB04WR.f"
		q2[j + i__ * q2_dim1] = q2[j - 1 + i__ * q2_dim1];
#line 325 "MB04WR.f"
/* L130: */
#line 325 "MB04WR.f"
	    }
#line 326 "MB04WR.f"
	    for (j = *ilo; j >= 1; --j) {
#line 327 "MB04WR.f"
		q2[j + i__ * q2_dim1] = 0.;
#line 328 "MB04WR.f"
/* L140: */
#line 328 "MB04WR.f"
	    }
#line 329 "MB04WR.f"
/* L150: */
#line 329 "MB04WR.f"
	}
#line 330 "MB04WR.f"
	nh = *n - *ilo;

#line 332 "MB04WR.f"
	if (nh > 0) {
#line 333 "MB04WR.f"
	    mb04wd_("No Transpose", "Transpose", &nh, &nh, &nh, &q1[*ilo + 1 
		    + (*ilo + 1) * q1_dim1], ldq1, &q2[*ilo + 1 + (*ilo + 1) *
		     q2_dim1], ldq2, &cs[*ilo], &tau[*ilo], &dwork[1], ldwork,
		     &ierr, (ftnlen)12, (ftnlen)9);
#line 336 "MB04WR.f"
	}
#line 337 "MB04WR.f"
    }
#line 338 "MB04WR.f"
    return 0;
/* *** Last line of MB04WR *** */
} /* mb04wr_ */

