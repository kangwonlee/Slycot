#line 1 "TF01RD.f"
/* TF01RD.f -- translated by f2c (version 20100827).
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

#line 1 "TF01RD.f"
/* Table of constant values */

static doublereal c_b8 = 1.;
static doublereal c_b9 = 0.;

/* Subroutine */ int tf01rd_(integer *na, integer *nb, integer *nc, integer *
	n, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *h__, integer *ldh, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, h_dim1, 
	    h_offset, i__1, i__2;

    /* Local variables */
    static integer i__, k, ldw;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer jwork;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
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

/*     To compute N Markov parameters M(1), M(2),..., M(N) from the */
/*     parameters (A,B,C) of a linear time-invariant system, where each */
/*     M(k) is an NC-by-NB matrix and k = 1,2,...,N. */

/*     All matrices are treated as dense, and hence TF01RD is not */
/*     intended for large sparse problems. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NA      (input) INTEGER */
/*             The order of the matrix A.  NA >= 0. */

/*     NB      (input) INTEGER */
/*             The number of system inputs.  NB >= 0. */

/*     NC      (input) INTEGER */
/*             The number of system outputs.  NC >= 0. */

/*     N       (input) INTEGER */
/*             The number of Markov parameters M(k) to be computed. */
/*             N >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,NA) */
/*             The leading NA-by-NA part of this array must contain the */
/*             state matrix A of the system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,NA). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,NB) */
/*             The leading NA-by-NB part of this array must contain the */
/*             input matrix B of the system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,NA). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,NA) */
/*             The leading NC-by-NA part of this array must contain the */
/*             output matrix C of the system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,NC). */

/*     H       (output) DOUBLE PRECISION array, dimension (LDH,N*NB) */
/*             The leading NC-by-N*NB part of this array contains the */
/*             multivariable parameters M(k), where each parameter M(k) */
/*             is an NC-by-NB matrix and k = 1,2,...,N. The Markov */
/*             parameters are stored such that H(i,(k-1)xNB+j) contains */
/*             the (i,j)-th element of M(k) for i = 1,2,...,NC and */
/*             j = 1,2,...,NB. */

/*     LDH     INTEGER */
/*             The leading dimension of array H.  LDH >= MAX(1,NC). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, 2*NA*NC). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     For the linear time-invariant discrete-time system */

/*            x(k+1) = A x(k) + B u(k) */
/*             y(k)  = C x(k) + D u(k), */

/*     the transfer function matrix G(z) is given by */
/*                            -1 */
/*              G(z) = C(zI-A)  B + D */
/*                             -1        -2     2   -3 */
/*                   = D + CB z   + CAB z   + CA B z   + ...          (1) */

/*     Using Markov parameters, G(z) can also be written as */
/*                                 -1        -2        -3 */
/*              G(z) = M(0) + M(1)z   + M(2)z   + M(3)z   + ...       (2) */

/*                                                               k-1 */
/*     Equating (1) and (2), we find that M(0) = D and M(k) = C A    B */
/*     for k > 0, from which the Markov parameters M(1),M(2)...,M(N) are */
/*     computed. */

/*     REFERENCES */

/*     [1] Chen, C.T. */
/*         Introduction to Linear System Theory. */
/*         H.R.W. Series in Electrical Engineering, Electronics and */
/*         Systems, Holt, Rinehart and Winston Inc., London, 1970. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately (NA + NB) x NA x NC x N */
/*     multiplications and additions. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TF01FD by S. Van Huffel, Katholieke */
/*     Univ. Leuven, Belgium. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Markov parameters, multivariable system, time-invariant system, */
/*     transfer function, transfer matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 161 "TF01RD.f"
    /* Parameter adjustments */
#line 161 "TF01RD.f"
    a_dim1 = *lda;
#line 161 "TF01RD.f"
    a_offset = 1 + a_dim1;
#line 161 "TF01RD.f"
    a -= a_offset;
#line 161 "TF01RD.f"
    b_dim1 = *ldb;
#line 161 "TF01RD.f"
    b_offset = 1 + b_dim1;
#line 161 "TF01RD.f"
    b -= b_offset;
#line 161 "TF01RD.f"
    c_dim1 = *ldc;
#line 161 "TF01RD.f"
    c_offset = 1 + c_dim1;
#line 161 "TF01RD.f"
    c__ -= c_offset;
#line 161 "TF01RD.f"
    h_dim1 = *ldh;
#line 161 "TF01RD.f"
    h_offset = 1 + h_dim1;
#line 161 "TF01RD.f"
    h__ -= h_offset;
#line 161 "TF01RD.f"
    --dwork;
#line 161 "TF01RD.f"

#line 161 "TF01RD.f"
    /* Function Body */
#line 161 "TF01RD.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 165 "TF01RD.f"
    if (*na < 0) {
#line 166 "TF01RD.f"
	*info = -1;
#line 167 "TF01RD.f"
    } else if (*nb < 0) {
#line 168 "TF01RD.f"
	*info = -2;
#line 169 "TF01RD.f"
    } else if (*nc < 0) {
#line 170 "TF01RD.f"
	*info = -3;
#line 171 "TF01RD.f"
    } else if (*n < 0) {
#line 172 "TF01RD.f"
	*info = -4;
#line 173 "TF01RD.f"
    } else if (*lda < max(1,*na)) {
#line 174 "TF01RD.f"
	*info = -6;
#line 175 "TF01RD.f"
    } else if (*ldb < max(1,*na)) {
#line 176 "TF01RD.f"
	*info = -8;
#line 177 "TF01RD.f"
    } else if (*ldc < max(1,*nc)) {
#line 178 "TF01RD.f"
	*info = -10;
#line 179 "TF01RD.f"
    } else if (*ldh < max(1,*nc)) {
#line 180 "TF01RD.f"
	*info = -12;
#line 181 "TF01RD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 181 "TF01RD.f"
	i__1 = 1, i__2 = (*na << 1) * *nc;
#line 181 "TF01RD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 182 "TF01RD.f"
	    *info = -14;
#line 183 "TF01RD.f"
	}
#line 183 "TF01RD.f"
    }

#line 185 "TF01RD.f"
    if (*info != 0) {

/*        Error return. */

#line 189 "TF01RD.f"
	i__1 = -(*info);
#line 189 "TF01RD.f"
	xerbla_("TF01RD", &i__1, (ftnlen)6);
#line 190 "TF01RD.f"
	return 0;
#line 191 "TF01RD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 195 "TF01RD.f"
    i__1 = min(*na,*nb), i__1 = min(i__1,*nc);
#line 195 "TF01RD.f"
    if (min(i__1,*n) == 0) {
#line 195 "TF01RD.f"
	return 0;
#line 195 "TF01RD.f"
    }

#line 198 "TF01RD.f"
    jwork = *nc * *na + 1;
#line 199 "TF01RD.f"
    ldw = max(1,*nc);
#line 200 "TF01RD.f"
    i__ = 1;

/*     Copy C in the workspace beginning from the position JWORK. */
/*     This workspace will contain the product C*A**(K-1), K = 1,2,...,N. */

#line 205 "TF01RD.f"
    dlacpy_("Full", nc, na, &c__[c_offset], ldc, &dwork[jwork], &ldw, (ftnlen)
	    4);

/*     Form M(1), M(2), ..., M(N). */

#line 209 "TF01RD.f"
    i__1 = *n;
#line 209 "TF01RD.f"
    for (k = 1; k <= i__1; ++k) {
#line 210 "TF01RD.f"
	dlacpy_("Full", nc, na, &dwork[jwork], &ldw, &dwork[1], &ldw, (ftnlen)
		4);

/*        Form (C * A**(K-1)) * B = M(K). */

#line 214 "TF01RD.f"
	dgemm_("No transpose", "No transpose", nc, nb, na, &c_b8, &dwork[1], &
		ldw, &b[b_offset], ldb, &c_b9, &h__[i__ * h_dim1 + 1], ldh, (
		ftnlen)12, (ftnlen)12);

#line 217 "TF01RD.f"
	if (k != *n) {

/*           Form C * A**K. */

#line 221 "TF01RD.f"
	    dgemm_("No transpose", "No transpose", nc, na, na, &c_b8, &dwork[
		    1], &ldw, &a[a_offset], lda, &c_b9, &dwork[jwork], &ldw, (
		    ftnlen)12, (ftnlen)12);

#line 224 "TF01RD.f"
	    i__ += *nb;
#line 225 "TF01RD.f"
	}
#line 226 "TF01RD.f"
/* L10: */
#line 226 "TF01RD.f"
    }

#line 228 "TF01RD.f"
    return 0;
/* *** Last line of TF01RD *** */
} /* tf01rd_ */

