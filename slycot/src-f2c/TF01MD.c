#line 1 "TF01MD.f"
/* TF01MD.f -- translated by f2c (version 20100827).
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

#line 1 "TF01MD.f"
/* Table of constant values */

static doublereal c_b4 = 0.;
static doublereal c_b8 = 1.;
static integer c__1 = 1;

/* Subroutine */ int tf01md_(integer *n, integer *m, integer *p, integer *ny, 
	doublereal *a, integer *lda, doublereal *b, integer *ldb, doublereal *
	c__, integer *ldc, doublereal *d__, integer *ldd, doublereal *u, 
	integer *ldu, doublereal *x, doublereal *y, integer *ldy, doublereal *
	dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, u_dim1, u_offset, y_dim1, y_offset, i__1;

    /* Local variables */
    static integer ik;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaset_(char *, integer *, integer *, 
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

/*     To compute the output sequence of a linear time-invariant */
/*     open-loop system given by its discrete-time state-space model */
/*     (A,B,C,D), where A is an N-by-N general matrix. */

/*     The initial state vector x(1) must be supplied by the user. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NY      (input) INTEGER */
/*             The number of output vectors y(k) to be computed. */
/*             NY >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state matrix A of the system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input matrix B of the system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             output matrix C of the system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             direct link matrix D of the system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,NY) */
/*             The leading M-by-NY part of this array must contain the */
/*             input vector sequence u(k), for k = 1,2,...,NY. */
/*             Specifically, the k-th column of U must contain u(k). */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,M). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the initial state vector */
/*             x(1) which consists of the N initial states of the system. */
/*             On exit, this array contains the final state vector */
/*             x(NY+1) of the N states of the system at instant NY. */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,NY) */
/*             The leading P-by-NY part of this array contains the output */
/*             vector sequence y(1),y(2),...,y(NY) such that the k-th */
/*             column of Y contains y(k) (the outputs at instant k), */
/*             for k = 1,2,...,NY. */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= MAX(1,P). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given an initial state vector x(1), the output vector sequence */
/*     y(1), y(2),..., y(NY) is obtained via the formulae */

/*        x(k+1) = A x(k) + B u(k) */
/*        y(k)   = C x(k) + D u(k), */

/*     where each element y(k) is a vector of length P containing the */
/*     outputs at instant k and k = 1,2,...,NY. */

/*     REFERENCES */

/*     [1] Luenberger, D.G. */
/*         Introduction to Dynamic Systems: Theory, Models and */
/*         Applications. */
/*         John Wiley & Sons, New York, 1979. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately (N + M) x (N + P) x NY */
/*     multiplications and additions. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Dec. 1996. */
/*     Supersedes Release 2.0 routine TF01AD by S. Van Huffel, Katholieke */
/*     Univ. Leuven, Belgium. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2003. */

/*     KEYWORDS */

/*     Discrete-time system, multivariable system, state-space model, */
/*     state-space representation, time response. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 165 "TF01MD.f"
    /* Parameter adjustments */
#line 165 "TF01MD.f"
    a_dim1 = *lda;
#line 165 "TF01MD.f"
    a_offset = 1 + a_dim1;
#line 165 "TF01MD.f"
    a -= a_offset;
#line 165 "TF01MD.f"
    b_dim1 = *ldb;
#line 165 "TF01MD.f"
    b_offset = 1 + b_dim1;
#line 165 "TF01MD.f"
    b -= b_offset;
#line 165 "TF01MD.f"
    c_dim1 = *ldc;
#line 165 "TF01MD.f"
    c_offset = 1 + c_dim1;
#line 165 "TF01MD.f"
    c__ -= c_offset;
#line 165 "TF01MD.f"
    d_dim1 = *ldd;
#line 165 "TF01MD.f"
    d_offset = 1 + d_dim1;
#line 165 "TF01MD.f"
    d__ -= d_offset;
#line 165 "TF01MD.f"
    u_dim1 = *ldu;
#line 165 "TF01MD.f"
    u_offset = 1 + u_dim1;
#line 165 "TF01MD.f"
    u -= u_offset;
#line 165 "TF01MD.f"
    --x;
#line 165 "TF01MD.f"
    y_dim1 = *ldy;
#line 165 "TF01MD.f"
    y_offset = 1 + y_dim1;
#line 165 "TF01MD.f"
    y -= y_offset;
#line 165 "TF01MD.f"
    --dwork;
#line 165 "TF01MD.f"

#line 165 "TF01MD.f"
    /* Function Body */
#line 165 "TF01MD.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 169 "TF01MD.f"
    if (*n < 0) {
#line 170 "TF01MD.f"
	*info = -1;
#line 171 "TF01MD.f"
    } else if (*m < 0) {
#line 172 "TF01MD.f"
	*info = -2;
#line 173 "TF01MD.f"
    } else if (*p < 0) {
#line 174 "TF01MD.f"
	*info = -3;
#line 175 "TF01MD.f"
    } else if (*ny < 0) {
#line 176 "TF01MD.f"
	*info = -4;
#line 177 "TF01MD.f"
    } else if (*lda < max(1,*n)) {
#line 178 "TF01MD.f"
	*info = -6;
#line 179 "TF01MD.f"
    } else if (*ldb < max(1,*n)) {
#line 180 "TF01MD.f"
	*info = -8;
#line 181 "TF01MD.f"
    } else if (*ldc < max(1,*p)) {
#line 182 "TF01MD.f"
	*info = -10;
#line 183 "TF01MD.f"
    } else if (*ldd < max(1,*p)) {
#line 184 "TF01MD.f"
	*info = -12;
#line 185 "TF01MD.f"
    } else if (*ldu < max(1,*m)) {
#line 186 "TF01MD.f"
	*info = -14;
#line 187 "TF01MD.f"
    } else if (*ldy < max(1,*p)) {
#line 188 "TF01MD.f"
	*info = -17;
#line 189 "TF01MD.f"
    }

#line 191 "TF01MD.f"
    if (*info != 0) {

/*        Error return. */

#line 195 "TF01MD.f"
	i__1 = -(*info);
#line 195 "TF01MD.f"
	xerbla_("TF01MD", &i__1, (ftnlen)6);
#line 196 "TF01MD.f"
	return 0;
#line 197 "TF01MD.f"
    }

/*     Quick return if possible. */

#line 201 "TF01MD.f"
    if (min(*p,*ny) == 0) {
#line 202 "TF01MD.f"
	return 0;
#line 203 "TF01MD.f"
    } else if (*n == 0) {

/*        Non-dynamic system: compute the output vectors. */

#line 207 "TF01MD.f"
	if (*m == 0) {
#line 208 "TF01MD.f"
	    dlaset_("Full", p, ny, &c_b4, &c_b4, &y[y_offset], ldy, (ftnlen)4)
		    ;
#line 209 "TF01MD.f"
	} else {
#line 210 "TF01MD.f"
	    dgemm_("No transpose", "No transpose", p, ny, m, &c_b8, &d__[
		    d_offset], ldd, &u[u_offset], ldu, &c_b4, &y[y_offset], 
		    ldy, (ftnlen)12, (ftnlen)12);
#line 212 "TF01MD.f"
	}
#line 213 "TF01MD.f"
	return 0;
#line 214 "TF01MD.f"
    }

#line 216 "TF01MD.f"
    i__1 = *ny;
#line 216 "TF01MD.f"
    for (ik = 1; ik <= i__1; ++ik) {
#line 217 "TF01MD.f"
	dgemv_("No transpose", p, n, &c_b8, &c__[c_offset], ldc, &x[1], &c__1,
		 &c_b4, &y[ik * y_dim1 + 1], &c__1, (ftnlen)12);

#line 220 "TF01MD.f"
	dgemv_("No transpose", n, n, &c_b8, &a[a_offset], lda, &x[1], &c__1, &
		c_b4, &dwork[1], &c__1, (ftnlen)12);
#line 222 "TF01MD.f"
	dgemv_("No transpose", n, m, &c_b8, &b[b_offset], ldb, &u[ik * u_dim1 
		+ 1], &c__1, &c_b8, &dwork[1], &c__1, (ftnlen)12);

#line 225 "TF01MD.f"
	dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);
#line 226 "TF01MD.f"
/* L10: */
#line 226 "TF01MD.f"
    }

#line 228 "TF01MD.f"
    dgemm_("No transpose", "No transpose", p, ny, m, &c_b8, &d__[d_offset], 
	    ldd, &u[u_offset], ldu, &c_b8, &y[y_offset], ldy, (ftnlen)12, (
	    ftnlen)12);

#line 231 "TF01MD.f"
    return 0;
/* *** Last line of TF01MD *** */
} /* tf01md_ */
