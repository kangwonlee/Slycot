#line 1 "MB04WP.f"
/* MB04WP.f -- translated by f2c (version 20100827).
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

#line 1 "MB04WP.f"
/* Table of constant values */

static doublereal c_b7 = 0.;
static doublereal c_b8 = 1.;

/* Subroutine */ int mb04wp_(integer *n, integer *ilo, doublereal *u1, 
	integer *ldu1, doublereal *u2, integer *ldu2, doublereal *cs, 
	doublereal *tau, doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer u1_dim1, u1_offset, u2_dim1, u2_offset, i__1, i__2;

    /* Local variables */
    static integer i__, j, nh, ierr;
    extern /* Subroutine */ int mb04wd_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dlaset_(char *, integer *, integer *, doublereal 
	    *, doublereal *, doublereal *, integer *, ftnlen), xerbla_(char *,
	     integer *, ftnlen);


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

/*     To generate an orthogonal symplectic matrix U, which is defined as */
/*     a product of symplectic reflectors and Givens rotators */

/*     U = diag( H(1),H(1) )      G(1)  diag( F(1),F(1) ) */
/*         diag( H(2),H(2) )      G(2)  diag( F(2),F(2) ) */
/*                                .... */
/*         diag( H(n-1),H(n-1) ) G(n-1) diag( F(n-1),F(n-1) ). */

/*     as returned by MB04PU. The matrix U is returned in terms of its */
/*     first N rows */

/*                      [  U1   U2 ] */
/*                  U = [          ]. */
/*                      [ -U2   U1 ] */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices U1 and U2.  N >= 0. */

/*     ILO     (input) INTEGER */
/*             ILO must have the same value as in the previous call of */
/*             MB04PU. U is equal to the unit matrix except in the */
/*             submatrix */
/*             U([ilo+1:n n+ilo+1:2*n], [ilo+1:n n+ilo+1:2*n]). */
/*             1 <= ILO <= N, if N > 0; ILO = 1, if N = 0. */

/*     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain in its i-th column the vector which defines the */
/*             elementary reflector F(i). */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix U1. */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1.  LDU1 >= MAX(1,N). */

/*     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain in its i-th column the vector which defines the */
/*             elementary reflector H(i) and, on the subdiagonal, the */
/*             scalar factor of H(i). */
/*             On exit, the leading N-by-N part of this array contains */
/*             the matrix U2. */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2.  LDU2 >= MAX(1,N). */

/*     CS      (input) DOUBLE PRECISION array, dimension (2N-2) */
/*             On entry, the first 2N-2 elements of this array must */
/*             contain the cosines and sines of the symplectic Givens */
/*             rotators G(i). */

/*     TAU     (input) DOUBLE PRECISION array, dimension (N-1) */
/*             On entry, the first N-1 elements of this array must */
/*             contain the scalar factors of the elementary reflectors */
/*             F(i). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -10,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. LDWORK >= MAX(1,2*(N-ILO)). */
/*             For optimum performance LDWORK should be larger. (See */
/*             SLICOT Library routine MB04WD). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires O(N**3) floating point operations and is */
/*     strongly backward stable. */

/*     REFERENCES */

/*     [1] C. F. VAN LOAN: */
/*         A symplectic method for approximating all the eigenvalues of */
/*         a Hamiltonian matrix. */
/*         Linear Algebra and its Applications, 61, pp. 233-251, 1984. */

/*     [2] D. KRESSNER: */
/*         Block algorithms for orthogonal symplectic factorizations. */
/*         BIT, 43 (4), pp. 775-790, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner (Technical Univ. Berlin, Germany) and */
/*     P. Benner (Technical Univ. Chemnitz, Germany), December 2003. */

/*     REVISIONS */

/*     V. Sima, Nov. 2008 (SLICOT version of the HAPACK routine DOSGPV). */

/*     KEYWORDS */

/*     Elementary matrix operations, orthogonal symplectic matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Check the scalar input parameters. */

#line 153 "MB04WP.f"
    /* Parameter adjustments */
#line 153 "MB04WP.f"
    u1_dim1 = *ldu1;
#line 153 "MB04WP.f"
    u1_offset = 1 + u1_dim1;
#line 153 "MB04WP.f"
    u1 -= u1_offset;
#line 153 "MB04WP.f"
    u2_dim1 = *ldu2;
#line 153 "MB04WP.f"
    u2_offset = 1 + u2_dim1;
#line 153 "MB04WP.f"
    u2 -= u2_offset;
#line 153 "MB04WP.f"
    --cs;
#line 153 "MB04WP.f"
    --tau;
#line 153 "MB04WP.f"
    --dwork;
#line 153 "MB04WP.f"

#line 153 "MB04WP.f"
    /* Function Body */
#line 153 "MB04WP.f"
    *info = 0;
#line 154 "MB04WP.f"
    if (*n < 0) {
#line 155 "MB04WP.f"
	*info = -1;
#line 156 "MB04WP.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 157 "MB04WP.f"
	*info = -2;
#line 158 "MB04WP.f"
    } else if (*ldu1 < max(1,*n)) {
#line 159 "MB04WP.f"
	*info = -4;
#line 160 "MB04WP.f"
    } else if (*ldu2 < max(1,*n)) {
#line 161 "MB04WP.f"
	*info = -6;
#line 162 "MB04WP.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
#line 162 "MB04WP.f"
	i__1 = 1, i__2 = *n - *ilo << 1;
#line 162 "MB04WP.f"
	if (*ldwork < max(i__1,i__2)) {
/* Computing MAX */
#line 163 "MB04WP.f"
	    i__1 = 1, i__2 = *n - *ilo << 1;
#line 163 "MB04WP.f"
	    dwork[1] = (doublereal) max(i__1,i__2);
#line 164 "MB04WP.f"
	    *info = -10;
#line 165 "MB04WP.f"
	}
#line 165 "MB04WP.f"
    }

/*     Return if there were illegal values. */

#line 169 "MB04WP.f"
    if (*info != 0) {
#line 170 "MB04WP.f"
	i__1 = -(*info);
#line 170 "MB04WP.f"
	xerbla_("MB04WP", &i__1, (ftnlen)6);
#line 171 "MB04WP.f"
	return 0;
#line 172 "MB04WP.f"
    }

/*     Quick return if possible. */

#line 176 "MB04WP.f"
    if (*n == 0) {
#line 177 "MB04WP.f"
	dwork[1] = 1.;
#line 178 "MB04WP.f"
	return 0;
#line 179 "MB04WP.f"
    }

/*     Shift the vectors which define the elementary reflectors one */
/*     column to the right, and set the first ilo rows and columns to */
/*     those of the unit matrix. */

#line 185 "MB04WP.f"
    i__1 = *ilo + 1;
#line 185 "MB04WP.f"
    for (j = *n; j >= i__1; --j) {
#line 186 "MB04WP.f"
	i__2 = j - 1;
#line 186 "MB04WP.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 187 "MB04WP.f"
	    u1[i__ + j * u1_dim1] = 0.;
#line 188 "MB04WP.f"
/* L10: */
#line 188 "MB04WP.f"
	}
#line 189 "MB04WP.f"
	i__2 = *n;
#line 189 "MB04WP.f"
	for (i__ = j + 1; i__ <= i__2; ++i__) {
#line 190 "MB04WP.f"
	    u1[i__ + j * u1_dim1] = u1[i__ + (j - 1) * u1_dim1];
#line 191 "MB04WP.f"
/* L20: */
#line 191 "MB04WP.f"
	}
#line 192 "MB04WP.f"
/* L30: */
#line 192 "MB04WP.f"
    }
#line 193 "MB04WP.f"
    dlaset_("All", n, ilo, &c_b7, &c_b8, &u1[u1_offset], ldu1, (ftnlen)3);
#line 194 "MB04WP.f"
    i__1 = *ilo + 1;
#line 194 "MB04WP.f"
    for (j = *n; j >= i__1; --j) {
#line 195 "MB04WP.f"
	i__2 = j - 1;
#line 195 "MB04WP.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 196 "MB04WP.f"
	    u2[i__ + j * u2_dim1] = 0.;
#line 197 "MB04WP.f"
/* L40: */
#line 197 "MB04WP.f"
	}
#line 198 "MB04WP.f"
	i__2 = *n;
#line 198 "MB04WP.f"
	for (i__ = j; i__ <= i__2; ++i__) {
#line 199 "MB04WP.f"
	    u2[i__ + j * u2_dim1] = u2[i__ + (j - 1) * u2_dim1];
#line 200 "MB04WP.f"
/* L50: */
#line 200 "MB04WP.f"
	}
#line 201 "MB04WP.f"
/* L60: */
#line 201 "MB04WP.f"
    }
#line 202 "MB04WP.f"
    dlaset_("All", n, ilo, &c_b7, &c_b7, &u2[u2_offset], ldu2, (ftnlen)3);
#line 203 "MB04WP.f"
    nh = *n - *ilo;
#line 204 "MB04WP.f"
    if (nh > 0) {
#line 205 "MB04WP.f"
	mb04wd_("No Transpose", "No Transpose", &nh, &nh, &nh, &u1[*ilo + 1 + 
		(*ilo + 1) * u1_dim1], ldu1, &u2[*ilo + 1 + (*ilo + 1) * 
		u2_dim1], ldu2, &cs[*ilo], &tau[*ilo], &dwork[1], ldwork, &
		ierr, (ftnlen)12, (ftnlen)12);
#line 208 "MB04WP.f"
    }
#line 209 "MB04WP.f"
    return 0;
/* *** Last line of MB04WP *** */
} /* mb04wp_ */

