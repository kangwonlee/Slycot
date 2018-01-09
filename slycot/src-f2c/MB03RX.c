#line 1 "MB03RX.f"
/* MB03RX.f -- translated by f2c (version 20100827).
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

#line 1 "MB03RX.f"
/* Subroutine */ int mb03rx_(char *jobv, integer *n, integer *kl, integer *ku,
	 doublereal *a, integer *lda, doublereal *x, integer *ldx, doublereal 
	*wr, doublereal *wi, doublereal *dwork, ftnlen jobv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, x_dim1, x_offset;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer l, ierr, ifst, ilst;
    extern /* Subroutine */ int dtrexc_(char *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen);


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

/*     To reorder the diagonal blocks of the principal submatrix between */
/*     the indices KL and KU (KU >= KL) of a real Schur form matrix A */
/*     together with their eigenvalues, using orthogonal similarity */
/*     transformations, such that the block specified by KU is moved in */
/*     the position KL. The transformations are optionally postmultiplied */
/*     in a given matrix X. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBV    CHARACTER*1 */
/*             Specifies whether or not the transformations are */
/*             accumulated, as follows: */
/*             = 'N':  The transformations are not accumulated; */
/*             = 'V':  The transformations are accumulated in X (the */
/*                     given matrix X is updated). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrices A and X.  N >= 0. */

/*     KL      (input) INTEGER */
/*             The lower boundary index for the rows and columns of the */
/*             principal submatrix of A whose diagonal blocks are to be */
/*             reordered, and also the target position for the block to */
/*             be moved.  1 <= KL <= KU <= N. */

/*     KU      (input/output) INTEGER */
/*             On entry, KU specifies the upper boundary index for the */
/*             rows and columns of the principal submatrix of A whose */
/*             diagonal blocks are to be reordered, and also the original */
/*             position for the block to be moved.  1 <= KL <= KU <= N. */
/*             On exit, KU specifies the upper boundary index for the */
/*             rows and columns of the principal submatrix of A whose */
/*             diagonal blocks have been reordered. The given value will */
/*             be increased by 1 if the moved block was 2-by-2 and it has */
/*             been replaced by two 1-by-1 blocks. Otherwise, its input */
/*             value is preserved. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the matrix A in real Schur canonical form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the ordered real Schur canonical form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (LDX,N) */
/*             On entry, if JOBV = 'V', the leading N-by-N part of this */
/*             array must contain a given matrix X. */
/*             On exit, if JOBV = 'V', the leading N-by-N part of this */
/*             array contains the product of the given matrix X and the */
/*             transformation matrix that performed the reordering of A. */
/*             If JOBV = 'N', this array is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of array X. */
/*             LDX >= 1,        if JOBV = 'N'; */
/*             LDX >= MAX(1,N), if JOBV = 'V'. */

/*     WR,     (input/output) DOUBLE PRECISION arrays, dimension (N) */
/*     WI      On entry, these arrays must contain the real and imaginary */
/*             parts, respectively, of the eigenvalues of the matrix A. */
/*             On exit, these arrays contain the real and imaginary */
/*             parts, respectively, of the eigenvalues of the matrix A, */
/*             possibly reordered. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     METHOD */

/*     An attempt is made to move the block in the position (KU,KU) to */
/*     the position (KL,KL) by a sequence of orthogonal similarity */
/*     transformations, each swapping two consecutive blocks. The */
/*     standard algorithm [1], [2] usually succeeds to perform this */
/*     reordering. A failure of this algorithm means that two consecutive */
/*     blocks (one of them being the desired block possibly moved) are */
/*     too close to swap. In such a case, the leading block of the two */
/*     is tried to be moved in the position (KL,KL) and the procedure is */
/*     repeated. */

/*     REFERENCES */

/*     [1] Stewart, G.W. */
/*         HQR3 and EXCHQZ: FORTRAN subroutines for calculating and */
/*         ordering the eigenvalues of a real upper Hessenberg matrix. */
/*         ACM TOMS, 2, pp. 275-280, 1976. */

/*     [2] Anderson, E., Bai, Z., Bischof, C., Demmel, J., Dongarra, J., */
/*         Du Croz, J., Greenbaum, A., Hammarling, S., McKenney, A., */
/*         Ostrouchov, S., and Sorensen, D. */
/*         LAPACK Users' Guide: Second Edition. */
/*         SIAM, Philadelphia, 1995. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable. If some eigenvalues are */
/*     ill-conditioned, their returned values could differ much from */
/*     their input values. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Eigenvalue, orthogonal transformation, real Schur form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 159 "MB03RX.f"
    /* Parameter adjustments */
#line 159 "MB03RX.f"
    a_dim1 = *lda;
#line 159 "MB03RX.f"
    a_offset = 1 + a_dim1;
#line 159 "MB03RX.f"
    a -= a_offset;
#line 159 "MB03RX.f"
    x_dim1 = *ldx;
#line 159 "MB03RX.f"
    x_offset = 1 + x_dim1;
#line 159 "MB03RX.f"
    x -= x_offset;
#line 159 "MB03RX.f"
    --wr;
#line 159 "MB03RX.f"
    --wi;
#line 159 "MB03RX.f"
    --dwork;
#line 159 "MB03RX.f"

#line 159 "MB03RX.f"
    /* Function Body */
#line 159 "MB03RX.f"
    if (*ku > *kl) {

/*        Try to move the block in position (KU,KU) to position (KL,KL). */

#line 163 "MB03RX.f"
	ifst = *ku;
/*        REPEAT */
#line 165 "MB03RX.f"
L10:
#line 166 "MB03RX.f"
	ilst = *kl;
#line 167 "MB03RX.f"
	dtrexc_(jobv, n, &a[a_offset], lda, &x[x_offset], ldx, &ifst, &ilst, &
		dwork[1], &ierr, (ftnlen)1);
#line 168 "MB03RX.f"
	if (ierr != 0) {

/*           During calculations, two adjacent blocks were too close */
/*           to swap; the desired block cannot be moved further, but the */
/*           block above it is suitable and is tried for moving. The */
/*           number of repeat cycles is usually 1, and at most the number */
/*           of blocks between the current position and the position KL. */

#line 176 "MB03RX.f"
	    ifst = ilst - 1;
#line 177 "MB03RX.f"
	    if (ifst > 1) {
#line 178 "MB03RX.f"
		if (a[ifst + (ifst - 1) * a_dim1] != 0.) {
#line 178 "MB03RX.f"
		    ifst = ilst - 2;
#line 178 "MB03RX.f"
		}
#line 180 "MB03RX.f"
	    }
#line 181 "MB03RX.f"
	    if (ilst > *kl) {
#line 181 "MB03RX.f"
		goto L10;
#line 181 "MB03RX.f"
	    }
#line 183 "MB03RX.f"
	}
/*        UNTIL ( ILST.EQ.KL on output from DTREXC ) */

/*        Recompute the eigenvalues for the modified part of A. */
/*        Note that KU must be incremented if the moved block was 2-by-2 */
/*        and it has been replaced by two 1-by-1 blocks. */

#line 190 "MB03RX.f"
	if (wi[*ku] != 0.) {
#line 191 "MB03RX.f"
	    if (a[*ku + 1 + *ku * a_dim1] == 0.) {
#line 191 "MB03RX.f"
		++(*ku);
#line 191 "MB03RX.f"
	    }
#line 193 "MB03RX.f"
	}

#line 195 "MB03RX.f"
	l = *kl;
/*        WHILE ( L.LT.KU .OR. ( L.EQ.KU .AND. L.LT.N ) ) DO */
#line 197 "MB03RX.f"
L20:
#line 197 "MB03RX.f"
	if (l < *ku || l == *ku && l < *n) {
#line 198 "MB03RX.f"
	    if (a[l + 1 + l * a_dim1] != 0.) {

/*              A 2x2 block. */

#line 202 "MB03RX.f"
		wr[l] = a[l + l * a_dim1];
#line 203 "MB03RX.f"
		wr[l + 1] = wr[l];
#line 204 "MB03RX.f"
		wi[l] = sqrt((d__1 = a[l + (l + 1) * a_dim1], abs(d__1))) * 
			sqrt((d__2 = a[l + 1 + l * a_dim1], abs(d__2)));
#line 206 "MB03RX.f"
		wi[l + 1] = -wi[l];
#line 207 "MB03RX.f"
		l += 2;
#line 208 "MB03RX.f"
	    } else {

/*              An 1x1 block. */

#line 212 "MB03RX.f"
		wr[l] = a[l + l * a_dim1];
#line 213 "MB03RX.f"
		wi[l] = 0.;
#line 214 "MB03RX.f"
		++l;
#line 215 "MB03RX.f"
	    }
#line 216 "MB03RX.f"
	    goto L20;
#line 217 "MB03RX.f"
	} else if (l == *n) {
#line 218 "MB03RX.f"
	    wr[l] = a[l + l * a_dim1];
#line 219 "MB03RX.f"
	    wi[l] = 0.;
#line 220 "MB03RX.f"
	}
/*        END WHILE 20 */
#line 222 "MB03RX.f"
    }

#line 224 "MB03RX.f"
    return 0;
/* *** Last line of MB03RX *** */
} /* mb03rx_ */
