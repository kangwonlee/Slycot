#line 1 "MB04TX.f"
/* MB04TX.f -- translated by f2c (version 20100827).
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

#line 1 "MB04TX.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb04tx_(logical *updatq, logical *updatz, integer *m, 
	integer *n, integer *nblcks, integer *inuk, integer *imuk, doublereal 
	*a, integer *lda, doublereal *e, integer *lde, doublereal *q, integer 
	*ldq, doublereal *z__, integer *ldz, integer *mnei)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer k, ca, ce, ra;
    static doublereal sc;
    static integer ip;
    static doublereal ss;
    static integer tp1, cja, cje, rje, muk, nuk, mup, nup, mup1, minf, ninf, 
	    sk1p1, tk1p1, meps, neps, mukp1;
    extern /* Subroutine */ int mb04tu_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), drotg_(
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer ismuk, isnuk;


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

/*     To separate the pencils s*E(eps)-A(eps) and s*E(inf)-A(inf) in */
/*     s*E(eps,inf)-A(eps,inf) using Algorithm 3.3.3 in [1]. */

/*     On entry, it is assumed that the M-by-N matrices A and E have */
/*     been obtained after applying the Algorithms 3.2.1 and 3.3.1 to */
/*     the pencil s*E - A as described in [1], i.e. */

/*                        | s*E(eps,inf)-A(eps,inf) |      X      | */
/*        Q'(s*E - A)Z  = |-------------------------|-------------| */
/*                        |             0           | s*E(r)-A(r) | */

/*     Here the pencil s*E(eps,inf)-A(eps,inf) is in staircase form. */
/*     This pencil contains all Kronecker column indices and infinite */
/*     elementary divisors of the pencil s*E - A. */
/*     The pencil s*E(r)-A(r) contains all Kronecker row indices and */
/*     finite elementary divisors of s*E - A. */
/*     Furthermore, the submatrices having full row and column rank in */
/*     the pencil s*E(eps,inf)-A(eps,inf) are assumed to be */
/*     triangularized. */

/*     On exit, the result then is */

/*                        Q'(s*E - A)Z = */

/*          | s*E(eps)-A(eps) |        X        |      X      | */
/*          |-----------------|-----------------|-------------| */
/*          |        0        | s*E(inf)-A(inf) |      X      | */
/*          |===================================|=============| */
/*          |                                   |             | */
/*          |                 0                 | s*E(r)-A(r) | */

/*     Note that the pencil s*E(r)-A(r) is not reduced further. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     UPDATQ  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Q the orthogonal row transformations, as follows: */
/*             = .FALSE.: Do not form Q; */
/*             = .TRUE.:  The given matrix Q is updated by the orthogonal */
/*                        row transformations used in the reduction. */

/*     UPDATZ  LOGICAL */
/*             Indicates whether the user wishes to accumulate in a */
/*             matrix Z the orthogonal column transformations, as */
/*             follows: */
/*             = .FALSE.: Do not form Z; */
/*             = .TRUE.:  The given matrix Z is updated by the orthogonal */
/*                        column transformations used in the reduction. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             Number of rows of A and E.  M >= 0. */

/*     N       (input) INTEGER */
/*             Number of columns of A and E.  N >= 0. */

/*     NBLCKS  (input/output) INTEGER */
/*             On entry, the number of submatrices having full row rank */
/*             (possibly zero) in A(eps,inf). */
/*             On exit, the input value has been reduced by one, if the */
/*             last submatrix is a 0-by-0 (empty) matrix. */

/*     INUK    (input/output) INTEGER array, dimension (NBLCKS) */
/*             On entry, this array contains the row dimensions nu(k), */
/*             (k=1, 2, ..., NBLCKS) of the submatrices having full row */
/*             rank in the pencil s*E(eps,inf)-A(eps,inf). */
/*             On exit, this array contains the row dimensions nu(k), */
/*             (k=1, 2, ..., NBLCKS) of the submatrices having full row */
/*             rank in the pencil s*E(eps)-A(eps). */

/*     IMUK    (input/output) INTEGER array, dimension (NBLCKS) */
/*             On entry, this array contains the column dimensions mu(k), */
/*             (k=1, 2, ..., NBLCKS) of the submatrices having full */
/*             column rank in the pencil s*E(eps,inf)-A(eps,inf). */
/*             On exit, this array contains the column dimensions mu(k), */
/*             (k=1, 2, ..., NBLCKS) of the submatrices having full */
/*             column rank in the pencil s*E(eps)-A(eps). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, this array contains the matrix A to be reduced. */
/*             On exit, it contains the transformed matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, this array contains the matrix E to be reduced. */
/*             On exit, it contains the transformed matrix E. */

/*     LDE     INTEGER */
/*             The leading dimension of array E.  LDE >= MAX(1,M). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,*) */
/*             On entry, if UPDATQ = .TRUE., then the leading M-by-M */
/*             part of this array must contain a given matrix Q (e.g. */
/*             from a previous call to another SLICOT routine), and on */
/*             exit, the leading M-by-M part of this array contains the */
/*             product of the input matrix Q and the row transformation */
/*             matrix that has transformed the rows of the matrices A */
/*             and E. */
/*             If UPDATQ = .FALSE., the array Q is not referenced and */
/*             can be supplied as a dummy array (i.e. set parameter */
/*             LDQ = 1 and declare this array to be Q(1,1) in the calling */
/*             program). */

/*     LDQ     INTEGER */
/*             The leading dimension of array Q. If UPDATQ = .TRUE., */
/*             LDQ >= MAX(1,M); if UPDATQ = .FALSE., LDQ >= 1. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,*) */
/*             On entry, if UPDATZ = .TRUE., then the leading N-by-N */
/*             part of this array must contain a given matrix Z (e.g. */
/*             from a previous call to another SLICOT routine), and on */
/*             exit, the leading N-by-N part of this array contains the */
/*             product of the input matrix Z and the column */
/*             transformation matrix that has transformed the columns of */
/*             the matrices A and E. */
/*             If UPDATZ = .FALSE., the array Z is not referenced and */
/*             can be supplied as a dummy array (i.e. set parameter */
/*             LDZ = 1 and declare this array to be Z(1,1) in the calling */
/*             program). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z. If UPDATZ = .TRUE., */
/*             LDZ >= MAX(1,N); if UPDATZ = .FALSE., LDZ >= 1. */

/*     MNEI    (output) INTEGER array, dimension (4) */
/*             MNEI(1) = MEPS = row    dimension of s*E(eps)-A(eps), */
/*             MNEI(2) = NEPS = column dimension of s*E(eps)-A(eps), */
/*             MNEI(3) = MINF = row    dimension of s*E(inf)-A(inf), */
/*             MNEI(4) = NINF = column dimension of s*E(inf)-A(inf). */

/*     REFERENCES */

/*     [1] Beelen, Th. */
/*         New Algorithms for Computing the Kronecker structure of a */
/*         Pencil with Applications to Systems and Control Theory. */
/*         Ph.D.Thesis, Eindhoven University of Technology, */
/*         The Netherlands, 1987. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Apr. 1997. */
/*     Supersedes Release 2.0 routine MB04FX by Th.G.J. Beelen, */
/*     Philips Glass Eindhoven, Holland. */

/*     REVISIONS */

/*     June 13, 1997, V. Sima. */
/*     November 24, 1997, A. Varga: initialization of MNEI to 0, instead */
/*                                  of ZERO. */

/*     KEYWORDS */

/*     Generalized eigenvalue problem, Kronecker indices, orthogonal */
/*     transformation, staircase form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Executable Statements .. */

#line 209 "MB04TX.f"
    /* Parameter adjustments */
#line 209 "MB04TX.f"
    --inuk;
#line 209 "MB04TX.f"
    --imuk;
#line 209 "MB04TX.f"
    a_dim1 = *lda;
#line 209 "MB04TX.f"
    a_offset = 1 + a_dim1;
#line 209 "MB04TX.f"
    a -= a_offset;
#line 209 "MB04TX.f"
    e_dim1 = *lde;
#line 209 "MB04TX.f"
    e_offset = 1 + e_dim1;
#line 209 "MB04TX.f"
    e -= e_offset;
#line 209 "MB04TX.f"
    q_dim1 = *ldq;
#line 209 "MB04TX.f"
    q_offset = 1 + q_dim1;
#line 209 "MB04TX.f"
    q -= q_offset;
#line 209 "MB04TX.f"
    z_dim1 = *ldz;
#line 209 "MB04TX.f"
    z_offset = 1 + z_dim1;
#line 209 "MB04TX.f"
    z__ -= z_offset;
#line 209 "MB04TX.f"
    --mnei;
#line 209 "MB04TX.f"

#line 209 "MB04TX.f"
    /* Function Body */
#line 209 "MB04TX.f"
    mnei[1] = 0;
#line 210 "MB04TX.f"
    mnei[2] = 0;
#line 211 "MB04TX.f"
    mnei[3] = 0;
#line 212 "MB04TX.f"
    mnei[4] = 0;
#line 213 "MB04TX.f"
    if (*m <= 0 || *n <= 0) {
#line 213 "MB04TX.f"
	return 0;
#line 213 "MB04TX.f"
    }

/*     Initialisation. */

#line 218 "MB04TX.f"
    ismuk = 0;
#line 219 "MB04TX.f"
    isnuk = 0;

#line 221 "MB04TX.f"
    i__1 = *nblcks;
#line 221 "MB04TX.f"
    for (k = 1; k <= i__1; ++k) {
#line 222 "MB04TX.f"
	ismuk += imuk[k];
#line 223 "MB04TX.f"
	isnuk += inuk[k];
#line 224 "MB04TX.f"
/* L20: */
#line 224 "MB04TX.f"
    }

/*     MEPS, NEPS are the dimensions of the pencil s*E(eps)-A(eps). */
/*     MEPS = Sum(k=1,...,nblcks) NU(k), */
/*     NEPS = Sum(k=1,...,nblcks) MU(k). */
/*     MINF, NINF are the dimensions of the pencil s*E(inf)-A(inf). */

#line 231 "MB04TX.f"
    meps = isnuk;
#line 232 "MB04TX.f"
    neps = ismuk;
#line 233 "MB04TX.f"
    minf = 0;
#line 234 "MB04TX.f"
    ninf = 0;

/*     MUKP1 = mu(k+1).  N.B. It is assumed that mu(NBLCKS + 1) = 0. */

#line 238 "MB04TX.f"
    mukp1 = 0;

#line 240 "MB04TX.f"
    for (k = *nblcks; k >= 1; --k) {
#line 241 "MB04TX.f"
	nuk = inuk[k];
#line 242 "MB04TX.f"
	muk = imuk[k];

/*        Reduce submatrix E(k,k+1) to square matrix. */
/*        NOTE that always NU(k) >= MU(k+1) >= 0. */

/*        WHILE ( NU(k) >  MU(k+1) ) DO */
#line 248 "MB04TX.f"
L40:
#line 248 "MB04TX.f"
	if (nuk > mukp1) {

/*           sk1p1 = sum(i=k+1,...,p-1) NU(i) */
/*           tk1p1 = sum(i=k+1,...,p-1) MU(i) */
/*           ismuk = sum(i=1,...,k) MU(i) */
/*           tp1   = sum(i=1,...,p-1) MU(i) = ismuk + tk1p1. */

#line 255 "MB04TX.f"
	    sk1p1 = 0;
#line 256 "MB04TX.f"
	    tk1p1 = 0;

#line 258 "MB04TX.f"
	    i__1 = *nblcks;
#line 258 "MB04TX.f"
	    for (ip = k + 1; ip <= i__1; ++ip) {

/*              Annihilate the elements originally present in the last */
/*              row of E(k,p+1) and A(k,p). */
/*              Start annihilating the first MU(p) - MU(p+1) elements by */
/*              applying column Givens rotations plus interchanging */
/*              elements. */
/*              Use original bottom diagonal element of A(k,k) as pivot. */
/*              Start position of pivot in A = (ra,ca). */

#line 268 "MB04TX.f"
		tp1 = ismuk + tk1p1;
#line 269 "MB04TX.f"
		ra = isnuk + sk1p1;
#line 270 "MB04TX.f"
		ca = tp1;

#line 272 "MB04TX.f"
		mup = imuk[ip];
#line 273 "MB04TX.f"
		nup = inuk[ip];
#line 274 "MB04TX.f"
		mup1 = nup;

#line 276 "MB04TX.f"
		i__2 = ca + mup - nup - 1;
#line 276 "MB04TX.f"
		for (cja = ca; cja <= i__2; ++cja) {

/*                 CJA = current column index of pivot in A. */

#line 280 "MB04TX.f"
		    drotg_(&a[ra + cja * a_dim1], &a[ra + (cja + 1) * a_dim1],
			     &sc, &ss);

/*                 Apply transformations to A- and E-matrix. */
/*                 Interchange columns simultaneously. */
/*                 Update column transformation matrix Z, if needed. */

#line 286 "MB04TX.f"
		    i__3 = ra - 1;
#line 286 "MB04TX.f"
		    mb04tu_(&i__3, &a[cja * a_dim1 + 1], &c__1, &a[(cja + 1) *
			     a_dim1 + 1], &c__1, &sc, &ss);
#line 288 "MB04TX.f"
		    a[ra + (cja + 1) * a_dim1] = a[ra + cja * a_dim1];
#line 289 "MB04TX.f"
		    a[ra + cja * a_dim1] = 0.;
#line 290 "MB04TX.f"
		    mb04tu_(&ra, &e[cja * e_dim1 + 1], &c__1, &e[(cja + 1) * 
			    e_dim1 + 1], &c__1, &sc, &ss);
#line 291 "MB04TX.f"
		    if (*updatz) {
#line 291 "MB04TX.f"
			mb04tu_(n, &z__[cja * z_dim1 + 1], &c__1, &z__[(cja + 
				1) * z_dim1 + 1], &c__1, &sc, &ss);
#line 291 "MB04TX.f"
		    }
#line 293 "MB04TX.f"
/* L60: */
#line 293 "MB04TX.f"
		}

/*              Annihilate the remaining elements originally present in */
/*              the last row of E(k,p+1) and A(k,p) by alternatingly */
/*              applying row and column rotations plus interchanging */
/*              elements. */
/*              Use diagonal elements of E(p,p+1) and original bottom */
/*              diagonal element of A(k,k) as pivots, respectively. */
/*              (re,ce) and (ra,ca) are the starting positions of the */
/*              pivots in E and A. */

#line 304 "MB04TX.f"
		ce = tp1 + mup;
#line 305 "MB04TX.f"
		ca = ce - mup1 - 1;

#line 307 "MB04TX.f"
		i__2 = ra + mup1;
#line 307 "MB04TX.f"
		for (rje = ra + 1; rje <= i__2; ++rje) {

/*                 (RJE,CJE) = current position pivot in E. */

#line 311 "MB04TX.f"
		    cje = ce + 1;
#line 312 "MB04TX.f"
		    cja = ca + 1;

/*                 Determine the row transformations. */
/*                 Apply these transformations to E- and A-matrix. */
/*                 Interchange the rows simultaneously. */
/*                 Update row transformation matrix Q, if needed. */

#line 319 "MB04TX.f"
		    drotg_(&e[rje + cje * e_dim1], &e[rje - 1 + cje * e_dim1],
			     &sc, &ss);
#line 320 "MB04TX.f"
		    i__3 = *n - cje;
#line 320 "MB04TX.f"
		    mb04tu_(&i__3, &e[rje + (cje + 1) * e_dim1], lde, &e[rje 
			    - 1 + (cje + 1) * e_dim1], lde, &sc, &ss);
#line 322 "MB04TX.f"
		    e[rje - 1 + cje * e_dim1] = e[rje + cje * e_dim1];
#line 323 "MB04TX.f"
		    e[rje + cje * e_dim1] = 0.;
#line 324 "MB04TX.f"
		    i__3 = *n - cja + 1;
#line 324 "MB04TX.f"
		    mb04tu_(&i__3, &a[rje + cja * a_dim1], lda, &a[rje - 1 + 
			    cja * a_dim1], lda, &sc, &ss);
#line 326 "MB04TX.f"
		    if (*updatq) {
#line 326 "MB04TX.f"
			mb04tu_(m, &q[rje * q_dim1 + 1], &c__1, &q[(rje - 1) *
				 q_dim1 + 1], &c__1, &sc, &ss);
#line 326 "MB04TX.f"
		    }

/*                 Determine the column transformations. */
/*                 Apply these transformations to A- and E-matrix. */
/*                 Interchange the columns simultaneously. */
/*                 Update column transformation matrix Z, if needed. */

#line 334 "MB04TX.f"
		    drotg_(&a[rje + cja * a_dim1], &a[rje + (cja + 1) * 
			    a_dim1], &sc, &ss);
#line 335 "MB04TX.f"
		    i__3 = rje - 1;
#line 335 "MB04TX.f"
		    mb04tu_(&i__3, &a[cja * a_dim1 + 1], &c__1, &a[(cja + 1) *
			     a_dim1 + 1], &c__1, &sc, &ss);
#line 337 "MB04TX.f"
		    a[rje + (cja + 1) * a_dim1] = a[rje + cja * a_dim1];
#line 338 "MB04TX.f"
		    a[rje + cja * a_dim1] = 0.;
#line 339 "MB04TX.f"
		    mb04tu_(&rje, &e[cja * e_dim1 + 1], &c__1, &e[(cja + 1) * 
			    e_dim1 + 1], &c__1, &sc, &ss);
#line 340 "MB04TX.f"
		    if (*updatz) {
#line 340 "MB04TX.f"
			mb04tu_(n, &z__[cja * z_dim1 + 1], &c__1, &z__[(cja + 
				1) * z_dim1 + 1], &c__1, &sc, &ss);
#line 340 "MB04TX.f"
		    }
#line 342 "MB04TX.f"
/* L80: */
#line 342 "MB04TX.f"
		}

#line 344 "MB04TX.f"
		sk1p1 += nup;
#line 345 "MB04TX.f"
		tk1p1 += mup;

#line 347 "MB04TX.f"
/* L100: */
#line 347 "MB04TX.f"
	    }

/*           Reduce A=A(eps,inf) and E=E(eps,inf) by ignoring their last */
/*           row and right most column. The row and column ignored */
/*           belong to the pencil s*E(inf)-A(inf). */
/*           Redefine blocks in new A and E. */

#line 354 "MB04TX.f"
	    --muk;
#line 355 "MB04TX.f"
	    --nuk;
#line 356 "MB04TX.f"
	    --ismuk;
#line 357 "MB04TX.f"
	    --isnuk;
#line 358 "MB04TX.f"
	    --meps;
#line 359 "MB04TX.f"
	    --neps;
#line 360 "MB04TX.f"
	    ++minf;
#line 361 "MB04TX.f"
	    ++ninf;

#line 363 "MB04TX.f"
	    goto L40;
#line 364 "MB04TX.f"
	}
/*        END WHILE 40 */

#line 367 "MB04TX.f"
	imuk[k] = muk;
#line 368 "MB04TX.f"
	inuk[k] = nuk;

/*        Now submatrix E(k,k+1) is square. */

/*        Consider next submatrix (k:=k-1). */

#line 374 "MB04TX.f"
	isnuk -= nuk;
#line 375 "MB04TX.f"
	ismuk -= muk;
#line 376 "MB04TX.f"
	mukp1 = muk;
#line 377 "MB04TX.f"
/* L120: */
#line 377 "MB04TX.f"
    }

/*     If mu(NBLCKS) = 0, then the last submatrix counted in NBLCKS is */
/*     a 0-by-0 (empty) matrix. This "matrix" must be removed. */

#line 382 "MB04TX.f"
    if (imuk[*nblcks] == 0) {
#line 382 "MB04TX.f"
	--(*nblcks);
#line 382 "MB04TX.f"
    }

/*     Store dimensions of the pencils s*E(eps)-A(eps) and */
/*     s*E(inf)-A(inf) in array MNEI. */

#line 387 "MB04TX.f"
    mnei[1] = meps;
#line 388 "MB04TX.f"
    mnei[2] = neps;
#line 389 "MB04TX.f"
    mnei[3] = minf;
#line 390 "MB04TX.f"
    mnei[4] = ninf;

#line 392 "MB04TX.f"
    return 0;
/* *** Last line of MB04TX *** */
} /* mb04tx_ */

