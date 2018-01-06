#line 1 "MB04VX.f"
/* MB04VX.f -- translated by f2c (version 20100827).
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

#line 1 "MB04VX.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int mb04vx_(logical *updatq, logical *updatz, integer *m, 
	integer *n, integer *nblcks, integer *inuk, integer *imuk, doublereal 
	*a, integer *lda, doublereal *e, integer *lde, doublereal *q, integer 
	*ldq, doublereal *z__, integer *ldz, integer *mnei)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2, i__3;

    /* Local variables */
    static integer k, ca, ra;
    static doublereal sc;
    static integer ip;
    static doublereal ss;
    static integer tp1, cja, cje, rje, muk, nuk, mup, nup, mup1, minf, sk1p1, 
	    tk1p1, meps, neps, mukp1;
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

/*     NBLCKS  (input) INTEGER */
/*             The number of submatrices having full row rank (possibly */
/*             zero) in A(eps,inf). */

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

/*     MNEI    (output) INTEGER array, dimension (3) */
/*             MNEI(1) = MEPS =    row dimension of sE(eps)-A(eps); */
/*             MNEI(2) = NEPS = column dimension of sE(eps)-A(eps); */
/*             MNEI(3) = MINF = order of the regular pencil */
/*                              sE(inf)-A(inf). */

/*     REFERENCES */

/*     [1] Beelen, Th. */
/*         New Algorithms for Computing the Kronecker structure of a */
/*         Pencil with Applications to Systems and Control Theory. */
/*         Ph.D.Thesis, Eindhoven University of Technology, */
/*         The Netherlands, 1987. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Jan. 1998. */
/*     Based on Release 3.0 routine MB04TX modified by A. Varga, */
/*     German Aerospace Research Establishment, Oberpfaffenhofen, */
/*     Germany, Nov. 1997, as follows: */
/*     1) NBLCKS is only an input variable; */
/*     2) the significance of MNEI is changed. */

/*     REVISIONS */

/*     A. Varga, DLR Oberpfaffenhofen, March 2002. */

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

#line 208 "MB04VX.f"
    /* Parameter adjustments */
#line 208 "MB04VX.f"
    --inuk;
#line 208 "MB04VX.f"
    --imuk;
#line 208 "MB04VX.f"
    a_dim1 = *lda;
#line 208 "MB04VX.f"
    a_offset = 1 + a_dim1;
#line 208 "MB04VX.f"
    a -= a_offset;
#line 208 "MB04VX.f"
    e_dim1 = *lde;
#line 208 "MB04VX.f"
    e_offset = 1 + e_dim1;
#line 208 "MB04VX.f"
    e -= e_offset;
#line 208 "MB04VX.f"
    q_dim1 = *ldq;
#line 208 "MB04VX.f"
    q_offset = 1 + q_dim1;
#line 208 "MB04VX.f"
    q -= q_offset;
#line 208 "MB04VX.f"
    z_dim1 = *ldz;
#line 208 "MB04VX.f"
    z_offset = 1 + z_dim1;
#line 208 "MB04VX.f"
    z__ -= z_offset;
#line 208 "MB04VX.f"
    --mnei;
#line 208 "MB04VX.f"

#line 208 "MB04VX.f"
    /* Function Body */
#line 208 "MB04VX.f"
    mnei[1] = 0;
#line 209 "MB04VX.f"
    mnei[2] = 0;
#line 210 "MB04VX.f"
    mnei[3] = 0;
#line 211 "MB04VX.f"
    if (*m <= 0 || *n <= 0) {
#line 211 "MB04VX.f"
	return 0;
#line 211 "MB04VX.f"
    }

/*     Initialisation. */

#line 216 "MB04VX.f"
    ismuk = 0;
#line 217 "MB04VX.f"
    isnuk = 0;

#line 219 "MB04VX.f"
    i__1 = *nblcks;
#line 219 "MB04VX.f"
    for (k = 1; k <= i__1; ++k) {
#line 220 "MB04VX.f"
	ismuk += imuk[k];
#line 221 "MB04VX.f"
	isnuk += inuk[k];
#line 222 "MB04VX.f"
/* L20: */
#line 222 "MB04VX.f"
    }

/*     MEPS, NEPS are the dimensions of the pencil s*E(eps)-A(eps). */
/*     MEPS = Sum(k=1,...,nblcks) NU(k), */
/*     NEPS = Sum(k=1,...,nblcks) MU(k). */
/*     MINF is the order of the regular pencil s*E(inf)-A(inf). */

#line 229 "MB04VX.f"
    meps = isnuk;
#line 230 "MB04VX.f"
    neps = ismuk;
#line 231 "MB04VX.f"
    minf = 0;

/*     MUKP1 = mu(k+1).  N.B. It is assumed that mu(NBLCKS + 1) = 0. */

#line 235 "MB04VX.f"
    mukp1 = 0;

#line 237 "MB04VX.f"
    for (k = *nblcks; k >= 1; --k) {
#line 238 "MB04VX.f"
	nuk = inuk[k];
#line 239 "MB04VX.f"
	muk = imuk[k];

/*        Reduce submatrix E(k,k+1) to square matrix. */
/*        NOTE that always NU(k) >= MU(k+1) >= 0. */

/*        WHILE ( NU(k) >  MU(k+1) ) DO */
#line 245 "MB04VX.f"
L40:
#line 245 "MB04VX.f"
	if (nuk > mukp1) {

/*           sk1p1 = sum(i=k+1,...,p-1) NU(i) */
/*           tk1p1 = sum(i=k+1,...,p-1) MU(i) */
/*           ismuk = sum(i=1,...,k) MU(i) */
/*           tp1   = sum(i=1,...,p-1) MU(i) = ismuk + tk1p1. */

#line 252 "MB04VX.f"
	    sk1p1 = 0;
#line 253 "MB04VX.f"
	    tk1p1 = 0;

#line 255 "MB04VX.f"
	    i__1 = *nblcks;
#line 255 "MB04VX.f"
	    for (ip = k + 1; ip <= i__1; ++ip) {

/*              Annihilate the elements originally present in the last */
/*              row of E(k,p+1) and A(k,p). */
/*              Start annihilating the first MU(p) - MU(p+1) elements by */
/*              applying column Givens rotations plus interchanging */
/*              elements. */
/*              Use original bottom diagonal element of A(k,k) as pivot. */
/*              Start position of pivot in A = (ra,ca). */

#line 265 "MB04VX.f"
		tp1 = ismuk + tk1p1;
#line 266 "MB04VX.f"
		ra = isnuk + sk1p1;
#line 267 "MB04VX.f"
		ca = tp1;

#line 269 "MB04VX.f"
		mup = imuk[ip];
#line 270 "MB04VX.f"
		nup = inuk[ip];
#line 271 "MB04VX.f"
		mup1 = nup;

#line 273 "MB04VX.f"
		i__2 = ca + mup - nup - 1;
#line 273 "MB04VX.f"
		for (cja = ca; cja <= i__2; ++cja) {

/*                 CJA = current column index of pivot in A. */

#line 277 "MB04VX.f"
		    drotg_(&a[ra + cja * a_dim1], &a[ra + (cja + 1) * a_dim1],
			     &sc, &ss);

/*                 Apply transformations to A- and E-matrix. */
/*                 Interchange columns simultaneously. */
/*                 Update column transformation matrix Z, if needed. */

#line 283 "MB04VX.f"
		    i__3 = ra - 1;
#line 283 "MB04VX.f"
		    mb04tu_(&i__3, &a[cja * a_dim1 + 1], &c__1, &a[(cja + 1) *
			     a_dim1 + 1], &c__1, &sc, &ss);
#line 285 "MB04VX.f"
		    a[ra + (cja + 1) * a_dim1] = a[ra + cja * a_dim1];
#line 286 "MB04VX.f"
		    a[ra + cja * a_dim1] = 0.;
#line 287 "MB04VX.f"
		    mb04tu_(&ra, &e[cja * e_dim1 + 1], &c__1, &e[(cja + 1) * 
			    e_dim1 + 1], &c__1, &sc, &ss);
#line 288 "MB04VX.f"
		    if (*updatz) {
#line 288 "MB04VX.f"
			mb04tu_(n, &z__[cja * z_dim1 + 1], &c__1, &z__[(cja + 
				1) * z_dim1 + 1], &c__1, &sc, &ss);
#line 288 "MB04VX.f"
		    }
#line 290 "MB04VX.f"
/* L60: */
#line 290 "MB04VX.f"
		}

/*              Annihilate the remaining elements originally present in */
/*              the last row of E(k,p+1) and A(k,p) by alternatingly */
/*              applying row and column rotations plus interchanging */
/*              elements. */
/*              Use diagonal elements of E(p,p+1) and original bottom */
/*              diagonal element of A(k,k) as pivots, respectively. */
/*              (re,ce) and (ra,ca) are the starting positions of the */
/*              pivots in E and A. */

#line 301 "MB04VX.f"
		cje = tp1 + mup;
#line 302 "MB04VX.f"
		cja = cje - mup1 - 1;

#line 304 "MB04VX.f"
		i__2 = ra + mup1;
#line 304 "MB04VX.f"
		for (rje = ra + 1; rje <= i__2; ++rje) {

/*                 (RJE,CJE) = current position pivot in E. */

#line 308 "MB04VX.f"
		    ++cje;
#line 309 "MB04VX.f"
		    ++cja;

/*                 Determine the row transformations. */
/*                 Apply these transformations to E- and A-matrix. */
/*                 Interchange the rows simultaneously. */
/*                 Update row transformation matrix Q, if needed. */

#line 316 "MB04VX.f"
		    drotg_(&e[rje + cje * e_dim1], &e[rje - 1 + cje * e_dim1],
			     &sc, &ss);
#line 317 "MB04VX.f"
		    i__3 = *n - cje;
#line 317 "MB04VX.f"
		    mb04tu_(&i__3, &e[rje + (cje + 1) * e_dim1], lde, &e[rje 
			    - 1 + (cje + 1) * e_dim1], lde, &sc, &ss);
#line 319 "MB04VX.f"
		    e[rje - 1 + cje * e_dim1] = e[rje + cje * e_dim1];
#line 320 "MB04VX.f"
		    e[rje + cje * e_dim1] = 0.;
#line 321 "MB04VX.f"
		    i__3 = *n - cja + 1;
#line 321 "MB04VX.f"
		    mb04tu_(&i__3, &a[rje + cja * a_dim1], lda, &a[rje - 1 + 
			    cja * a_dim1], lda, &sc, &ss);
#line 323 "MB04VX.f"
		    if (*updatq) {
#line 323 "MB04VX.f"
			mb04tu_(m, &q[rje * q_dim1 + 1], &c__1, &q[(rje - 1) *
				 q_dim1 + 1], &c__1, &sc, &ss);
#line 323 "MB04VX.f"
		    }

/*                 Determine the column transformations. */
/*                 Apply these transformations to A- and E-matrix. */
/*                 Interchange the columns simultaneously. */
/*                 Update column transformation matrix Z, if needed. */

#line 331 "MB04VX.f"
		    drotg_(&a[rje + cja * a_dim1], &a[rje + (cja + 1) * 
			    a_dim1], &sc, &ss);
#line 332 "MB04VX.f"
		    i__3 = rje - 1;
#line 332 "MB04VX.f"
		    mb04tu_(&i__3, &a[cja * a_dim1 + 1], &c__1, &a[(cja + 1) *
			     a_dim1 + 1], &c__1, &sc, &ss);
#line 334 "MB04VX.f"
		    a[rje + (cja + 1) * a_dim1] = a[rje + cja * a_dim1];
#line 335 "MB04VX.f"
		    a[rje + cja * a_dim1] = 0.;
#line 336 "MB04VX.f"
		    mb04tu_(&rje, &e[cja * e_dim1 + 1], &c__1, &e[(cja + 1) * 
			    e_dim1 + 1], &c__1, &sc, &ss);
#line 337 "MB04VX.f"
		    if (*updatz) {
#line 337 "MB04VX.f"
			mb04tu_(n, &z__[cja * z_dim1 + 1], &c__1, &z__[(cja + 
				1) * z_dim1 + 1], &c__1, &sc, &ss);
#line 337 "MB04VX.f"
		    }
#line 339 "MB04VX.f"
/* L80: */
#line 339 "MB04VX.f"
		}

#line 341 "MB04VX.f"
		sk1p1 += nup;
#line 342 "MB04VX.f"
		tk1p1 += mup;

#line 344 "MB04VX.f"
/* L100: */
#line 344 "MB04VX.f"
	    }

/*           Reduce A=A(eps,inf) and E=E(eps,inf) by ignoring their last */
/*           row and right most column. The row and column ignored */
/*           belong to the pencil s*E(inf)-A(inf). */
/*           Redefine blocks in new A and E. */

#line 351 "MB04VX.f"
	    --muk;
#line 352 "MB04VX.f"
	    --nuk;
#line 353 "MB04VX.f"
	    --ismuk;
#line 354 "MB04VX.f"
	    --isnuk;
#line 355 "MB04VX.f"
	    --meps;
#line 356 "MB04VX.f"
	    --neps;
#line 357 "MB04VX.f"
	    ++minf;

#line 359 "MB04VX.f"
	    goto L40;
#line 360 "MB04VX.f"
	}
/*        END WHILE 40 */

#line 363 "MB04VX.f"
	imuk[k] = muk;
#line 364 "MB04VX.f"
	inuk[k] = nuk;

/*        Now submatrix E(k,k+1) is square. */

/*        Consider next submatrix (k:=k-1). */

#line 370 "MB04VX.f"
	isnuk -= nuk;
#line 371 "MB04VX.f"
	ismuk -= muk;
#line 372 "MB04VX.f"
	mukp1 = muk;
#line 373 "MB04VX.f"
/* L120: */
#line 373 "MB04VX.f"
    }

/*     Store dimensions of the pencils s*E(eps)-A(eps) and */
/*     s*E(inf)-A(inf) in array MNEI. */

#line 378 "MB04VX.f"
    mnei[1] = meps;
#line 379 "MB04VX.f"
    mnei[2] = neps;
#line 380 "MB04VX.f"
    mnei[3] = minf;

#line 382 "MB04VX.f"
    return 0;
/* *** Last line of MB04VX *** */
} /* mb04vx_ */

