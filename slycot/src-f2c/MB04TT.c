#line 1 "MB04TT.f"
/* MB04TT.f -- translated by f2c (version 20100827).
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

#line 1 "MB04TT.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b7 = 0.;
static logical c_false = FALSE_;

/* Subroutine */ int mb04tt_(logical *updatq, logical *updatz, integer *m, 
	integer *n, integer *ifira, integer *ifica, integer *nca, doublereal *
	a, integer *lda, doublereal *e, integer *lde, doublereal *q, integer *
	ldq, doublereal *z__, integer *ldz, integer *istair, integer *rank, 
	doublereal *tol, integer *iwork)
{
    /* System generated locals */
    integer a_dim1, a_offset, e_dim1, e_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2;
    doublereal d__1;

    /* Local variables */
    static integer i__, k, l, ii, kk, mj, ll, ip, nj;
    static doublereal sc, ss;
    static integer jc1, jc2, mk1;
    static doublereal bmx;
    static integer ist1, ist2, lsav;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer jpvt;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), drotg_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer itype;
    static logical lzero;
    static integer ifica1, ifira1;
    extern integer idamax_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlapmt_(logical *, integer *, integer *, doublereal *, integer *, 
	    integer *);
    static integer mxrank;
    static doublereal eijpvt, bmxnrm;
    static integer istpvt;


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

/*     Let A and E be M-by-N matrices with E in column echelon form. */
/*     Let AA and EE be the following submatrices of A and E: */
/*       AA := A(IFIRA : M ; IFICA : N) */
/*       EE := E(IFIRA : M ; IFICA : N). */
/*     Let Aj and Ej be the following submatrices of AA and EE: */
/*       Aj := A(IFIRA : M ; IFICA : IFICA + NCA - 1) and */
/*       Ej := E(IFIRA : M ; IFICA + NCA : N). */

/*     To transform (AA,EE) such that Aj is row compressed while keeping */
/*     matrix Ej in column echelon form (which may be different from the */
/*     form on entry). */
/*     In fact the routine performs the j-th step of Algorithm 3.2.1 in */
/*     [1]. Furthermore, it determines the rank RANK of the submatrix Ej, */
/*     which is equal to the number of corner points in submatrix Ej. */

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
/*             M is the number of rows of the matrices A, E and Q. */
/*             M >= 0. */

/*     N       (input) INTEGER */
/*             N is the number of columns of the matrices A, E and Z. */
/*             N >= 0. */

/*     IFIRA   (input) INTEGER */
/*             IFIRA is the first row index of the submatrices Aj and Ej */
/*             in the matrices A and E, respectively. */

/*     IFICA   (input) INTEGER */
/*             IFICA and IFICA + NCA are the first column indices of the */
/*             submatrices Aj and Ej in the matrices A and E, */
/*             respectively. */

/*     NCA     (input) INTEGER */
/*             NCA is the number of columns of the submatrix Aj in A. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, A(IFIRA : M ; IFICA : IFICA + NCA - 1) contains */
/*             the matrix Aj. */
/*             On exit, it contains the matrix A with AA that has been */
/*             row compressed while keeping EE in column echelon form. */

/*     LDA     INTEGER */
/*             The leading dimension of array A. LDA >= MAX(1,M). */

/*     E       (input/output) DOUBLE PRECISION array, dimension (LDE,N) */
/*             On entry, E(IFIRA : M ; IFICA + NCA : N) contains the */
/*             matrix Ej which is in column echelon form. */
/*             On exit, it contains the transformed matrix EE which is */
/*             kept in column echelon form. */

/*     LDE     INTEGER */
/*             The leading dimension of array E. LDE >= MAX(1,M). */

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

/*     ISTAIR  (input/output) INTEGER array, dimension (M) */
/*             On entry, ISTAIR contains information on the column */
/*             echelon form of the input matrix E as follows: */
/*             ISTAIR(i) = +j: the boundary element E(i,j) is a corner */
/*                             point; */
/*                         -j: the boundary element E(i,j) is not a */
/*                             corner point (where i=1,...,M). */
/*             On exit, ISTAIR contains the same information for the */
/*             transformed matrix E. */

/*     RANK    (output) INTEGER */
/*             Numerical rank of the submatrix Aj in A (based on TOL). */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance used when considering matrix elements */
/*             to be zero. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     REFERENCES */

/*     [1] Beelen, Th. */
/*         New Algorithms for Computing the Kronecker structure of a */
/*         Pencil with Applications to Systems and Control Theory. */
/*         Ph.D.Thesis, Eindhoven University of Technology, */
/*         The Netherlands, 1987. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Mar. 1997. */
/*     Supersedes Release 2.0 routine MB04FZ by Th.G.J. Beelen, */
/*     Philips Glass Eindhoven, Holland. */

/*     REVISIONS */

/*     June 13, 1997, V. Sima. */
/*     November 24, 1997, A. Varga: array starting point A(KK,LL) */
/*                                  correctly set when calling DLASET. */

/*     KEYWORDS */

/*     Echelon form, orthogonal transformation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 212 "MB04TT.f"
    /* Parameter adjustments */
#line 212 "MB04TT.f"
    a_dim1 = *lda;
#line 212 "MB04TT.f"
    a_offset = 1 + a_dim1;
#line 212 "MB04TT.f"
    a -= a_offset;
#line 212 "MB04TT.f"
    e_dim1 = *lde;
#line 212 "MB04TT.f"
    e_offset = 1 + e_dim1;
#line 212 "MB04TT.f"
    e -= e_offset;
#line 212 "MB04TT.f"
    q_dim1 = *ldq;
#line 212 "MB04TT.f"
    q_offset = 1 + q_dim1;
#line 212 "MB04TT.f"
    q -= q_offset;
#line 212 "MB04TT.f"
    z_dim1 = *ldz;
#line 212 "MB04TT.f"
    z_offset = 1 + z_dim1;
#line 212 "MB04TT.f"
    z__ -= z_offset;
#line 212 "MB04TT.f"
    --istair;
#line 212 "MB04TT.f"
    --iwork;
#line 212 "MB04TT.f"

#line 212 "MB04TT.f"
    /* Function Body */
#line 212 "MB04TT.f"
    *rank = 0;
#line 213 "MB04TT.f"
    if (*m <= 0 || *n <= 0) {
#line 213 "MB04TT.f"
	return 0;
#line 213 "MB04TT.f"
    }

/*     Initialisation. */

/*     NJ = number of columns in submatrix Aj, */
/*     MJ = number of rows in submatrices Aj and Ej. */

#line 221 "MB04TT.f"
    nj = *nca;
#line 222 "MB04TT.f"
    mj = *m + 1 - *ifira;
#line 223 "MB04TT.f"
    ifira1 = *ifira - 1;
#line 224 "MB04TT.f"
    ifica1 = *ifica - 1;

#line 226 "MB04TT.f"
    i__1 = nj;
#line 226 "MB04TT.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 227 "MB04TT.f"
	iwork[i__] = i__;
#line 228 "MB04TT.f"
/* L20: */
#line 228 "MB04TT.f"
    }

#line 230 "MB04TT.f"
    k = 1;
#line 231 "MB04TT.f"
    lzero = FALSE_;
#line 232 "MB04TT.f"
    *rank = min(nj,mj);
#line 233 "MB04TT.f"
    mxrank = *rank;

/*     WHILE ( K <= MXRANK ) and ( LZERO = FALSE ) DO */
#line 236 "MB04TT.f"
L40:
#line 236 "MB04TT.f"
    if (k <= mxrank && ! lzero) {

/*        Determine column in Aj with largest max-norm. */

#line 240 "MB04TT.f"
	bmxnrm = 0.;
#line 241 "MB04TT.f"
	lsav = k;
#line 242 "MB04TT.f"
	kk = ifira1 + k;

#line 244 "MB04TT.f"
	i__1 = nj;
#line 244 "MB04TT.f"
	for (l = k; l <= i__1; ++l) {

/*           IDAMAX call gives the relative index in column L of Aj where */
/*           max element is found. */
/*           Note: the first element in column L is in row K of */
/*                 matrix Aj. */

#line 251 "MB04TT.f"
	    ll = ifica1 + l;
#line 252 "MB04TT.f"
	    i__2 = mj - k + 1;
#line 252 "MB04TT.f"
	    bmx = (d__1 = a[idamax_(&i__2, &a[kk + ll * a_dim1], &c__1) + kk 
		    - 1 + ll * a_dim1], abs(d__1));
#line 253 "MB04TT.f"
	    if (bmx > bmxnrm) {
#line 254 "MB04TT.f"
		bmxnrm = bmx;
#line 255 "MB04TT.f"
		lsav = l;
#line 256 "MB04TT.f"
	    }
#line 257 "MB04TT.f"
/* L60: */
#line 257 "MB04TT.f"
	}

#line 259 "MB04TT.f"
	ll = ifica1 + k;
#line 260 "MB04TT.f"
	if (bmxnrm < *tol) {

/*           Set submatrix of Aj to zero. */

#line 264 "MB04TT.f"
	    i__1 = mj - k + 1;
#line 264 "MB04TT.f"
	    i__2 = nj - k + 1;
#line 264 "MB04TT.f"
	    dlaset_("Full", &i__1, &i__2, &c_b7, &c_b7, &a[kk + ll * a_dim1], 
		    lda, (ftnlen)4);
#line 266 "MB04TT.f"
	    lzero = TRUE_;
#line 267 "MB04TT.f"
	    *rank = k - 1;
#line 268 "MB04TT.f"
	} else {

/*           Check whether columns have to be interchanged. */

#line 272 "MB04TT.f"
	    if (lsav != k) {

/*              Interchange the columns in A which correspond to the */
/*              columns lsav and k in Aj. Store the permutation in IWORK. */

#line 277 "MB04TT.f"
		dswap_(m, &a[ll * a_dim1 + 1], &c__1, &a[(ifica1 + lsav) * 
			a_dim1 + 1], &c__1);
#line 278 "MB04TT.f"
		ip = iwork[lsav];
#line 279 "MB04TT.f"
		iwork[lsav] = iwork[k];
#line 280 "MB04TT.f"
		iwork[k] = ip;
#line 281 "MB04TT.f"
	    }

#line 283 "MB04TT.f"
	    ++k;
#line 284 "MB04TT.f"
	    mk1 = *n - ll + 1;

#line 286 "MB04TT.f"
	    i__1 = k;
#line 286 "MB04TT.f"
	    for (i__ = mj; i__ >= i__1; --i__) {

/*              II = absolute row number in A corresponding to row i in */
/*                   Aj. */

#line 291 "MB04TT.f"
		ii = ifira1 + i__;

/*              Construct Givens transformation to annihilate Aj(i,k). */
/*              Apply the row transformation to whole matrix A */
/*              (NOT only to Aj). */
/*              Update row transformation matrix Q, if needed. */

#line 298 "MB04TT.f"
		drotg_(&a[ii - 1 + ll * a_dim1], &a[ii + ll * a_dim1], &sc, &
			ss);
#line 299 "MB04TT.f"
		i__2 = mk1 - 1;
#line 299 "MB04TT.f"
		drot_(&i__2, &a[ii - 1 + (ll + 1) * a_dim1], lda, &a[ii + (ll 
			+ 1) * a_dim1], lda, &sc, &ss);
#line 301 "MB04TT.f"
		a[ii + ll * a_dim1] = 0.;
#line 302 "MB04TT.f"
		if (*updatq) {
#line 302 "MB04TT.f"
		    drot_(m, &q[(ii - 1) * q_dim1 + 1], &c__1, &q[ii * q_dim1 
			    + 1], &c__1, &sc, &ss);
#line 302 "MB04TT.f"
		}

/*              Determine boundary type of matrix E at rows II-1 and II. */

#line 307 "MB04TT.f"
		ist1 = istair[ii - 1];
#line 308 "MB04TT.f"
		ist2 = istair[ii];
#line 309 "MB04TT.f"
		if (ist1 * ist2 > 0) {
#line 310 "MB04TT.f"
		    if (ist1 > 0) {

/*                    boundary form = (* x) */
/*                                    (0 *) */

#line 315 "MB04TT.f"
			itype = 1;
#line 316 "MB04TT.f"
		    } else {

/*                    boundary form = (x x) */
/*                                    (x x) */

#line 321 "MB04TT.f"
			itype = 3;
#line 322 "MB04TT.f"
		    }
#line 323 "MB04TT.f"
		} else {
#line 324 "MB04TT.f"
		    if (ist1 < 0) {

/*                    boundary form = (x x) */
/*                                    (* x) */

#line 329 "MB04TT.f"
			itype = 2;
#line 330 "MB04TT.f"
		    } else {

/*                    boundary form = (* x) */
/*                                    (0 x) */

#line 335 "MB04TT.f"
			itype = 4;
#line 336 "MB04TT.f"
		    }
#line 337 "MB04TT.f"
		}

/*              Apply row transformation also to matrix E. */

/*              JC1 = absolute number of the column in E in which stair */
/*                    element of row i-1 of Ej is present. */
/*              JC2 = absolute number of the column in E in which stair */
/*                    element of row i of Ej is present. */

/*              Note: JC1 < JC2   if ITYPE = 1. */
/*                    JC1 = JC2   if ITYPE = 2, 3 or 4. */

#line 349 "MB04TT.f"
		jc1 = abs(ist1);
#line 350 "MB04TT.f"
		jc2 = abs(ist2);
#line 351 "MB04TT.f"
		jpvt = min(jc1,jc2);

#line 353 "MB04TT.f"
		i__2 = *n - jpvt + 1;
#line 353 "MB04TT.f"
		drot_(&i__2, &e[ii - 1 + jpvt * e_dim1], lde, &e[ii + jpvt * 
			e_dim1], lde, &sc, &ss);
#line 355 "MB04TT.f"
		eijpvt = e[ii + jpvt * e_dim1];

#line 357 "MB04TT.f"
		if (itype == 1) {

/*                 Construct column Givens transformation to annihilate */
/*                 E(ii,jpvt). */
/*                 Apply column Givens transformation to matrix E */
/*                 (NOT only to Ej). */

#line 364 "MB04TT.f"
		    drotg_(&e[ii + (jpvt + 1) * e_dim1], &e[ii + jpvt * 
			    e_dim1], &sc, &ss);
#line 365 "MB04TT.f"
		    i__2 = ii - 1;
#line 365 "MB04TT.f"
		    drot_(&i__2, &e[(jpvt + 1) * e_dim1 + 1], &c__1, &e[jpvt *
			     e_dim1 + 1], &c__1, &sc, &ss);
#line 367 "MB04TT.f"
		    e[ii + jpvt * e_dim1] = 0.;

/*                 Apply this transformation also to matrix A */
/*                 (NOT only to Aj). */
/*                 Update column transformation matrix Z, if needed. */

#line 373 "MB04TT.f"
		    drot_(m, &a[(jpvt + 1) * a_dim1 + 1], &c__1, &a[jpvt * 
			    a_dim1 + 1], &c__1, &sc, &ss);
#line 374 "MB04TT.f"
		    if (*updatz) {
#line 374 "MB04TT.f"
			drot_(n, &z__[(jpvt + 1) * z_dim1 + 1], &c__1, &z__[
				jpvt * z_dim1 + 1], &c__1, &sc, &ss);
#line 374 "MB04TT.f"
		    }

#line 377 "MB04TT.f"
		} else if (itype == 2) {
#line 378 "MB04TT.f"
		    if (abs(eijpvt) < *tol) {

/*                                                        (x x)    (* x) */
/*                    Boundary form has been changed from (* x) to (0 x). */

#line 383 "MB04TT.f"
			istpvt = istair[ii];
#line 384 "MB04TT.f"
			istair[ii - 1] = istpvt;
#line 385 "MB04TT.f"
			istair[ii] = -(istpvt + 1);
#line 386 "MB04TT.f"
			e[ii + jpvt * e_dim1] = 0.;
#line 387 "MB04TT.f"
		    }

#line 389 "MB04TT.f"
		} else if (itype == 4) {
#line 390 "MB04TT.f"
		    if (abs(eijpvt) >= *tol) {

/*                                                        (* x)    (x x) */
/*                    Boundary form has been changed from (0 x) to (* x). */

#line 395 "MB04TT.f"
			istpvt = istair[ii - 1];
#line 396 "MB04TT.f"
			istair[ii - 1] = -istpvt;
#line 397 "MB04TT.f"
			istair[ii] = istpvt;
#line 398 "MB04TT.f"
		    }
#line 399 "MB04TT.f"
		}
#line 400 "MB04TT.f"
/* L80: */
#line 400 "MB04TT.f"
	    }

#line 402 "MB04TT.f"
	}
#line 403 "MB04TT.f"
	goto L40;
#line 404 "MB04TT.f"
    }
/*     END WHILE 40 */

/*     Permute columns of Aj to original order. */

#line 409 "MB04TT.f"
    i__1 = ifira1 + *rank;
#line 409 "MB04TT.f"
    dlapmt_(&c_false, &i__1, &nj, &a[*ifica * a_dim1 + 1], lda, &iwork[1]);

#line 411 "MB04TT.f"
    return 0;
/* *** Last line of MB04TT *** */
} /* mb04tt_ */

