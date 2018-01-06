#line 1 "TB01PD.f"
/* TB01PD.f -- translated by f2c (version 20100827).
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

#line 1 "TB01PD.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int tb01pd_(char *job, char *equil, integer *n, integer *m, 
	integer *p, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, integer *nr, doublereal *tol, integer *
	iwork, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	job_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, kl, iz, itau;
    extern /* Subroutine */ int ab07md_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen), tb01id_(
	    char *, integer *, integer *, integer *, doublereal *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), tb01ud_(char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, integer *, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01xd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen);
    static integer ncont, maxmp, jwork;
    static logical lnjobc;
    static integer indcon;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical lnjobo, lequil;
    static integer wrkopt;


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

/*     To find a reduced (controllable, observable, or minimal) state- */
/*     space representation (Ar,Br,Cr) for any original state-space */
/*     representation (A,B,C). The matrix Ar is in upper block */
/*     Hessenberg form. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates whether the user wishes to remove the */
/*             uncontrollable and/or unobservable parts as follows: */
/*             = 'M':  Remove both the uncontrollable and unobservable */
/*                     parts to get a minimal state-space representation; */
/*             = 'C':  Remove the uncontrollable part only to get a */
/*                     controllable state-space representation; */
/*             = 'O':  Remove the unobservable part only to get an */
/*                     observable state-space representation. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily balance */
/*             the triplet (A,B,C) as follows: */
/*             = 'S':  Perform balancing (scaling); */
/*             = 'N':  Do not perform balancing. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, i.e. */
/*             the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.   P >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the original state dynamics matrix A. */
/*             On exit, the leading NR-by-NR part of this array contains */
/*             the upper block Hessenberg state dynamics matrix Ar of a */
/*             minimal, controllable, or observable realization for the */
/*             original system, depending on the value of JOB, JOB = 'M', */
/*             JOB = 'C', or JOB = 'O', respectively. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M), */
/*             if JOB = 'C', or (LDB,MAX(M,P)), otherwise. */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B; if JOB = 'M', */
/*             or JOB = 'O', the remainder of the leading N-by-MAX(M,P) */
/*             part is used as internal workspace. */
/*             On exit, the leading NR-by-M part of this array contains */
/*             the transformed input/state matrix Br of a minimal, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'M', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             If JOB = 'C', only the first IWORK(1) rows of B are */
/*             nonzero. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C; if JOB = 'M', */
/*             or JOB = 'O', the remainder of the leading MAX(M,P)-by-N */
/*             part is used as internal workspace. */
/*             On exit, the leading P-by-NR part of this array contains */
/*             the transformed state/output matrix Cr of a minimal, */
/*             controllable, or observable realization for the original */
/*             system, depending on the value of JOB, JOB = 'M', */
/*             JOB = 'C', or JOB = 'O', respectively. */
/*             If JOB = 'M', or JOB = 'O', only the last IWORK(1) columns */
/*             (in the first NR columns) of C are nonzero. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= MAX(1,M,P) if N > 0. */
/*             LDC >= 1          if N = 0. */

/*     NR      (output) INTEGER */
/*             The order of the reduced state-space representation */
/*             (Ar,Br,Cr) of a minimal, controllable, or observable */
/*             realization for the original system, depending on */
/*             JOB = 'M', JOB = 'C', or JOB = 'O'. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used in rank determination when */
/*             transforming (A, B, C). If the user sets TOL > 0, then */
/*             the given value of TOL is used as a lower bound for the */
/*             reciprocal condition number (see the description of the */
/*             argument RCOND in the SLICOT routine MB03OD);  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance */
/*             (determined by the SLICOT routine TB01UD) is used instead. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N+MAX(M,P)) */
/*             On exit, if INFO = 0, the first nonzero elements of */
/*             IWORK(1:N) return the orders of the diagonal blocks of A. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N + MAX(N, 3*M, 3*P)). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     If JOB = 'M', the matrices A and B are operated on by orthogonal */
/*     similarity transformations (made up of products of Householder */
/*     transformations) so as to produce an upper block Hessenberg matrix */
/*     A1 and a matrix B1 with all but its first rank(B) rows zero; this */
/*     separates out the controllable part of the original system. */
/*     Applying the same algorithm to the dual of this subsystem, */
/*     therefore separates out the controllable and observable (i.e. */
/*     minimal) part of the original system representation, with the */
/*     final Ar upper block Hessenberg (after using pertransposition). */
/*     If JOB = 'C', or JOB = 'O', only the corresponding part of the */
/*     above procedure is applied. */

/*     REFERENCES */

/*     [1] Van Dooren, P. */
/*         The Generalized Eigenstructure Problem in Linear System */
/*         Theory. (Algorithm 1) */
/*         IEEE Trans. Auto. Contr., AC-26, pp. 111-129, 1981. */

/*     NUMERICAL ASPECTS */
/*                               3 */
/*     The algorithm requires 0(N ) operations and is backward stable. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 1998. */

/*     REVISIONS */

/*     A. Varga, DLR Oberpfaffenhofen, July 1998. */
/*     A. Varga, DLR Oberpfaffenhofen, April 28, 1999. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2004. */

/*     KEYWORDS */

/*     Hessenberg form, minimal realization, multivariable system, */
/*     orthogonal transformation, state-space model, state-space */
/*     representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 218 "TB01PD.f"
    /* Parameter adjustments */
#line 218 "TB01PD.f"
    a_dim1 = *lda;
#line 218 "TB01PD.f"
    a_offset = 1 + a_dim1;
#line 218 "TB01PD.f"
    a -= a_offset;
#line 218 "TB01PD.f"
    b_dim1 = *ldb;
#line 218 "TB01PD.f"
    b_offset = 1 + b_dim1;
#line 218 "TB01PD.f"
    b -= b_offset;
#line 218 "TB01PD.f"
    c_dim1 = *ldc;
#line 218 "TB01PD.f"
    c_offset = 1 + c_dim1;
#line 218 "TB01PD.f"
    c__ -= c_offset;
#line 218 "TB01PD.f"
    --iwork;
#line 218 "TB01PD.f"
    --dwork;
#line 218 "TB01PD.f"

#line 218 "TB01PD.f"
    /* Function Body */
#line 218 "TB01PD.f"
    *info = 0;
#line 219 "TB01PD.f"
    maxmp = max(*m,*p);
#line 220 "TB01PD.f"
    lnjobc = ! lsame_(job, "C", (ftnlen)1, (ftnlen)1);
#line 221 "TB01PD.f"
    lnjobo = ! lsame_(job, "O", (ftnlen)1, (ftnlen)1);
#line 222 "TB01PD.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 226 "TB01PD.f"
    if (lnjobc && lnjobo && ! lsame_(job, "M", (ftnlen)1, (ftnlen)1)) {
#line 227 "TB01PD.f"
	*info = -1;
#line 228 "TB01PD.f"
    } else if (! lequil && ! lsame_(equil, "N", (ftnlen)1, (ftnlen)1)) {
#line 229 "TB01PD.f"
	*info = -2;
#line 230 "TB01PD.f"
    } else if (*n < 0) {
#line 231 "TB01PD.f"
	*info = -3;
#line 232 "TB01PD.f"
    } else if (*m < 0) {
#line 233 "TB01PD.f"
	*info = -4;
#line 234 "TB01PD.f"
    } else if (*p < 0) {
#line 235 "TB01PD.f"
	*info = -5;
#line 236 "TB01PD.f"
    } else if (*lda < max(1,*n)) {
#line 237 "TB01PD.f"
	*info = -7;
#line 238 "TB01PD.f"
    } else if (*ldb < max(1,*n)) {
#line 239 "TB01PD.f"
	*info = -9;
#line 240 "TB01PD.f"
    } else if (*ldc < 1 || *n > 0 && *ldc < maxmp) {
#line 241 "TB01PD.f"
	*info = -11;
#line 242 "TB01PD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 242 "TB01PD.f"
	i__3 = *n, i__4 = maxmp * 3;
#line 242 "TB01PD.f"
	i__1 = 1, i__2 = *n + max(i__3,i__4);
#line 242 "TB01PD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 243 "TB01PD.f"
	    *info = -16;
#line 244 "TB01PD.f"
	}
#line 244 "TB01PD.f"
    }

#line 246 "TB01PD.f"
    if (*info != 0) {

/*        Error return. */

#line 250 "TB01PD.f"
	i__1 = -(*info);
#line 250 "TB01PD.f"
	xerbla_("TB01PD", &i__1, (ftnlen)6);
#line 251 "TB01PD.f"
	return 0;
#line 252 "TB01PD.f"
    }

/*     Quick return if possible. */

#line 256 "TB01PD.f"
    if (*n == 0 || lnjobc && min(*n,*p) == 0 || lnjobo && min(*n,*m) == 0) {
#line 258 "TB01PD.f"
	*nr = 0;

#line 260 "TB01PD.f"
	i__1 = *n;
#line 260 "TB01PD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 261 "TB01PD.f"
	    iwork[i__] = 0;
#line 262 "TB01PD.f"
/* L5: */
#line 262 "TB01PD.f"
	}

#line 264 "TB01PD.f"
	dwork[1] = 1.;
#line 265 "TB01PD.f"
	return 0;
#line 266 "TB01PD.f"
    }

/*     If required, balance the triplet (A,B,C) (default MAXRED). */
/*     Workspace: need N. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the code, */
/*     as well as the preferred amount for good performance.) */

#line 275 "TB01PD.f"
    if (lequil) {
#line 276 "TB01PD.f"
	maxred = 0.;
#line 277 "TB01PD.f"
	tb01id_("A", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb, &
		c__[c_offset], ldc, &dwork[1], info, (ftnlen)1);
#line 279 "TB01PD.f"
	wrkopt = *n;
#line 280 "TB01PD.f"
    } else {
#line 281 "TB01PD.f"
	wrkopt = 1;
#line 282 "TB01PD.f"
    }

#line 284 "TB01PD.f"
    iz = 1;
#line 285 "TB01PD.f"
    itau = 1;
#line 286 "TB01PD.f"
    jwork = itau + *n;
#line 287 "TB01PD.f"
    if (lnjobo) {

/*        Separate out controllable subsystem (of order NCONT): */
/*        A <-- Z'*A*Z,  B <-- Z'*B,  C <-- C*Z. */

/*        Workspace: need   N + MAX(N, 3*M, P). */
/*                   prefer larger. */

#line 295 "TB01PD.f"
	i__1 = *ldwork - jwork + 1;
#line 295 "TB01PD.f"
	tb01ud_("No Z", n, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[
		c_offset], ldc, &ncont, &indcon, &iwork[1], &dwork[iz], &c__1,
		 &dwork[itau], tol, &iwork[*n + 1], &dwork[jwork], &i__1, 
		info, (ftnlen)4);

#line 299 "TB01PD.f"
	wrkopt = (integer) dwork[jwork] + jwork - 1;
#line 300 "TB01PD.f"
    } else {
#line 301 "TB01PD.f"
	ncont = *n;
#line 302 "TB01PD.f"
    }

#line 304 "TB01PD.f"
    if (lnjobc) {

/*        Separate out the observable subsystem (of order NR): */
/*        Form the dual of the subsystem of order NCONT (which is */
/*        controllable, if JOB = 'M'), leaving rest as it is. */

#line 310 "TB01PD.f"
	ab07md_("Z", &ncont, m, p, &a[a_offset], lda, &b[b_offset], ldb, &c__[
		c_offset], ldc, &dwork[1], &c__1, info, (ftnlen)1);

/*        And separate out the controllable part of this dual subsystem. */

/*        Workspace: need   NCONT + MAX(NCONT, 3*P, M). */
/*                   prefer larger. */

#line 318 "TB01PD.f"
	i__1 = *ldwork - jwork + 1;
#line 318 "TB01PD.f"
	tb01ud_("No Z", &ncont, p, m, &a[a_offset], lda, &b[b_offset], ldb, &
		c__[c_offset], ldc, nr, &indcon, &iwork[1], &dwork[iz], &c__1,
		 &dwork[itau], tol, &iwork[*n + 1], &dwork[jwork], &i__1, 
		info, (ftnlen)4);

/* Computing MAX */
#line 322 "TB01PD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 322 "TB01PD.f"
	wrkopt = max(i__1,i__2);

/*        Transpose and reorder (to get a block upper Hessenberg */
/*        matrix A), giving, for JOB = 'M', the controllable and */
/*        observable (i.e., minimal) part of original system. */

#line 328 "TB01PD.f"
	if (indcon > 0) {
#line 329 "TB01PD.f"
	    kl = iwork[1] - 1;
#line 330 "TB01PD.f"
	    if (indcon >= 2) {
#line 330 "TB01PD.f"
		kl += iwork[2];
#line 330 "TB01PD.f"
	    }
#line 332 "TB01PD.f"
	} else {
#line 333 "TB01PD.f"
	    kl = 0;
#line 334 "TB01PD.f"
	}
/* Computing MAX */
#line 335 "TB01PD.f"
	i__2 = 0, i__3 = *nr - 1;
#line 335 "TB01PD.f"
	i__1 = max(i__2,i__3);
#line 335 "TB01PD.f"
	tb01xd_("Zero D", nr, p, m, &kl, &i__1, &a[a_offset], lda, &b[
		b_offset], ldb, &c__[c_offset], ldc, &dwork[1], &c__1, info, (
		ftnlen)6);
#line 337 "TB01PD.f"
    } else {
#line 338 "TB01PD.f"
	*nr = ncont;
#line 339 "TB01PD.f"
    }

/*     Annihilate the trailing components of IWORK(1:N). */

#line 343 "TB01PD.f"
    i__1 = *n;
#line 343 "TB01PD.f"
    for (i__ = indcon + 1; i__ <= i__1; ++i__) {
#line 344 "TB01PD.f"
	iwork[i__] = 0;
#line 345 "TB01PD.f"
/* L10: */
#line 345 "TB01PD.f"
    }

/*     Set optimal workspace dimension. */

#line 349 "TB01PD.f"
    dwork[1] = (doublereal) wrkopt;
#line 350 "TB01PD.f"
    return 0;
/* *** Last line of TB01PD *** */
} /* tb01pd_ */

