#line 1 "IB01ND.f"
/* IB01ND.f -- translated by f2c (version 20100827).
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

#line 1 "IB01ND.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b37 = .66666666666666663;
static integer c__0 = 0;
static doublereal c_b50 = 0.;

/* Subroutine */ int ib01nd_(char *meth, char *jobd, integer *nobr, integer *
	m, integer *l, doublereal *r__, integer *ldr, doublereal *sv, 
	doublereal *tol, integer *iwork, doublereal *dwork, integer *ldwork, 
	integer *iwarn, integer *info, ftnlen meth_len, ftnlen jobd_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, i__1, i__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, nr, nr2, nr3, nr4;
    static doublereal dum[1], eps;
    static integer rank, ierr, itau;
    static doublereal sval[3], toll;
    static integer rank1;
    static logical n4sid;
    static integer itau2, itau3;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    mb04id_(integer *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *), mb03od_(char *, integer *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen), mb04od_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), mb03ud_(char *, char *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     integer *, ftnlen, ftnlen);
    static logical jobdm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb04iy_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static integer lnobr, mnobr;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical moesp;
    static integer jwork;
    static doublereal rcond1, rcond2;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer llmnob, lmmnob, llnobr, lmnobr, mmnobr;
    extern /* Subroutine */ int dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal thresh;
    static integer nrsave;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal svlmax;


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

/*     To find the singular value decomposition (SVD) giving the system */
/*     order, using the triangular factor of the concatenated block */
/*     Hankel matrices. Related preliminary calculations needed for */
/*     computing the system matrices are also performed. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     METH    CHARACTER*1 */
/*             Specifies the subspace identification method to be used, */
/*             as follows: */
/*             = 'M':  MOESP  algorithm with past inputs and outputs; */
/*             = 'N':  N4SID  algorithm. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not the matrices B and D should later */
/*             be computed using the MOESP approach, as follows: */
/*             = 'M':  the matrices B and D should later be computed */
/*                     using the MOESP approach; */
/*             = 'N':  the matrices B and D should not be computed using */
/*                     the MOESP approach. */
/*             This parameter is not relevant for METH = 'N'. */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             The number of block rows,  s,  in the input and output */
/*             block Hankel matrices.  NOBR > 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     R       (input/output) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             On entry, the leading 2*(M+L)*NOBR-by-2*(M+L)*NOBR upper */
/*             triangular part of this array must contain the upper */
/*             triangular factor R from the QR factorization of the */
/*             concatenated block Hankel matrices. Denote  R_ij, */
/*             i,j = 1:4,  the ij submatrix of  R,  partitioned by */
/*             M*NOBR,  M*NOBR,  L*NOBR,  and  L*NOBR  rows and columns. */
/*             On exit, if INFO = 0, the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of this */
/*             array contains the matrix S, the processed upper */
/*             triangular factor R, as required by other subroutines. */
/*             Specifically, let  S_ij, i,j = 1:4,  be the ij submatrix */
/*             of  S,  partitioned by  M*NOBR,  L*NOBR,  M*NOBR,  and */
/*             L*NOBR  rows and columns. The submatrix  S_22  contains */
/*             the matrix of left singular vectors needed subsequently. */
/*             Useful information is stored in  S_11  and in the */
/*             block-column  S_14 : S_44.  For METH = 'M' and JOBD = 'M', */
/*             the upper triangular part of  S_31  contains the upper */
/*             triangular factor in the QR factorization of the matrix */
/*             R_1c = [ R_12'  R_22'  R_11' ]',  and  S_12  contains the */
/*             corresponding leading part of the transformed matrix */
/*             R_2c = [ R_13'  R_23'  R_14' ]'.  For  METH = 'N',  the */
/*             subarray  S_41 : S_43  contains the transpose of the */
/*             matrix contained in  S_14 : S_34. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= MAX( 2*(M+L)*NOBR, 3*M*NOBR ), */
/*                                  for METH = 'M' and JOBD = 'M'; */
/*             LDR >= 2*(M+L)*NOBR, for METH = 'M' and JOBD = 'N' or */
/*                                  for METH = 'N'. */

/*     SV      (output) DOUBLE PRECISION array, dimension ( L*NOBR ) */
/*             The singular values of the relevant part of the triangular */
/*             factor from the QR factorization of the concatenated block */
/*             Hankel matrices. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used for estimating the rank of */
/*             matrices. If the user sets  TOL > 0,  then the given value */
/*             of  TOL  is used as a lower bound for the reciprocal */
/*             condition number;  an m-by-n matrix whose estimated */
/*             condition number is less than  1/TOL  is considered to */
/*             be of full rank.  If the user sets  TOL <= 0,  then an */
/*             implicitly computed, default tolerance, defined by */
/*             TOLDEF = m*n*EPS,  is used instead, where  EPS  is the */
/*             relative machine precision (see LAPACK Library routine */
/*             DLAMCH). */
/*             This parameter is not used for  METH = 'M'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension ((M+L)*NOBR) */
/*             This parameter is not referenced for METH = 'M'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK,  and, for  METH = 'N',  DWORK(2)  and  DWORK(3) */
/*             contain the reciprocal condition numbers of the */
/*             triangular factors of the matrices  U_f  and  r_1  [6]. */
/*             On exit, if  INFO = -12,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( (2*M-1)*NOBR, (M+L)*NOBR, 5*L*NOBR ), */
/*                                         if METH = 'M' and JOBD = 'M'; */
/*             LDWORK >=  5*L*NOBR,        if METH = 'M' and JOBD = 'N'; */
/*             LDWORK >=  5*(M+L)*NOBR+1,  if METH = 'N'. */
/*             For good performance,  LDWORK  should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 4:  the least squares problems with coefficient matrix */
/*                   U_f,  used for computing the weighted oblique */
/*                   projection (for METH = 'N'), have a rank-deficient */
/*                   coefficient matrix; */
/*             = 5:  the least squares problem with coefficient matrix */
/*                   r_1  [6], used for computing the weighted oblique */
/*                   projection (for METH = 'N'), has a rank-deficient */
/*                   coefficient matrix. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  the singular value decomposition (SVD) algorithm did */
/*                   not converge. */

/*     METHOD */

/*     A singular value decomposition (SVD) of a certain matrix is */
/*     computed, which reveals the order  n  of the system as the number */
/*     of "non-zero" singular values. For the MOESP approach, this matrix */
/*     is  [ R_24'  R_34' ]' := R(ms+1:(2m+l)s,(2m+l)s+1:2(m+l)s), */
/*     where  R  is the upper triangular factor  R  constructed by SLICOT */
/*     Library routine  IB01MD.  For the N4SID approach, a weighted */
/*     oblique projection is computed from the upper triangular factor  R */
/*     and its SVD is then found. */

/*     REFERENCES */

/*     [1] Verhaegen M., and Dewilde, P. */
/*         Subspace Model Identification. Part 1: The output-error */
/*         state-space model identification class of algorithms. */
/*         Int. J. Control, 56, pp. 1187-1210, 1992. */

/*     [2] Verhaegen M. */
/*         Subspace Model Identification. Part 3: Analysis of the */
/*         ordinary output-error state-space model identification */
/*         algorithm. */
/*         Int. J. Control, 58, pp. 555-586, 1993. */

/*     [3] Verhaegen M. */
/*         Identification of the deterministic part of MIMO state space */
/*         models given in innovations form from input-output data. */
/*         Automatica, Vol.30, No.1, pp.61-74, 1994. */

/*     [4] Van Overschee, P., and De Moor, B. */
/*         N4SID: Subspace Algorithms for the Identification of */
/*         Combined Deterministic-Stochastic Systems. */
/*         Automatica, Vol.30, No.1, pp. 75-93, 1994. */

/*     [5] Van Overschee, P., and De Moor, B. */
/*         Subspace Identification for Linear Systems: Theory - */
/*         Implementation - Applications. */
/*         Kluwer Academic Publishers, Boston/London/Dordrecht, 1996. */

/*     [6] Sima, V. */
/*         Subspace-based Algorithms for Multivariable System */
/*         Identification. */
/*         Studies in Informatics and Control, 5, pp. 335-344, 1996. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */
/*                                      3 */
/*     The algorithm requires 0(((m+l)s) ) floating point operations. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999. */

/*     REVISIONS */

/*     Feb. 2000, Feb. 2001, Feb. 2004, March 2005. */

/*     KEYWORDS */

/*     Identification methods, multivariable systems, QR decomposition, */
/*     singular value decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 257 "IB01ND.f"
    /* Parameter adjustments */
#line 257 "IB01ND.f"
    r_dim1 = *ldr;
#line 257 "IB01ND.f"
    r_offset = 1 + r_dim1;
#line 257 "IB01ND.f"
    r__ -= r_offset;
#line 257 "IB01ND.f"
    --sv;
#line 257 "IB01ND.f"
    --iwork;
#line 257 "IB01ND.f"
    --dwork;
#line 257 "IB01ND.f"

#line 257 "IB01ND.f"
    /* Function Body */
#line 257 "IB01ND.f"
    moesp = lsame_(meth, "M", (ftnlen)1, (ftnlen)1);
#line 258 "IB01ND.f"
    n4sid = lsame_(meth, "N", (ftnlen)1, (ftnlen)1);
#line 259 "IB01ND.f"
    jobdm = lsame_(jobd, "M", (ftnlen)1, (ftnlen)1);
#line 260 "IB01ND.f"
    mnobr = *m * *nobr;
#line 261 "IB01ND.f"
    lnobr = *l * *nobr;
#line 262 "IB01ND.f"
    llnobr = lnobr + lnobr;
#line 263 "IB01ND.f"
    lmnobr = lnobr + mnobr;
#line 264 "IB01ND.f"
    mmnobr = mnobr + mnobr;
#line 265 "IB01ND.f"
    lmmnob = mmnobr + lnobr;
#line 266 "IB01ND.f"
    nr = lmnobr + lmnobr;
#line 267 "IB01ND.f"
    *iwarn = 0;
#line 268 "IB01ND.f"
    *info = 0;

/*     Check the scalar input parameters. */

#line 272 "IB01ND.f"
    if (! (moesp || n4sid)) {
#line 273 "IB01ND.f"
	*info = -1;
#line 274 "IB01ND.f"
    } else if (moesp && ! (jobdm || lsame_(jobd, "N", (ftnlen)1, (ftnlen)1))) 
	    {
#line 275 "IB01ND.f"
	*info = -2;
#line 276 "IB01ND.f"
    } else if (*nobr <= 0) {
#line 277 "IB01ND.f"
	*info = -3;
#line 278 "IB01ND.f"
    } else if (*m < 0) {
#line 279 "IB01ND.f"
	*info = -4;
#line 280 "IB01ND.f"
    } else if (*l <= 0) {
#line 281 "IB01ND.f"
	*info = -5;
#line 282 "IB01ND.f"
    } else if (*ldr < nr || moesp && jobdm && *ldr < mnobr * 3) {
#line 284 "IB01ND.f"
	*info = -7;
#line 285 "IB01ND.f"
    } else {

/*        Compute workspace. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*         minimal amount of workspace needed at that point in the code, */
/*         as well as the preferred amount for good performance. */
/*         NB refers to the optimal block size for the immediately */
/*         following subroutine, as returned by ILAENV.) */

#line 294 "IB01ND.f"
	minwrk = 1;
#line 295 "IB01ND.f"
	if (*ldwork >= 1) {
#line 296 "IB01ND.f"
	    if (moesp) {
#line 297 "IB01ND.f"
		minwrk = lnobr * 5;
#line 298 "IB01ND.f"
		if (jobdm) {
/* Computing MAX */
#line 298 "IB01ND.f"
		    i__1 = mmnobr - *nobr, i__1 = max(i__1,lmnobr);
#line 298 "IB01ND.f"
		    minwrk = max(i__1,minwrk);
#line 298 "IB01ND.f"
		}
#line 300 "IB01ND.f"
		maxwrk = lnobr + lnobr * ilaenv_(&c__1, "DGEQRF", " ", &
			lmnobr, &lnobr, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 302 "IB01ND.f"
	    } else {

/* Computing MAX */
#line 304 "IB01ND.f"
		i__1 = minwrk, i__2 = lmnobr * 5 + 1;
#line 304 "IB01ND.f"
		minwrk = max(i__1,i__2);
/* Computing MAX */
#line 305 "IB01ND.f"
		i__1 = mnobr + mnobr * ilaenv_(&c__1, "DGEQRF", " ", &mmnobr, 
			&mnobr, &c_n1, &c_n1, (ftnlen)6, (ftnlen)1), i__2 = 
			mnobr + llnobr * ilaenv_(&c__1, "DORMQR", "LT", &
			mmnobr, &llnobr, &mnobr, &c_n1, (ftnlen)6, (ftnlen)2);
#line 305 "IB01ND.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 309 "IB01ND.f"
		i__1 = maxwrk, i__2 = mnobr + lnobr * ilaenv_(&c__1, "DORMQR",
			 "LN", &mmnobr, &lnobr, &mnobr, &c_n1, (ftnlen)6, (
			ftnlen)2);
#line 309 "IB01ND.f"
		maxwrk = max(i__1,i__2);
/* Computing MAX */
#line 312 "IB01ND.f"
		i__1 = maxwrk, i__2 = lnobr + lnobr * ilaenv_(&c__1, "DGEQRF",
			 " ", &lmmnob, &lnobr, &c_n1, &c_n1, (ftnlen)6, (
			ftnlen)1);
#line 312 "IB01ND.f"
		maxwrk = max(i__1,i__2);
#line 314 "IB01ND.f"
	    }
#line 315 "IB01ND.f"
	    maxwrk = max(minwrk,maxwrk);
#line 316 "IB01ND.f"
	}

#line 318 "IB01ND.f"
	if (*ldwork < minwrk) {
#line 319 "IB01ND.f"
	    *info = -12;
#line 320 "IB01ND.f"
	    dwork[1] = (doublereal) minwrk;
#line 321 "IB01ND.f"
	}
#line 322 "IB01ND.f"
    }

/*     Return if there are illegal arguments. */

#line 326 "IB01ND.f"
    if (*info != 0) {
#line 327 "IB01ND.f"
	i__1 = -(*info);
#line 327 "IB01ND.f"
	xerbla_("IB01ND", &i__1, (ftnlen)6);
#line 328 "IB01ND.f"
	return 0;
#line 329 "IB01ND.f"
    }

/*     Compute pointers to the needed blocks of  R. */

#line 333 "IB01ND.f"
    nr2 = mnobr + 1;
#line 334 "IB01ND.f"
    nr3 = mmnobr + 1;
#line 335 "IB01ND.f"
    nr4 = lmmnob + 1;
#line 336 "IB01ND.f"
    itau = 1;
#line 337 "IB01ND.f"
    jwork = itau + mnobr;

#line 339 "IB01ND.f"
    if (moesp) {

/*        MOESP approach. */

#line 343 "IB01ND.f"
	if (*m > 0 && jobdm) {

/*           Rearrange the blocks of  R: */
/*           Copy the (1,1) block into the position (3,2) and */
/*           copy the (1,4) block into (3,3). */

#line 349 "IB01ND.f"
	    dlacpy_("Upper", &mnobr, &mnobr, &r__[r_offset], ldr, &r__[nr3 + 
		    nr2 * r_dim1], ldr, (ftnlen)5);
#line 351 "IB01ND.f"
	    dlacpy_("Full", &mnobr, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &r__[
		    nr3 + nr3 * r_dim1], ldr, (ftnlen)4);

/*           Using structure, triangularize the matrix */
/*              R_1c = [ R_12'  R_22'  R_11' ]' */
/*           and then apply the transformations to the matrix */
/*              R_2c = [ R_13'  R_23'  R_14' ]'. */
/*           Workspace: need M*NOBR + MAX(M-1,L)*NOBR. */

#line 360 "IB01ND.f"
	    mb04od_("Upper", &mnobr, &lnobr, &mnobr, &r__[nr2 + nr2 * r_dim1],
		     ldr, &r__[nr3 + nr2 * r_dim1], ldr, &r__[nr2 + nr3 * 
		    r_dim1], ldr, &r__[nr3 + nr3 * r_dim1], ldr, &dwork[itau],
		     &dwork[jwork], (ftnlen)5);
#line 363 "IB01ND.f"
	    i__1 = mnobr - 1;
#line 363 "IB01ND.f"
	    i__2 = *ldwork - jwork + 1;
#line 363 "IB01ND.f"
	    mb04id_(&mmnobr, &mnobr, &i__1, &lnobr, &r__[nr2 * r_dim1 + 1], 
		    ldr, &r__[nr3 * r_dim1 + 1], ldr, &dwork[itau], &dwork[
		    jwork], &i__2, &ierr);
/* Computing MAX */
#line 366 "IB01ND.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 366 "IB01ND.f"
	    maxwrk = max(i__1,i__2);

/*           Copy the leading  M*NOBR x M*NOBR  and  M*NOBR x L*NOBR */
/*           submatrices of  R_1c  and  R_2c,  respectively, into their */
/*           final positions, required by SLICOT Library routine  IB01PD. */

#line 372 "IB01ND.f"
	    dlacpy_("Upper", &mnobr, &mnobr, &r__[nr2 * r_dim1 + 1], ldr, &
		    r__[lmnobr + 1 + r_dim1], ldr, (ftnlen)5);
#line 374 "IB01ND.f"
	    dlacpy_("Full", &mnobr, &lnobr, &r__[nr3 * r_dim1 + 1], ldr, &r__[
		    nr2 * r_dim1 + 1], ldr, (ftnlen)4);
#line 376 "IB01ND.f"
	}

/*        Copy [ R_24'  R_34' ]'  in  [ R_22'  R_32' ]'. */

#line 380 "IB01ND.f"
	dlacpy_("Full", &lmnobr, &lnobr, &r__[nr2 + nr4 * r_dim1], ldr, &r__[
		nr2 + nr2 * r_dim1], ldr, (ftnlen)4);

/*        Triangularize the matrix in  [ R_22'  R_32' ]'. */
/*        Workspace: need 2*L*NOBR; prefer L*NOBR + L*NOBR*NB. */

#line 386 "IB01ND.f"
	jwork = itau + lnobr;
#line 387 "IB01ND.f"
	i__1 = *ldwork - jwork + 1;
#line 387 "IB01ND.f"
	dgeqrf_(&lmnobr, &lnobr, &r__[nr2 + nr2 * r_dim1], ldr, &dwork[itau], 
		&dwork[jwork], &i__1, &ierr);

#line 390 "IB01ND.f"
    } else {

/*        N4SID approach. */

#line 394 "IB01ND.f"
	dum[0] = 0.;
#line 395 "IB01ND.f"
	llmnob = llnobr + mnobr;

/*        Set the precision parameters. A threshold value  EPS**(2/3)  is */
/*        used for deciding to use pivoting or not, where  EPS  is the */
/*        relative machine precision (see LAPACK Library routine DLAMCH). */

#line 401 "IB01ND.f"
	toll = *tol;
#line 402 "IB01ND.f"
	eps = dlamch_("Precision", (ftnlen)9);
#line 403 "IB01ND.f"
	thresh = pow_dd(&eps, &c_b37);

#line 405 "IB01ND.f"
	if (*m > 0) {

/*           For efficiency of later calculations, interchange the first */
/*           two block-columns. The corresponding submatrices are */
/*           redefined according to their new position. */

#line 411 "IB01ND.f"
	    i__1 = mnobr;
#line 411 "IB01ND.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 412 "IB01ND.f"
		dswap_(&i__, &r__[i__ * r_dim1 + 1], &c__1, &r__[(mnobr + i__)
			 * r_dim1 + 1], &c__1);
#line 413 "IB01ND.f"
		dcopy_(&mnobr, &r__[i__ + 1 + (mnobr + i__) * r_dim1], &c__1, 
			&r__[i__ + 1 + i__ * r_dim1], &c__1);
#line 414 "IB01ND.f"
		i__2 = mmnobr - i__;
#line 414 "IB01ND.f"
		dcopy_(&i__2, dum, &c__0, &r__[i__ + 1 + (mnobr + i__) * 
			r_dim1], &c__1);
#line 415 "IB01ND.f"
/* L10: */
#line 415 "IB01ND.f"
	    }

/*           Now, */

/*           U_f = [ R_11'  R_21'    0      0   ]', */
/*           U_p = [ R_12'    0      0      0   ]', */
/*           Y_p = [ R_13'  R_23'  R_33'    0   ]',  and */
/*           Y_f = [ R_14'  R_24'  R_34'  R_44' ]', */

/*           where  R_21,  R_12,  R_33,  and  R_44  are upper triangular. */
/*           Define  W_p := [ U_p  Y_p ]. */

/*           Prepare the computation of residuals of the two least */
/*           squares problems giving the weighted oblique projection P: */

/*           r_1 = W_p - U_f X_1,   X_1 = arg min || U_f X - W_p ||, */
/*           r_2 = Y_f - U_f X_2,   X_2 = arg min || U_f X - Y_f ||, */

/*           P = (arg min || r_1 X - r_2 ||)' r_1'.                   (1) */

/*           Alternately,  P'  is given by the projection */
/*              P' = Q_1 (Q_1)' r_2, */
/*           where  Q_1  contains the first  k  columns of the orthogonal */
/*           matrix in the  QR  factorization of  r_1,  k := rank(r_1). */

/*           Triangularize the matrix  U_f = q r  (using structure), and */
/*           apply the transformation  q'  to the corresponding part of */
/*           the matrices  W_p,  and  Y_f. */
/*           Workspace: need 2*(M+L)*NOBR. */

#line 445 "IB01ND.f"
	    i__1 = mnobr - 1;
#line 445 "IB01ND.f"
	    i__2 = *ldwork - jwork + 1;
#line 445 "IB01ND.f"
	    mb04id_(&mmnobr, &mnobr, &i__1, &llmnob, &r__[r_offset], ldr, &
		    r__[nr2 * r_dim1 + 1], ldr, &dwork[itau], &dwork[jwork], &
		    i__2, &ierr);
/* Computing MAX */
#line 448 "IB01ND.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 448 "IB01ND.f"
	    maxwrk = max(i__1,i__2);

/*           Save updated  Y_f  (transposed) in the last block-row of  R. */

#line 452 "IB01ND.f"
	    ma02ad_("Full", &lmmnob, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &
		    r__[nr4 + r_dim1], ldr, (ftnlen)4);

/*           Check the condition of the triangular factor  r  and decide */
/*           to use pivoting or not. */
/*           Workspace: need 4*M*NOBR. */

#line 459 "IB01ND.f"
	    dtrcon_("1-norm", "Upper", "NonUnit", &mnobr, &r__[r_offset], ldr,
		     &rcond1, &dwork[jwork], &iwork[1], &ierr, (ftnlen)6, (
		    ftnlen)5, (ftnlen)7);

#line 462 "IB01ND.f"
	    if (toll <= 0.) {
#line 462 "IB01ND.f"
		toll = mnobr * mnobr * eps;
#line 462 "IB01ND.f"
	    }
#line 464 "IB01ND.f"
	    if (rcond1 > max(toll,thresh)) {

/*              U_f is considered full rank and no pivoting is used. */

#line 468 "IB01ND.f"
		dlaset_("Full", &mnobr, &llmnob, &c_b50, &c_b50, &r__[nr2 * 
			r_dim1 + 1], ldr, (ftnlen)4);
#line 470 "IB01ND.f"
	    } else {

/*              Save information about  q  in the (2,1) block of  R. */
/*              Use QR factorization with column pivoting,  r P = Q R. */
/*              Information on  Q  is stored in the strict lower triangle */
/*              of R_11  and in  DWORK(ITAU2). */

#line 477 "IB01ND.f"
		i__1 = mnobr - 1;
#line 477 "IB01ND.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 478 "IB01ND.f"
		    i__2 = nr2;
#line 478 "IB01ND.f"
		    for (j = mmnobr; j >= i__2; --j) {
#line 479 "IB01ND.f"
			r__[j + i__ * r_dim1] = r__[j - mnobr + i__ + i__ * 
				r_dim1];
#line 480 "IB01ND.f"
/* L15: */
#line 480 "IB01ND.f"
		    }
#line 481 "IB01ND.f"
		    i__2 = mnobr - i__;
#line 481 "IB01ND.f"
		    dcopy_(&i__2, dum, &c__0, &r__[i__ + 1 + i__ * r_dim1], &
			    c__1);
#line 482 "IB01ND.f"
		    iwork[i__] = 0;
#line 483 "IB01ND.f"
/* L20: */
#line 483 "IB01ND.f"
		}

#line 485 "IB01ND.f"
		iwork[mnobr] = 0;

/*              Workspace: need   5*M*NOBR+1. */
/*                         prefer 4*M*NOBR + (M*NOBR+1)*NB. */

#line 490 "IB01ND.f"
		itau2 = jwork;
#line 491 "IB01ND.f"
		jwork = itau2 + mnobr;
#line 492 "IB01ND.f"
		svlmax = 0.;
#line 493 "IB01ND.f"
		i__1 = *ldwork - jwork + 1;
#line 493 "IB01ND.f"
		mb03od_("QR", &mnobr, &mnobr, &r__[r_offset], ldr, &iwork[1], 
			&toll, &svlmax, &dwork[itau2], &rank, sval, &dwork[
			jwork], &i__1, &ierr, (ftnlen)2);
/* Computing MAX */
#line 496 "IB01ND.f"
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 496 "IB01ND.f"
		maxwrk = max(i__1,i__2);

/*              Workspace: need   2*M*NOBR + (M+2*L)*NOBR; */
/*                         prefer 2*M*NOBR + (M+2*L)*NOBR*NB. */

#line 501 "IB01ND.f"
		i__1 = *ldwork - jwork + 1;
#line 501 "IB01ND.f"
		dormqr_("Left", "Transpose", &mnobr, &llmnob, &mnobr, &r__[
			r_offset], ldr, &dwork[itau2], &r__[nr2 * r_dim1 + 1],
			 ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)
			9);
/* Computing MAX */
#line 504 "IB01ND.f"
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 504 "IB01ND.f"
		maxwrk = max(i__1,i__2);
#line 505 "IB01ND.f"
		if (rank < mnobr) {

/*                 The least squares problem is rank-deficient. */

#line 509 "IB01ND.f"
		    *iwarn = 4;
#line 510 "IB01ND.f"
		}

/*              Determine residuals r_1 and r_2: premultiply by  Q  and */
/*              then by  q. */
/*              Workspace: need   2*M*NOBR + (M+2*L)*NOBR); */
/*                         prefer 2*M*NOBR + (M+2*L)*NOBR*NB. */

#line 517 "IB01ND.f"
		dlaset_("Full", &rank, &llmnob, &c_b50, &c_b50, &r__[nr2 * 
			r_dim1 + 1], ldr, (ftnlen)4);
#line 519 "IB01ND.f"
		i__1 = *ldwork - jwork + 1;
#line 519 "IB01ND.f"
		dormqr_("Left", "NoTranspose", &mnobr, &llmnob, &mnobr, &r__[
			r_offset], ldr, &dwork[itau2], &r__[nr2 * r_dim1 + 1],
			 ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)
			11);
/* Computing MAX */
#line 522 "IB01ND.f"
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 522 "IB01ND.f"
		maxwrk = max(i__1,i__2);
#line 523 "IB01ND.f"
		jwork = itau2;

/*              Restore the transformation  q. */

#line 527 "IB01ND.f"
		i__1 = mnobr - 1;
#line 527 "IB01ND.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 528 "IB01ND.f"
		    i__2 = mmnobr;
#line 528 "IB01ND.f"
		    for (j = nr2; j <= i__2; ++j) {
#line 529 "IB01ND.f"
			r__[j - mnobr + i__ + i__ * r_dim1] = r__[j + i__ * 
				r_dim1];
#line 530 "IB01ND.f"
/* L25: */
#line 530 "IB01ND.f"
		    }
#line 531 "IB01ND.f"
/* L30: */
#line 531 "IB01ND.f"
		}

#line 533 "IB01ND.f"
	    }

/*           Premultiply by the transformation  q  (apply transformations */
/*           in backward order). */
/*           Workspace: need   M*NOBR + (M+2*L)*NOBR; */
/*                      prefer larger. */

#line 540 "IB01ND.f"
	    i__1 = mnobr - 1;
#line 540 "IB01ND.f"
	    i__2 = *ldwork - jwork + 1;
#line 540 "IB01ND.f"
	    mb04iy_("Left", "NoTranspose", &mmnobr, &llmnob, &mnobr, &i__1, &
		    r__[r_offset], ldr, &dwork[itau], &r__[nr2 * r_dim1 + 1], 
		    ldr, &dwork[jwork], &i__2, &ierr, (ftnlen)4, (ftnlen)11);
/* Computing MAX */
#line 543 "IB01ND.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 543 "IB01ND.f"
	    maxwrk = max(i__1,i__2);

#line 545 "IB01ND.f"
	} else {

/*           Save  Y_f  (transposed) in the last block-row of  R. */

#line 549 "IB01ND.f"
	    ma02ad_("Full", &lmmnob, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &
		    r__[nr4 + r_dim1], ldr, (ftnlen)4);
#line 551 "IB01ND.f"
	    rcond1 = 1.;
#line 552 "IB01ND.f"
	}

/*        Triangularize the matrix  r_1  for determining the oblique */
/*        projection  P  in least squares problem in (1).  Exploit the */
/*        fact that the third block-row of r_1  has the structure */
/*        [ 0  T ],  where  T  is an upper triangular matrix.  Then apply */
/*        the corresponding transformations  Q'  to the matrix  r_2. */
/*        Workspace: need   2*M*NOBR; */
/*                   prefer   M*NOBR + M*NOBR*NB. */

#line 562 "IB01ND.f"
	i__1 = *ldwork - jwork + 1;
#line 562 "IB01ND.f"
	dgeqrf_(&mmnobr, &mnobr, &r__[nr2 * r_dim1 + 1], ldr, &dwork[itau], &
		dwork[jwork], &i__1, &ierr);

/*        Workspace: need   M*NOBR + 2*L*NOBR; */
/*                   prefer M*NOBR + 2*L*NOBR*NB. */

#line 568 "IB01ND.f"
	i__1 = *ldwork - jwork + 1;
#line 568 "IB01ND.f"
	dormqr_("Left", "Transpose", &mmnobr, &llnobr, &mnobr, &r__[nr2 * 
		r_dim1 + 1], ldr, &dwork[itau], &r__[nr3 * r_dim1 + 1], ldr, &
		dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)9);
#line 571 "IB01ND.f"
	nrsave = nr2;

#line 573 "IB01ND.f"
	itau2 = jwork;
#line 574 "IB01ND.f"
	jwork = itau2 + lnobr;
#line 575 "IB01ND.f"
	i__1 = lnobr - 1;
#line 575 "IB01ND.f"
	i__2 = *ldwork - jwork + 1;
#line 575 "IB01ND.f"
	mb04id_(&lmnobr, &lnobr, &i__1, &lnobr, &r__[nr2 + nr3 * r_dim1], ldr,
		 &r__[nr2 + nr4 * r_dim1], ldr, &dwork[itau2], &dwork[jwork], 
		&i__2, &ierr);
/* Computing MAX */
#line 578 "IB01ND.f"
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 578 "IB01ND.f"
	maxwrk = max(i__1,i__2);

/*        Check the condition of the triangular matrix of order  (m+l)*s */
/*        just determined, and decide to use pivoting or not. */
/*        Workspace: need 4*(M+L)*NOBR. */

#line 584 "IB01ND.f"
	dtrcon_("1-norm", "Upper", "NonUnit", &lmnobr, &r__[nr2 * r_dim1 + 1],
		 ldr, &rcond2, &dwork[jwork], &iwork[1], &ierr, (ftnlen)6, (
		ftnlen)5, (ftnlen)7);

#line 587 "IB01ND.f"
	if (*tol <= 0.) {
#line 587 "IB01ND.f"
	    toll = lmnobr * lmnobr * eps;
#line 587 "IB01ND.f"
	}
#line 589 "IB01ND.f"
	if (rcond2 <= max(toll,thresh)) {
#line 590 "IB01ND.f"
	    if (*m > 0) {

/*              Save information about  Q  in  R_11  (in the strict lower */
/*              triangle),  R_21  and  R_31  (transposed information). */

#line 595 "IB01ND.f"
		i__1 = mmnobr - 1;
#line 595 "IB01ND.f"
		dlacpy_("Lower", &i__1, &mnobr, &r__[nr2 * r_dim1 + 2], ldr, &
			r__[r_dim1 + 2], ldr, (ftnlen)5);
#line 597 "IB01ND.f"
		nrsave = 1;

#line 599 "IB01ND.f"
		i__1 = lmnobr;
#line 599 "IB01ND.f"
		for (i__ = nr2; i__ <= i__1; ++i__) {
#line 600 "IB01ND.f"
		    dcopy_(&mnobr, &r__[i__ + 1 + (mnobr + i__) * r_dim1], &
			    c__1, &r__[mnobr + i__ + r_dim1], ldr);
#line 602 "IB01ND.f"
/* L40: */
#line 602 "IB01ND.f"
		}

#line 604 "IB01ND.f"
	    }

#line 606 "IB01ND.f"
	    i__1 = lmnobr - 1;
#line 606 "IB01ND.f"
	    i__2 = lmnobr - 1;
#line 606 "IB01ND.f"
	    dlaset_("Lower", &i__1, &i__2, &c_b50, &c_b50, &r__[nr2 * r_dim1 
		    + 2], ldr, (ftnlen)5);

/*           Use QR factorization with column pivoting. */
/*           Workspace: need   5*(M+L)*NOBR+1. */
/*                      prefer 4*(M+L)*NOBR + ((M+L)*NOBR+1)*NB. */

#line 613 "IB01ND.f"
	    i__1 = lmnobr;
#line 613 "IB01ND.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 614 "IB01ND.f"
		iwork[i__] = 0;
#line 615 "IB01ND.f"
/* L50: */
#line 615 "IB01ND.f"
	    }

#line 617 "IB01ND.f"
	    itau3 = jwork;
#line 618 "IB01ND.f"
	    jwork = itau3 + lmnobr;
#line 619 "IB01ND.f"
	    svlmax = 0.;
#line 620 "IB01ND.f"
	    i__1 = *ldwork - jwork + 1;
#line 620 "IB01ND.f"
	    mb03od_("QR", &lmnobr, &lmnobr, &r__[nr2 * r_dim1 + 1], ldr, &
		    iwork[1], &toll, &svlmax, &dwork[itau3], &rank1, sval, &
		    dwork[jwork], &i__1, &ierr, (ftnlen)2);
/* Computing MAX */
#line 623 "IB01ND.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 623 "IB01ND.f"
	    maxwrk = max(i__1,i__2);

/*           Workspace: need   2*(M+L)*NOBR + L*NOBR; */
/*                      prefer 2*(M+L)*NOBR + L*NOBR*NB. */

#line 628 "IB01ND.f"
	    i__1 = *ldwork - jwork + 1;
#line 628 "IB01ND.f"
	    dormqr_("Left", "Transpose", &lmnobr, &lnobr, &lmnobr, &r__[nr2 * 
		    r_dim1 + 1], ldr, &dwork[itau3], &r__[nr4 * r_dim1 + 1], 
		    ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 631 "IB01ND.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 631 "IB01ND.f"
	    maxwrk = max(i__1,i__2);
#line 632 "IB01ND.f"
	    if (rank1 < lmnobr) {

/*              The least squares problem is rank-deficient. */

#line 636 "IB01ND.f"
		*iwarn = 5;
#line 637 "IB01ND.f"
	    }

/*           Apply the orthogonal transformations, in backward order, to */
/*           [r_2(1:rank(r_1),:)' 0]',  to obtain  P'. */
/*           Workspace: need   2*(M+L)*NOBR + L*NOBR; */
/*                      prefer 2*(M+L)*NOBR + L*NOBR*NB. */

#line 644 "IB01ND.f"
	    i__1 = lmnobr - rank1;
#line 644 "IB01ND.f"
	    dlaset_("Full", &i__1, &lnobr, &c_b50, &c_b50, &r__[rank1 + 1 + 
		    nr4 * r_dim1], ldr, (ftnlen)4);
#line 646 "IB01ND.f"
	    i__1 = *ldwork - jwork + 1;
#line 646 "IB01ND.f"
	    dormqr_("Left", "NoTranspose", &lmnobr, &lnobr, &lmnobr, &r__[nr2 
		    * r_dim1 + 1], ldr, &dwork[itau3], &r__[nr4 * r_dim1 + 1],
		     ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)11);
/* Computing MAX */
#line 649 "IB01ND.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 649 "IB01ND.f"
	    maxwrk = max(i__1,i__2);
#line 650 "IB01ND.f"
	    jwork = itau3;

#line 652 "IB01ND.f"
	    if (*m > 0) {

/*              Restore the saved transpose matrix from  R_31. */

#line 656 "IB01ND.f"
		i__1 = lmnobr;
#line 656 "IB01ND.f"
		for (i__ = nr2; i__ <= i__1; ++i__) {
#line 657 "IB01ND.f"
		    dcopy_(&mnobr, &r__[mnobr + i__ + r_dim1], ldr, &r__[i__ 
			    + 1 + (mnobr + i__) * r_dim1], &c__1);
#line 659 "IB01ND.f"
/* L60: */
#line 659 "IB01ND.f"
		}

#line 661 "IB01ND.f"
	    }

#line 663 "IB01ND.f"
	}

/*        Workspace: need   M*NOBR + L*NOBR; */
/*                   prefer larger. */

#line 668 "IB01ND.f"
	i__1 = lnobr - 1;
#line 668 "IB01ND.f"
	i__2 = *ldwork - jwork + 1;
#line 668 "IB01ND.f"
	mb04iy_("Left", "NoTranspose", &lmnobr, &lnobr, &lnobr, &i__1, &r__[
		nr2 + nr3 * r_dim1], ldr, &dwork[itau2], &r__[nr2 + nr4 * 
		r_dim1], ldr, &dwork[jwork], &i__2, &ierr, (ftnlen)4, (ftnlen)
		11);
/* Computing MAX */
#line 672 "IB01ND.f"
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 672 "IB01ND.f"
	maxwrk = max(i__1,i__2);

/*        Workspace: need   M*NOBR + L*NOBR; */
/*                   prefer M*NOBR + L*NOBR*NB. */

#line 677 "IB01ND.f"
	jwork = itau2;
#line 678 "IB01ND.f"
	i__1 = *ldwork - jwork + 1;
#line 678 "IB01ND.f"
	dormqr_("Left", "NoTranspose", &mmnobr, &lnobr, &mnobr, &r__[nrsave * 
		r_dim1 + 1], ldr, &dwork[itau], &r__[nr4 * r_dim1 + 1], ldr, &
		dwork[jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)11);

/*        Now, the matrix  P'  is available in  R_14 : R_34. */
/*        Triangularize the matrix  P'. */
/*        Workspace: need   2*L*NOBR; */
/*                   prefer   L*NOBR + L*NOBR*NB. */

#line 687 "IB01ND.f"
	jwork = itau + lnobr;
#line 688 "IB01ND.f"
	i__1 = *ldwork - jwork + 1;
#line 688 "IB01ND.f"
	dgeqrf_(&lmmnob, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &dwork[itau], &
		dwork[jwork], &i__1, &ierr);

/*        Copy the triangular factor to its final position,  R_22. */

#line 693 "IB01ND.f"
	dlacpy_("Upper", &lnobr, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &r__[
		nr2 + nr2 * r_dim1], ldr, (ftnlen)5);

/*        Restore  Y_f. */

#line 698 "IB01ND.f"
	ma02ad_("Full", &lnobr, &lmmnob, &r__[nr4 + r_dim1], ldr, &r__[nr4 * 
		r_dim1 + 1], ldr, (ftnlen)4);
#line 700 "IB01ND.f"
    }

/*     Find the singular value decomposition of  R_22. */
/*     Workspace: need 5*L*NOBR. */

#line 705 "IB01ND.f"
    mb03ud_("NoVectors", "Vectors", &lnobr, &r__[nr2 + nr2 * r_dim1], ldr, 
	    dum, &c__1, &sv[1], &dwork[1], ldwork, &ierr, (ftnlen)9, (ftnlen)
	    7);
#line 707 "IB01ND.f"
    if (ierr != 0) {
#line 708 "IB01ND.f"
	*info = 2;
#line 709 "IB01ND.f"
	return 0;
#line 710 "IB01ND.f"
    }
/* Computing MAX */
#line 711 "IB01ND.f"
    i__1 = maxwrk, i__2 = (integer) dwork[1];
#line 711 "IB01ND.f"
    maxwrk = max(i__1,i__2);

/*     Transpose  R(m*s+1:(m+L)*s,m*s+1:(m+L)*s)  in-situ; its */
/*     columns will then be the singular vectors needed subsequently. */

#line 716 "IB01ND.f"
    i__1 = lmnobr;
#line 716 "IB01ND.f"
    for (i__ = nr2 + 1; i__ <= i__1; ++i__) {
#line 717 "IB01ND.f"
	i__2 = lmnobr - i__ + 1;
#line 717 "IB01ND.f"
	dswap_(&i__2, &r__[i__ + (i__ - 1) * r_dim1], &c__1, &r__[i__ - 1 + 
		i__ * r_dim1], ldr);
#line 718 "IB01ND.f"
/* L70: */
#line 718 "IB01ND.f"
    }

/*     Return optimal workspace in  DWORK(1)  and reciprocal condition */
/*     numbers, if  METH = 'N'. */

#line 723 "IB01ND.f"
    dwork[1] = (doublereal) maxwrk;
#line 724 "IB01ND.f"
    if (n4sid) {
#line 725 "IB01ND.f"
	dwork[2] = rcond1;
#line 726 "IB01ND.f"
	dwork[3] = rcond2;
#line 727 "IB01ND.f"
    }
#line 728 "IB01ND.f"
    return 0;

/* *** Last line of IB01ND *** */
} /* ib01nd_ */

