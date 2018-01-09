#line 1 "IB01MY.f"
/* IB01MY.f -- translated by f2c (version 20100827).
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

#line 1 "IB01MY.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b18 = 1.;
static integer c__0 = 0;
static doublereal c_b111 = 0.;

/* Subroutine */ int ib01my_(char *meth, char *batch, char *conct, integer *
	nobr, integer *m, integer *l, integer *nsmp, doublereal *u, integer *
	ldu, doublereal *y, integer *ldy, doublereal *r__, integer *ldr, 
	integer *iwork, doublereal *dwork, integer *ldwork, integer *iwarn, 
	integer *info, ftnlen meth_len, ftnlen batch_len, ftnlen conct_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j, k, jd;
    static doublereal cs;
    static integer nr;
    static doublereal sn;
    static integer ns, icj, ing, ipg, jds;
    static doublereal dum[1];
    static integer nrg;
    static doublereal upd, tau;
    static integer nsm, ipy;
    static doublereal beta;
    static integer ingc, ipgc, icol, imax, ingp, ierr, itau, lnrg, mnrg, irev;
    static logical last, n4sid;
    static integer nobr2;
    extern /* Subroutine */ int ma02ed_(char *, integer *, doublereal *, 
	    integer *, ftnlen), ma02fd_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *), mb04id_(integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *), mb04od_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *),
	     dlarf_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, ftnlen), 
	    dgemm_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer nobr21, iconn, lnobr, mnobr;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical moesp, first;
    static integer jwork;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    static logical onebch;
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern integer idamax_(integer *, doublereal *, integer *);
    static logical connec;
    static integer icycle;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen), xerbla_(char *, integer *, ftnlen);
    static integer llnobr, mmnobr;
    static logical interm;
    extern /* Subroutine */ int dormqr_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen);
    static integer ldrwrk, minwrk, maxwrk, nsmpsm;


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

/*     To construct an upper triangular factor  R  of the concatenated */
/*     block Hankel matrices using input-output data, via a fast QR */
/*     algorithm based on displacement rank.  The input-output data can, */
/*     optionally, be processed sequentially. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     METH    CHARACTER*1 */
/*             Specifies the subspace identification method to be used, */
/*             as follows: */
/*             = 'M':  MOESP  algorithm with past inputs and outputs; */
/*             = 'N':  N4SID  algorithm. */

/*     BATCH   CHARACTER*1 */
/*             Specifies whether or not sequential data processing is to */
/*             be used, and, for sequential processing, whether or not */
/*             the current data block is the first block, an intermediate */
/*             block, or the last block, as follows: */
/*             = 'F':  the first block in sequential data processing; */
/*             = 'I':  an intermediate block in sequential data */
/*                     processing; */
/*             = 'L':  the last block in sequential data processing; */
/*             = 'O':  one block only (non-sequential data processing). */
/*             NOTE that when  100  cycles of sequential data processing */
/*                  are completed for  BATCH = 'I',  a warning is */
/*                  issued, to prevent for an infinite loop. */

/*     CONCT   CHARACTER*1 */
/*             Specifies whether or not the successive data blocks in */
/*             sequential data processing belong to a single experiment, */
/*             as follows: */
/*             = 'C':  the current data block is a continuation of the */
/*                     previous data block and/or it will be continued */
/*                     by the next data block; */
/*             = 'N':  there is no connection between the current data */
/*                     block and the previous and/or the next ones. */
/*             This parameter is not used if BATCH = 'O'. */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             The number of block rows,  s,  in the input and output */
/*             block Hankel matrices to be processed.  NOBR > 0. */
/*             (In the MOESP theory,  NOBR  should be larger than  n, the */
/*             estimated dimension of state vector.) */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */
/*             When M = 0, no system inputs are processed. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     NSMP    (input) INTEGER */
/*             The number of rows of matrices  U  and  Y  (number of */
/*             samples,  t). (When sequential data processing is used, */
/*             NSMP  is the number of samples of the current data */
/*             block.) */
/*             NSMP >= 2*(M+L+1)*NOBR - 1,  for non-sequential */
/*                                          processing; */
/*             NSMP >= 2*NOBR,  for sequential processing. */
/*             The total number of samples when calling the routine with */
/*             BATCH = 'L'  should be at least  2*(M+L+1)*NOBR - 1. */
/*             The  NSMP  argument may vary from a cycle to another in */
/*             sequential data processing, but  NOBR, M,  and  L  should */
/*             be kept constant. For efficiency, it is advisable to use */
/*             NSMP  as large as possible. */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,M) */
/*             The leading NSMP-by-M part of this array must contain the */
/*             t-by-m input-data sequence matrix  U, */
/*             U = [u_1 u_2 ... u_m].  Column  j  of  U  contains the */
/*             NSMP  values of the j-th input component for consecutive */
/*             time increments. */
/*             If M = 0, this array is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= NSMP, if M > 0; */
/*             LDU >= 1,    if M = 0. */

/*     Y       (input) DOUBLE PRECISION array, dimension (LDY,L) */
/*             The leading NSMP-by-L part of this array must contain the */
/*             t-by-l output-data sequence matrix  Y, */
/*             Y = [y_1 y_2 ... y_l].  Column  j  of  Y  contains the */
/*             NSMP  values of the j-th output component for consecutive */
/*             time increments. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= NSMP. */

/*     R       (output) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             If INFO = 0 and BATCH = 'L' or 'O', the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of this */
/*             array contains the upper triangular factor R from the */
/*             QR factorization of the concatenated block Hankel */
/*             matrices. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= 2*(M+L)*NOBR. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (M+L) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -16,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */
/*             The first (M+L)*2*NOBR*(M+L+c) elements of  DWORK  should */
/*             be preserved during successive calls of the routine */
/*             with  BATCH = 'F'  or  'I',  till the final call with */
/*             BATCH = 'L',  where */
/*             c = 1,  if the successive data blocks do not belong to a */
/*                     single experiment  (CONCT = 'N'); */
/*             c = 2,  if the successive data blocks belong to a single */
/*                     experiment  (CONCT = 'C'). */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+3), */
/*                              if BATCH <> 'O' and CONCT = 'C'; */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+1), */
/*                              if BATCH = 'F' or 'I' and CONCT = 'N'; */
/*             LDWORK >= (M+L)*4*NOBR*(M+L+1)+(M+L)*2*NOBR, */
/*                              if BATCH = 'L' and CONCT = 'N', */
/*                              or BATCH = 'O'. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  the number of 100 cycles in sequential data */
/*                   processing has been exhausted without signaling */
/*                   that the last block of data was get. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the fast QR factorization algorithm failed. The */
/*                   matrix H'*H is not (numerically) positive definite. */

/*     METHOD */

/*     Consider the  t x 2(m+l)s  matrix H of concatenated block Hankel */
/*     matrices */

/*          H = [ Uf'         Up'      Y'      ],  for METH = 'M', */
/*                  s+1,2s,t    1,s,t   1,2s,t */

/*          H = [ U'       Y'      ],              for METH = 'N', */
/*                 1,2s,t   1,2s,t */

/*     where  Up     , Uf        , U      , and  Y        are block */
/*              1,s,t    s+1,2s,t   1,2s,t        1,2s,t */
/*     Hankel matrices defined in terms of the input and output data [3]. */
/*     The fast QR algorithm uses a factorization of H'*H which exploits */
/*     the block-Hankel structure, via a displacement rank technique [5]. */

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

/*     [5] Kressner, D., Mastronardi, N., Sima, V., Van Dooren, P., and */
/*         Van Huffel, S. */
/*         A Fast Algorithm for Subspace State-space System */
/*         Identification via Exploitation of the Displacement Structure. */
/*         J. Comput. Appl. Math., Vol.132, No.1, pp. 71-81, 2001. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is reliable and efficient. Numerical */
/*     difficulties are possible when the matrix H'*H is nearly rank */
/*     defficient. The method cannot be used if the matrix H'*H is not */
/*     numerically positive definite. */
/*                                     2           3 2 */
/*     The algorithm requires 0(2t(m+l) s)+0(4(m+l) s ) floating point */
/*     operations. */

/*     CONTRIBUTORS */

/*     V. Sima, Katholieke Universiteit Leuven, June 2000. */
/*     Partly based on Matlab codes developed by N. Mastronardi, */
/*     Katholieke Universiteit Leuven, February 2000. */

/*     REVISIONS */

/*     V. Sima, July 2000, August 2000, Feb. 2004, May 2009. */

/*     KEYWORDS */

/*     Displacement rank, Hankel matrix, Householder transformation, */
/*     identification methods, multivariable systems. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Save Statement .. */
/*        ICYCLE  is used to count the cycles for  BATCH = 'I'. */
/*        MAXWRK  is used to store the optimal workspace. */
/*        NSMPSM  is used to sum up the  NSMP  values for  BATCH <> 'O'. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 291 "IB01MY.f"
    /* Parameter adjustments */
#line 291 "IB01MY.f"
    u_dim1 = *ldu;
#line 291 "IB01MY.f"
    u_offset = 1 + u_dim1;
#line 291 "IB01MY.f"
    u -= u_offset;
#line 291 "IB01MY.f"
    y_dim1 = *ldy;
#line 291 "IB01MY.f"
    y_offset = 1 + y_dim1;
#line 291 "IB01MY.f"
    y -= y_offset;
#line 291 "IB01MY.f"
    r_dim1 = *ldr;
#line 291 "IB01MY.f"
    r_offset = 1 + r_dim1;
#line 291 "IB01MY.f"
    r__ -= r_offset;
#line 291 "IB01MY.f"
    --iwork;
#line 291 "IB01MY.f"
    --dwork;
#line 291 "IB01MY.f"

#line 291 "IB01MY.f"
    /* Function Body */
#line 291 "IB01MY.f"
    moesp = lsame_(meth, "M", (ftnlen)1, (ftnlen)1);
#line 292 "IB01MY.f"
    n4sid = lsame_(meth, "N", (ftnlen)1, (ftnlen)1);
#line 293 "IB01MY.f"
    onebch = lsame_(batch, "O", (ftnlen)1, (ftnlen)1);
#line 294 "IB01MY.f"
    first = lsame_(batch, "F", (ftnlen)1, (ftnlen)1) || onebch;
#line 295 "IB01MY.f"
    interm = lsame_(batch, "I", (ftnlen)1, (ftnlen)1);
#line 296 "IB01MY.f"
    last = lsame_(batch, "L", (ftnlen)1, (ftnlen)1) || onebch;
#line 297 "IB01MY.f"
    if (! onebch) {
#line 298 "IB01MY.f"
	connec = lsame_(conct, "C", (ftnlen)1, (ftnlen)1);
#line 299 "IB01MY.f"
    } else {
#line 300 "IB01MY.f"
	connec = FALSE_;
#line 301 "IB01MY.f"
    }
#line 302 "IB01MY.f"
    mnobr = *m * *nobr;
#line 303 "IB01MY.f"
    lnobr = *l * *nobr;
#line 304 "IB01MY.f"
    mmnobr = mnobr + mnobr;
#line 305 "IB01MY.f"
    llnobr = lnobr + lnobr;
#line 306 "IB01MY.f"
    nobr2 = *nobr << 1;
#line 307 "IB01MY.f"
    nobr21 = nobr2 - 1;
#line 308 "IB01MY.f"
    *iwarn = 0;
#line 309 "IB01MY.f"
    *info = 0;
#line 310 "IB01MY.f"
    if (first) {
#line 311 "IB01MY.f"
	icycle = 1;
#line 312 "IB01MY.f"
	maxwrk = 1;
#line 313 "IB01MY.f"
	nsmpsm = 0;
#line 314 "IB01MY.f"
    }
#line 315 "IB01MY.f"
    nsmpsm += *nsmp;
#line 316 "IB01MY.f"
    nr = mmnobr + llnobr;

/*     Check the scalar input parameters. */

#line 320 "IB01MY.f"
    if (! (moesp || n4sid)) {
#line 321 "IB01MY.f"
	*info = -1;
#line 322 "IB01MY.f"
    } else if (! (first || interm || last)) {
#line 323 "IB01MY.f"
	*info = -2;
#line 324 "IB01MY.f"
    } else if (! onebch) {
#line 325 "IB01MY.f"
	if (! (connec || lsame_(conct, "N", (ftnlen)1, (ftnlen)1))) {
#line 325 "IB01MY.f"
	    *info = -3;
#line 325 "IB01MY.f"
	}
#line 327 "IB01MY.f"
    }
#line 328 "IB01MY.f"
    if (*info == 0) {
#line 329 "IB01MY.f"
	if (*nobr <= 0) {
#line 330 "IB01MY.f"
	    *info = -4;
#line 331 "IB01MY.f"
	} else if (*m < 0) {
#line 332 "IB01MY.f"
	    *info = -5;
#line 333 "IB01MY.f"
	} else if (*l <= 0) {
#line 334 "IB01MY.f"
	    *info = -6;
#line 335 "IB01MY.f"
	} else if (*nsmp < nobr2 || last && nsmpsm < nr + nobr21) {
#line 337 "IB01MY.f"
	    *info = -7;
#line 338 "IB01MY.f"
	} else if (*ldu < 1 || *m > 0 && *ldu < *nsmp) {
#line 339 "IB01MY.f"
	    *info = -9;
#line 340 "IB01MY.f"
	} else if (*ldy < *nsmp) {
#line 341 "IB01MY.f"
	    *info = -11;
#line 342 "IB01MY.f"
	} else if (*ldr < nr) {
#line 343 "IB01MY.f"
	    *info = -13;
#line 344 "IB01MY.f"
	} else {

/*           Compute workspace. */
/*           NRG is the number of positive (or negative) generators. */

#line 349 "IB01MY.f"
	    nrg = *m + *l + 1;
#line 350 "IB01MY.f"
	    if (! onebch && connec) {
#line 351 "IB01MY.f"
		minwrk = nr * (nrg + 2);
#line 352 "IB01MY.f"
	    } else if (first || interm) {
#line 353 "IB01MY.f"
		minwrk = nr * nrg;
#line 354 "IB01MY.f"
	    } else {
#line 355 "IB01MY.f"
		minwrk = (nr << 1) * nrg + nr;
#line 356 "IB01MY.f"
	    }
#line 357 "IB01MY.f"
	    maxwrk = max(minwrk,maxwrk);

#line 359 "IB01MY.f"
	    if (*ldwork < minwrk) {
#line 359 "IB01MY.f"
		*info = -16;
#line 359 "IB01MY.f"
	    }
#line 361 "IB01MY.f"
	}
#line 362 "IB01MY.f"
    }

/*     Return if there are illegal arguments. */

#line 366 "IB01MY.f"
    if (*info != 0) {
#line 367 "IB01MY.f"
	nsmpsm = 0;
#line 368 "IB01MY.f"
	if (*info == -16) {
#line 368 "IB01MY.f"
	    dwork[1] = (doublereal) minwrk;
#line 368 "IB01MY.f"
	}
#line 370 "IB01MY.f"
	i__1 = -(*info);
#line 370 "IB01MY.f"
	xerbla_("IB01MY", &i__1, (ftnlen)6);
#line 371 "IB01MY.f"
	return 0;
#line 372 "IB01MY.f"
    }

/*     Compute the  R  factor from a fast QR factorization of the */
/*     matrix  H,  a concatenation of two block Hankel matrices. */
/*     Specifically, a displacement rank technique is applied to */
/*     the block Toeplitz matrix,  G = (P*H)'*(P*H),  where  P  is a */
/*     2-by-2 block diagonal matrix, having as diagonal blocks identity */
/*     matrices with columns taken in the reverse order. */
/*     The technique builds and processes the generators of  G.  The */
/*     matrices  G  and  G1 = H'*H  have the same  R  factor. */

/*     Set the parameters for constructing the correlations of the */
/*     current block. */
/*     NSM is the number of processed samples in U and Y, t - 2s. */
/*     IPG and ING are pointers to the "positive" and "negative" */
/*     generators, stored row-wise in the workspace. All "positive" */
/*     generators are stored before any "negative" generators. */
/*     If BATCH <> 'O' and CONCT = 'C', the "connection" elements of */
/*     two successive batches are stored in the same workspace as the */
/*     "negative" generators (which will be computed later on). */
/*     IPY is a pointer to the Y part of the "positive" generators. */
/*     LDRWRK is used as a leading dimension for the workspace part used */
/*     to store the "connection" elements. */

#line 396 "IB01MY.f"
    ns = *nsmp - nobr21;
#line 397 "IB01MY.f"
    nsm = ns - 1;
#line 398 "IB01MY.f"
    mnrg = *m * nrg;
#line 399 "IB01MY.f"
    lnrg = *l * nrg;

#line 401 "IB01MY.f"
    ldrwrk = nobr2 << 1;
#line 402 "IB01MY.f"
    if (first) {
#line 403 "IB01MY.f"
	upd = 0.;
#line 404 "IB01MY.f"
    } else {
#line 405 "IB01MY.f"
	upd = 1.;
#line 406 "IB01MY.f"
    }
#line 407 "IB01MY.f"
    dum[0] = 0.;

#line 409 "IB01MY.f"
    ipg = 1;
#line 410 "IB01MY.f"
    ipy = ipg + *m;
#line 411 "IB01MY.f"
    ing = ipg + nrg * nr;
#line 412 "IB01MY.f"
    iconn = ing;

#line 414 "IB01MY.f"
    if (! first && connec) {

/*        Restore the saved (M+L)*2*NOBR "connection" elements of */
/*        U  and  Y  into their appropriate position in sequential */
/*        processing. The process is performed column-wise, in */
/*        reverse order, first for  Y  and then for  U. */
/*        ICONN is a pointer to the first saved "connection" element. */
/*        Workspace: need   (M+L)*2*NOBR*(M+L+3). */

#line 423 "IB01MY.f"
	irev = iconn + nr;
#line 424 "IB01MY.f"
	icol = iconn + (nr << 1);

#line 426 "IB01MY.f"
	i__1 = *m + *l;
#line 426 "IB01MY.f"
	for (i__ = 2; i__ <= i__1; ++i__) {
#line 427 "IB01MY.f"
	    irev -= nobr2;
#line 428 "IB01MY.f"
	    icol -= ldrwrk;
#line 429 "IB01MY.f"
	    dcopy_(&nobr2, &dwork[irev], &c__1, &dwork[icol], &c__1);
#line 430 "IB01MY.f"
/* L10: */
#line 430 "IB01MY.f"
	}

#line 432 "IB01MY.f"
	if (*m > 0) {
#line 432 "IB01MY.f"
	    dlacpy_("Full", &nobr2, m, &u[u_offset], ldu, &dwork[iconn + 
		    nobr2], &ldrwrk, (ftnlen)4);
#line 432 "IB01MY.f"
	}
#line 435 "IB01MY.f"
	dlacpy_("Full", &nobr2, l, &y[y_offset], ldy, &dwork[iconn + ldrwrk * 
		*m + nobr2], &ldrwrk, (ftnlen)4);
#line 437 "IB01MY.f"
    }

#line 439 "IB01MY.f"
    if (*m > 0) {

/*        Let  Guu(i,j) = Guu0(i,j) + u_i*u_j' + u_(i+1)*u_(j+1)' + */
/*                              ... + u_(i+NSM-1)*u_(j+NSM-1)', */
/*        where  u_i'  is the i-th row of  U,  j = 1 : 2s,  i = 1 : j, */
/*        NSM = NSMP - 2s,  and  Guu0(i,j)  is a zero matrix for */
/*        BATCH = 'O' or 'F', and it is the matrix Guu(i,j) computed */
/*        till the current block for BATCH = 'I' or 'L'. The matrix */
/*        Guu(i,j)  is  m-by-m,  and  Guu(j,j)  is symmetric. The */
/*        submatrices of the first block-row, Guu(1,j), are needed only. */

/*        Compute/update  Guu(1,1). */

#line 452 "IB01MY.f"
	if (! first && connec) {
#line 452 "IB01MY.f"
	    dsyrk_("Upper", "Transpose", m, &nobr2, &c_b18, &dwork[iconn], &
		    ldrwrk, &upd, &dwork[ipg], &nrg, (ftnlen)5, (ftnlen)9);
#line 452 "IB01MY.f"
	}
#line 455 "IB01MY.f"
	dsyrk_("Upper", "Transpose", m, &nsm, &c_b18, &u[u_offset], ldu, &upd,
		 &dwork[ipg], &nrg, (ftnlen)5, (ftnlen)9);
#line 457 "IB01MY.f"
	ma02ed_("Upper", m, &dwork[ipg], &nrg, (ftnlen)5);

#line 459 "IB01MY.f"
	jd = 1;

#line 461 "IB01MY.f"
	if (first || ! connec) {

#line 463 "IB01MY.f"
	    i__1 = nobr2;
#line 463 "IB01MY.f"
	    for (j = 2; j <= i__1; ++j) {
#line 464 "IB01MY.f"
		jd += *m;

/*              Compute/update  Guu(1,j). */

#line 468 "IB01MY.f"
		dgemm_("Transpose", "NoTranspose", m, m, &nsm, &c_b18, &u[
			u_offset], ldu, &u[j + u_dim1], ldu, &upd, &dwork[ipg 
			+ (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
#line 471 "IB01MY.f"
/* L20: */
#line 471 "IB01MY.f"
	    }

#line 473 "IB01MY.f"
	} else {

#line 475 "IB01MY.f"
	    i__1 = nobr2;
#line 475 "IB01MY.f"
	    for (j = 2; j <= i__1; ++j) {
#line 476 "IB01MY.f"
		jd += *m;

/*              Compute/update  Guu(1,j)  for sequential processing */
/*              with connected blocks. */

#line 481 "IB01MY.f"
		dgemm_("Transpose", "NoTranspose", m, m, &nobr2, &c_b18, &
			dwork[iconn], &ldrwrk, &dwork[iconn + j - 1], &ldrwrk,
			 &upd, &dwork[ipg + (jd - 1) * nrg], &nrg, (ftnlen)9, 
			(ftnlen)11);
#line 484 "IB01MY.f"
		dgemm_("Transpose", "NoTranspose", m, m, &nsm, &c_b18, &u[
			u_offset], ldu, &u[j + u_dim1], ldu, &c_b18, &dwork[
			ipg + (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
#line 487 "IB01MY.f"
/* L30: */
#line 487 "IB01MY.f"
	    }

#line 489 "IB01MY.f"
	}

/*        Let  Guy(i,j) = Guy0(i,j) + u_i*y_j' + u_(i+1)*y_(j+1)' + */
/*                              ... + u_(i+NSM-1)*y_(j+NSM-1)', */
/*        where  u_i'  is the i-th row of  U,  y_j'  is the j-th row */
/*        of  Y,  j = 1 : 2s,  i = 1 : 2s,  NSM = NSMP - 2s,  and */
/*        Guy0(i,j)  is a zero matrix for  BATCH = 'O' or 'F', and it */
/*        is the matrix Guy(i,j) computed till the current block for */
/*        BATCH = 'I' or 'L'.  Guy(i,j) is m-by-L. The submatrices */
/*        of the first block-row, Guy(1,j), as well as the transposes */
/*        of the submatrices of the first block-column, i.e., Gyu(1,j), */
/*        are needed only. */

#line 502 "IB01MY.f"
	jd = mmnobr + 1;

#line 504 "IB01MY.f"
	if (first || ! connec) {

#line 506 "IB01MY.f"
	    i__1 = nobr2;
#line 506 "IB01MY.f"
	    for (j = 1; j <= i__1; ++j) {

/*              Compute/update  Guy(1,j). */

#line 510 "IB01MY.f"
		dgemm_("Transpose", "NoTranspose", m, l, &nsm, &c_b18, &u[
			u_offset], ldu, &y[j + y_dim1], ldy, &upd, &dwork[ipg 
			+ (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
#line 513 "IB01MY.f"
		jd += *l;
#line 514 "IB01MY.f"
/* L40: */
#line 514 "IB01MY.f"
	    }

#line 516 "IB01MY.f"
	} else {

#line 518 "IB01MY.f"
	    i__1 = nobr2;
#line 518 "IB01MY.f"
	    for (j = 1; j <= i__1; ++j) {

/*              Compute/update  Guy(1,j)  for sequential processing */
/*              with connected blocks. */

#line 523 "IB01MY.f"
		dgemm_("Transpose", "NoTranspose", m, l, &nobr2, &c_b18, &
			dwork[iconn], &ldrwrk, &dwork[iconn + ldrwrk * *m + j 
			- 1], &ldrwrk, &upd, &dwork[ipg + (jd - 1) * nrg], &
			nrg, (ftnlen)9, (ftnlen)11);
#line 527 "IB01MY.f"
		dgemm_("Transpose", "NoTranspose", m, l, &nsm, &c_b18, &u[
			u_offset], ldu, &y[j + y_dim1], ldy, &c_b18, &dwork[
			ipg + (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
#line 530 "IB01MY.f"
		jd += *l;
#line 531 "IB01MY.f"
/* L50: */
#line 531 "IB01MY.f"
	    }

#line 533 "IB01MY.f"
	}

/*        Now, the first M "positive" generators have been built. */
/*        Transpose  Guy(1,1)  in the first block of the  Y  part of the */
/*        "positive" generators. */

#line 539 "IB01MY.f"
	i__1 = *l;
#line 539 "IB01MY.f"
	for (j = 1; j <= i__1; ++j) {
#line 540 "IB01MY.f"
	    dcopy_(m, &dwork[ipg + (mmnobr + j - 1) * nrg], &c__1, &dwork[ipy 
		    + j - 1], &nrg);
#line 542 "IB01MY.f"
/* L60: */
#line 542 "IB01MY.f"
	}

#line 544 "IB01MY.f"
	jd = 1;

#line 546 "IB01MY.f"
	if (first || ! connec) {

#line 548 "IB01MY.f"
	    i__1 = nobr2;
#line 548 "IB01MY.f"
	    for (j = 2; j <= i__1; ++j) {
#line 549 "IB01MY.f"
		jd += *m;

/*              Compute/update  Gyu(1,j). */

#line 553 "IB01MY.f"
		dgemm_("Transpose", "NoTranspose", l, m, &nsm, &c_b18, &y[
			y_offset], ldy, &u[j + u_dim1], ldu, &upd, &dwork[ipy 
			+ (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
#line 556 "IB01MY.f"
/* L70: */
#line 556 "IB01MY.f"
	    }

#line 558 "IB01MY.f"
	} else {

#line 560 "IB01MY.f"
	    i__1 = nobr2;
#line 560 "IB01MY.f"
	    for (j = 2; j <= i__1; ++j) {
#line 561 "IB01MY.f"
		jd += *m;

/*              Compute/update  Gyu(1,j)  for sequential processing */
/*              with connected blocks. */

#line 566 "IB01MY.f"
		dgemm_("Transpose", "NoTranspose", l, m, &nobr2, &c_b18, &
			dwork[iconn + ldrwrk * *m], &ldrwrk, &dwork[iconn + j 
			- 1], &ldrwrk, &upd, &dwork[ipy + (jd - 1) * nrg], &
			nrg, (ftnlen)9, (ftnlen)11);
#line 570 "IB01MY.f"
		dgemm_("Transpose", "NoTranspose", l, m, &nsm, &c_b18, &y[
			y_offset], ldy, &u[j + u_dim1], ldu, &c_b18, &dwork[
			ipy + (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
#line 573 "IB01MY.f"
/* L80: */
#line 573 "IB01MY.f"
	    }

#line 575 "IB01MY.f"
	}

#line 577 "IB01MY.f"
    }

/*     Let  Gyy(i,j) = Gyy0(i,j) + y_i*y_i' + y_(i+1)*y_(i+1)' + ... + */
/*                                 y_(i+NSM-1)*y_(i+NSM-1)', */
/*     where  y_i'  is the i-th row of  Y,  j = 1 : 2s,  i = 1 : j, */
/*     NSM = NSMP - 2s,  and  Gyy0(i,j)  is a zero matrix for */
/*     BATCH = 'O' or 'F', and it is the matrix Gyy(i,j) computed till */
/*     the current block for BATCH = 'I' or 'L'.  Gyy(i,j) is L-by-L, */
/*     and  Gyy(j,j)  is symmetric. The submatrices of the first */
/*     block-row, Gyy(1,j), are needed only. */

#line 588 "IB01MY.f"
    jd = mmnobr + 1;

/*     Compute/update  Gyy(1,1). */

#line 592 "IB01MY.f"
    if (! first && connec) {
#line 592 "IB01MY.f"
	dsyrk_("Upper", "Transpose", l, &nobr2, &c_b18, &dwork[iconn + ldrwrk 
		* *m], &ldrwrk, &upd, &dwork[ipy + mmnobr * nrg], &nrg, (
		ftnlen)5, (ftnlen)9);
#line 592 "IB01MY.f"
    }
#line 596 "IB01MY.f"
    dsyrk_("Upper", "Transpose", l, &nsm, &c_b18, &y[y_offset], ldy, &upd, &
	    dwork[ipy + mmnobr * nrg], &nrg, (ftnlen)5, (ftnlen)9);
#line 598 "IB01MY.f"
    ma02ed_("Upper", l, &dwork[ipy + mmnobr * nrg], &nrg, (ftnlen)5);

#line 600 "IB01MY.f"
    if (first || ! connec) {

#line 602 "IB01MY.f"
	i__1 = nobr2;
#line 602 "IB01MY.f"
	for (j = 2; j <= i__1; ++j) {
#line 603 "IB01MY.f"
	    jd += *l;

/*           Compute/update  Gyy(1,j). */

#line 607 "IB01MY.f"
	    dgemm_("Transpose", "NoTranspose", l, l, &nsm, &c_b18, &y[
		    y_offset], ldy, &y[j + y_dim1], ldy, &upd, &dwork[ipy + (
		    jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
#line 610 "IB01MY.f"
/* L90: */
#line 610 "IB01MY.f"
	}

#line 612 "IB01MY.f"
    } else {

#line 614 "IB01MY.f"
	i__1 = nobr2;
#line 614 "IB01MY.f"
	for (j = 2; j <= i__1; ++j) {
#line 615 "IB01MY.f"
	    jd += *l;

/*           Compute/update  Gyy(1,j)  for sequential processing with */
/*           connected blocks. */

#line 620 "IB01MY.f"
	    dgemm_("Transpose", "NoTranspose", l, l, &nobr2, &c_b18, &dwork[
		    iconn + ldrwrk * *m], &ldrwrk, &dwork[iconn + ldrwrk * *m 
		    + j - 1], &ldrwrk, &upd, &dwork[ipy + (jd - 1) * nrg], &
		    nrg, (ftnlen)9, (ftnlen)11);
#line 624 "IB01MY.f"
	    dgemm_("Transpose", "NoTranspose", l, l, &nsm, &c_b18, &y[
		    y_offset], ldy, &y[j + y_dim1], ldy, &c_b18, &dwork[ipy + 
		    (jd - 1) * nrg], &nrg, (ftnlen)9, (ftnlen)11);
#line 627 "IB01MY.f"
/* L100: */
#line 627 "IB01MY.f"
	}

#line 629 "IB01MY.f"
    }

#line 631 "IB01MY.f"
    if (! last) {
#line 632 "IB01MY.f"
	if (first) {

/*           For sequential processing, save the first 2*NOBR-1 rows of */
/*           the first block of  U  and  Y  in the appropriate */
/*           (M+L)*(2*NOBR-1) locations of  DWORK  starting at (1+M)*NRG. */
/*           These will be used to construct the last negative generator. */

#line 639 "IB01MY.f"
	    jd = nrg;
#line 640 "IB01MY.f"
	    if (*m > 0) {
#line 641 "IB01MY.f"
		dcopy_(m, dum, &c__0, &dwork[jd], &nrg);

#line 643 "IB01MY.f"
		i__1 = nobr21;
#line 643 "IB01MY.f"
		for (j = 1; j <= i__1; ++j) {
#line 644 "IB01MY.f"
		    jd += mnrg;
#line 645 "IB01MY.f"
		    dcopy_(m, &u[j + u_dim1], ldu, &dwork[jd], &nrg);
#line 646 "IB01MY.f"
/* L110: */
#line 646 "IB01MY.f"
		}

#line 648 "IB01MY.f"
		jd += mnrg;
#line 649 "IB01MY.f"
	    }
#line 650 "IB01MY.f"
	    dcopy_(l, dum, &c__0, &dwork[jd], &nrg);

#line 652 "IB01MY.f"
	    i__1 = nobr21;
#line 652 "IB01MY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 653 "IB01MY.f"
		jd += lnrg;
#line 654 "IB01MY.f"
		dcopy_(l, &y[j + y_dim1], ldy, &dwork[jd], &nrg);
#line 655 "IB01MY.f"
/* L120: */
#line 655 "IB01MY.f"
	    }

#line 657 "IB01MY.f"
	}

#line 659 "IB01MY.f"
	if (connec) {

/*           For sequential processing with connected data blocks, */
/*           save the remaining ("connection") elements of  U  and  Y */
/*           in (M+L)*2*NOBR locations of  DWORK  starting at ICONN. */

#line 665 "IB01MY.f"
	    if (*m > 0) {
#line 665 "IB01MY.f"
		dlacpy_("Full", &nobr2, m, &u[ns + u_dim1], ldu, &dwork[iconn]
			, &nobr2, (ftnlen)4);
#line 665 "IB01MY.f"
	    }
#line 668 "IB01MY.f"
	    dlacpy_("Full", &nobr2, l, &y[ns + y_dim1], ldy, &dwork[iconn + 
		    mmnobr], &nobr2, (ftnlen)4);
#line 670 "IB01MY.f"
	}

/*        Return to get new data. */

#line 674 "IB01MY.f"
	++icycle;
#line 675 "IB01MY.f"
	if (icycle > 100) {
#line 675 "IB01MY.f"
	    *iwarn = 1;
#line 675 "IB01MY.f"
	}
#line 677 "IB01MY.f"
	return 0;
#line 678 "IB01MY.f"
    }

#line 680 "IB01MY.f"
    if (last) {

/*        Try to compute the R factor. */

/*        Scale the first M+L positive generators and set the first */
/*        M+L negative generators. */
/*        Workspace: need   (M+L)*4*NOBR*(M+L+1)+M+L. */

#line 688 "IB01MY.f"
	jwork = (nrg << 1) * nr + 1;
#line 689 "IB01MY.f"
	i__1 = nrg + 1;
#line 689 "IB01MY.f"
	dcopy_(m, &dwork[ipg], &i__1, &dwork[jwork], &c__1);
#line 690 "IB01MY.f"
	i__1 = nrg + 1;
#line 690 "IB01MY.f"
	dcopy_(l, &dwork[ipy + mmnobr * nrg], &i__1, &dwork[jwork + *m], &
		c__1);

#line 693 "IB01MY.f"
	i__1 = *m + *l;
#line 693 "IB01MY.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 694 "IB01MY.f"
	    i__2 = *m + *l;
#line 694 "IB01MY.f"
	    iwork[i__] = idamax_(&i__2, &dwork[jwork], &c__1);
#line 695 "IB01MY.f"
	    dwork[jwork + iwork[i__] - 1] = 0.;
#line 696 "IB01MY.f"
/* L130: */
#line 696 "IB01MY.f"
	}

#line 698 "IB01MY.f"
	i__1 = *m + *l;
#line 698 "IB01MY.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 699 "IB01MY.f"
	    imax = iwork[i__];
#line 700 "IB01MY.f"
	    if (imax <= *m) {
#line 701 "IB01MY.f"
		icol = imax;
#line 702 "IB01MY.f"
	    } else {
#line 703 "IB01MY.f"
		icol = mmnobr - *m + imax;
#line 704 "IB01MY.f"
	    }
#line 705 "IB01MY.f"
	    beta = sqrt((d__1 = dwork[ipg + imax - 1 + (icol - 1) * nrg], abs(
		    d__1)));
#line 706 "IB01MY.f"
	    if (beta == 0.) {

/*              Error exit. */

#line 710 "IB01MY.f"
		*info = 1;
#line 711 "IB01MY.f"
		return 0;
#line 712 "IB01MY.f"
	    }
#line 713 "IB01MY.f"
	    d__1 = 1. / beta;
#line 713 "IB01MY.f"
	    dscal_(&nr, &d__1, &dwork[ipg + imax - 1], &nrg);
#line 714 "IB01MY.f"
	    dcopy_(&nr, &dwork[ipg + imax - 1], &nrg, &dwork[ing + imax - 1], 
		    &nrg);
#line 716 "IB01MY.f"
	    dwork[ipg + imax - 1 + (icol - 1) * nrg] = beta;
#line 717 "IB01MY.f"
	    dwork[ing + imax - 1 + (icol - 1) * nrg] = 0.;

#line 719 "IB01MY.f"
	    i__2 = *m + *l;
#line 719 "IB01MY.f"
	    for (j = i__ + 1; j <= i__2; ++j) {
#line 720 "IB01MY.f"
		dwork[ipg + iwork[j] - 1 + (icol - 1) * nrg] = 0.;
#line 721 "IB01MY.f"
/* L140: */
#line 721 "IB01MY.f"
	    }

#line 723 "IB01MY.f"
/* L150: */
#line 723 "IB01MY.f"
	}

/*        Compute the last two generators. */

#line 727 "IB01MY.f"
	if (! first) {

/*           For sequential processing, move the stored last negative */
/*           generator. */

#line 732 "IB01MY.f"
	    dcopy_(&nr, &dwork[nrg], &nrg, &dwork[ing + nrg - 1], &nrg);
#line 733 "IB01MY.f"
	}

#line 735 "IB01MY.f"
	jd = nrg;
#line 736 "IB01MY.f"
	if (*m > 0) {

#line 738 "IB01MY.f"
	    i__1 = *nsmp;
#line 738 "IB01MY.f"
	    for (j = ns; j <= i__1; ++j) {
#line 739 "IB01MY.f"
		dcopy_(m, &u[j + u_dim1], ldu, &dwork[jd], &nrg);
#line 740 "IB01MY.f"
		jd += mnrg;
#line 741 "IB01MY.f"
/* L160: */
#line 741 "IB01MY.f"
	    }

#line 743 "IB01MY.f"
	}

#line 745 "IB01MY.f"
	i__1 = *nsmp;
#line 745 "IB01MY.f"
	for (j = ns; j <= i__1; ++j) {
#line 746 "IB01MY.f"
	    dcopy_(l, &y[j + y_dim1], ldy, &dwork[jd], &nrg);
#line 747 "IB01MY.f"
	    jd += lnrg;
#line 748 "IB01MY.f"
/* L170: */
#line 748 "IB01MY.f"
	}

#line 750 "IB01MY.f"
	if (first) {
#line 751 "IB01MY.f"
	    if (*m > 0) {
#line 752 "IB01MY.f"
		dcopy_(m, dum, &c__0, &dwork[jd], &nrg);

#line 754 "IB01MY.f"
		i__1 = nobr21;
#line 754 "IB01MY.f"
		for (j = 1; j <= i__1; ++j) {
#line 755 "IB01MY.f"
		    jd += mnrg;
#line 756 "IB01MY.f"
		    dcopy_(m, &u[j + u_dim1], ldu, &dwork[jd], &nrg);
#line 757 "IB01MY.f"
/* L180: */
#line 757 "IB01MY.f"
		}

#line 759 "IB01MY.f"
		jd += mnrg;
#line 760 "IB01MY.f"
	    }
#line 761 "IB01MY.f"
	    dcopy_(l, dum, &c__0, &dwork[jd], &nrg);

#line 763 "IB01MY.f"
	    i__1 = nobr21;
#line 763 "IB01MY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 764 "IB01MY.f"
		jd += lnrg;
#line 765 "IB01MY.f"
		dcopy_(l, &y[j + y_dim1], ldy, &dwork[jd], &nrg);
#line 766 "IB01MY.f"
/* L190: */
#line 766 "IB01MY.f"
	    }

#line 768 "IB01MY.f"
	}

#line 770 "IB01MY.f"
	itau = jwork;
#line 771 "IB01MY.f"
	ipgc = ipg + mmnobr * nrg;

#line 773 "IB01MY.f"
	if (*m > 0) {

/*           Process the input part of the generators. */

#line 777 "IB01MY.f"
	    jwork = itau + *m;

/*           Reduce the first M columns of the matrix G1 of positive */
/*           generators to an upper triangular form. */
/*           Workspace: need   (M+L)*4*NOBR*(M+L+1)+2*M; */
/*                   prefer (M+L)*4*NOBR*(M+L+1)+M+M*NB. */

#line 784 "IB01MY.f"
	    ingc = ing;
#line 785 "IB01MY.f"
	    i__1 = *ldwork - jwork + 1;
#line 785 "IB01MY.f"
	    dgeqrf_(&nrg, m, &dwork[ipg], &nrg, &dwork[itau], &dwork[jwork], &
		    i__1, &ierr);
/* Computing MAX */
#line 787 "IB01MY.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 787 "IB01MY.f"
	    maxwrk = max(i__1,i__2);

/*           Workspace: need   (M+L)*4*NOBR*(M+L+1)+(M+L)*2*NOBR; */
/*                      prefer (M+L)*4*NOBR*(M+L+1)+M+ */
/*                                                 ((M+L)*2*NOBR-M)*NB. */

#line 793 "IB01MY.f"
	    i__1 = nr - *m;
#line 793 "IB01MY.f"
	    i__2 = *ldwork - jwork + 1;
#line 793 "IB01MY.f"
	    dormqr_("Left", "Transpose", &nrg, &i__1, m, &dwork[ipg], &nrg, &
		    dwork[itau], &dwork[ipg + mnrg], &nrg, &dwork[jwork], &
		    i__2, &ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 796 "IB01MY.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 796 "IB01MY.f"
	    maxwrk = max(i__1,i__2);

/*           Annihilate, column by column, the first M columns of the */
/*           matrix G2 of negative generators, using Householder */
/*           transformations and modified hyperbolic plane rotations. */
/*           In the DLARF calls, ITAU is a pointer to the workspace */
/*           array. */

#line 804 "IB01MY.f"
	    i__1 = *m;
#line 804 "IB01MY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 805 "IB01MY.f"
		dlarfg_(&nrg, &dwork[ingc], &dwork[ingc + 1], &c__1, &tau);
#line 806 "IB01MY.f"
		beta = dwork[ingc];
#line 807 "IB01MY.f"
		dwork[ingc] = 1.;
#line 808 "IB01MY.f"
		ingp = ingc + nrg;
#line 809 "IB01MY.f"
		i__2 = nr - j;
#line 809 "IB01MY.f"
		dlarf_("Left", &nrg, &i__2, &dwork[ingc], &c__1, &tau, &dwork[
			ingp], &nrg, &dwork[itau], (ftnlen)4);
#line 811 "IB01MY.f"
		dwork[ingc] = beta;

/*              Compute the coefficients of the modified hyperbolic */
/*              rotation. */

#line 816 "IB01MY.f"
		ma02fd_(&dwork[ipg + (j - 1) * (nrg + 1)], &dwork[ingc], &cs, 
			&sn, &ierr);
#line 818 "IB01MY.f"
		if (ierr != 0) {

/*                 Error return: the matrix H'*H is not (numerically) */
/*                 positive definite. */

#line 823 "IB01MY.f"
		    *info = 1;
#line 824 "IB01MY.f"
		    return 0;
#line 825 "IB01MY.f"
		}

#line 827 "IB01MY.f"
		i__2 = (nr - 1) * nrg;
#line 827 "IB01MY.f"
		i__3 = nrg;
#line 827 "IB01MY.f"
		for (i__ = j * nrg; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ 
			+= i__3) {
#line 828 "IB01MY.f"
		    dwork[ipg + j - 1 + i__] = (dwork[ipg + j - 1 + i__] - sn 
			    * dwork[ing + i__]) / cs;
#line 830 "IB01MY.f"
		    dwork[ing + i__] = -sn * dwork[ipg + j - 1 + i__] + cs * 
			    dwork[ing + i__];
#line 832 "IB01MY.f"
/* L200: */
#line 832 "IB01MY.f"
		}

#line 834 "IB01MY.f"
		ingc = ingp;
#line 835 "IB01MY.f"
/* L210: */
#line 835 "IB01MY.f"
	    }

/*           Save one block row of R, and shift the generators for the */
/*           calculation of the following row. */

#line 840 "IB01MY.f"
	    dlacpy_("Upper", m, &nr, &dwork[ipg], &nrg, &r__[r_offset], ldr, (
		    ftnlen)5);

#line 842 "IB01MY.f"
	    i__1 = mnrg;
#line 842 "IB01MY.f"
	    i__3 = -mnrg;
#line 842 "IB01MY.f"
	    for (i__ = (mmnobr - *m) * nrg; i__3 < 0 ? i__ >= i__1 : i__ <= 
		    i__1; i__ += i__3) {
#line 843 "IB01MY.f"
		dlacpy_("Full", m, m, &dwork[ipg + i__ - mnrg], &nrg, &dwork[
			ipg + i__], &nrg, (ftnlen)4);
#line 845 "IB01MY.f"
/* L220: */
#line 845 "IB01MY.f"
	    }

#line 847 "IB01MY.f"
	    i__3 = (mmnobr + *l) * nrg;
#line 847 "IB01MY.f"
	    i__1 = -lnrg;
#line 847 "IB01MY.f"
	    for (i__ = (nr - *l) * nrg; i__1 < 0 ? i__ >= i__3 : i__ <= i__3; 
		    i__ += i__1) {
#line 848 "IB01MY.f"
		dlacpy_("Full", m, l, &dwork[ipg + i__ - lnrg], &nrg, &dwork[
			ipg + i__], &nrg, (ftnlen)4);
#line 850 "IB01MY.f"
/* L230: */
#line 850 "IB01MY.f"
	    }

#line 852 "IB01MY.f"
	    dlaset_("Full", m, l, &c_b111, &c_b111, &dwork[ipgc], &nrg, (
		    ftnlen)4);

/*           Update the input part of generators using Schur algorithm. */
/*           Workspace: need   (M+L)*4*NOBR*(M+L+1)+2*NOBR*(M+L)-M. */

#line 857 "IB01MY.f"
	    jds = mnrg;
#line 858 "IB01MY.f"
	    icol = *m;

#line 860 "IB01MY.f"
	    i__1 = nobr2;
#line 860 "IB01MY.f"
	    for (k = 2; k <= i__1; ++k) {
#line 861 "IB01MY.f"
		i__3 = nr - icol - *m;
#line 861 "IB01MY.f"
		i__2 = *l + 1;
#line 861 "IB01MY.f"
		mb04od_("Full", m, &i__3, &i__2, &dwork[ipg + jds], &nrg, &
			dwork[ipy + jds], &nrg, &dwork[ipg + jds + mnrg], &
			nrg, &dwork[ipy + jds + mnrg], &nrg, &dwork[itau], &
			dwork[jwork], (ftnlen)4);

#line 867 "IB01MY.f"
		i__3 = *m;
#line 867 "IB01MY.f"
		for (j = 1; j <= i__3; ++j) {
#line 868 "IB01MY.f"
		    icj = icol + j;
#line 869 "IB01MY.f"
		    dlarfg_(&nrg, &dwork[ingc], &dwork[ingc + 1], &c__1, &tau)
			    ;
#line 870 "IB01MY.f"
		    beta = dwork[ingc];
#line 871 "IB01MY.f"
		    dwork[ingc] = 1.;
#line 872 "IB01MY.f"
		    ingp = ingc + nrg;
#line 873 "IB01MY.f"
		    i__2 = nr - icj;
#line 873 "IB01MY.f"
		    dlarf_("Left", &nrg, &i__2, &dwork[ingc], &c__1, &tau, &
			    dwork[ingp], &nrg, &dwork[itau], (ftnlen)4);
#line 875 "IB01MY.f"
		    dwork[ingc] = beta;

/*                 Compute the coefficients of the modified hyperbolic */
/*                 rotation. */

#line 880 "IB01MY.f"
		    ma02fd_(&dwork[ipg + j - 1 + (icj - 1) * nrg], &dwork[
			    ingc], &cs, &sn, &ierr);
#line 882 "IB01MY.f"
		    if (ierr != 0) {

/*                    Error return: the matrix H'*H is not (numerically) */
/*                    positive definite. */

#line 887 "IB01MY.f"
			*info = 1;
#line 888 "IB01MY.f"
			return 0;
#line 889 "IB01MY.f"
		    }

#line 891 "IB01MY.f"
		    i__2 = (nr - 1) * nrg;
#line 891 "IB01MY.f"
		    i__4 = nrg;
#line 891 "IB01MY.f"
		    for (i__ = icj * nrg; i__4 < 0 ? i__ >= i__2 : i__ <= 
			    i__2; i__ += i__4) {
#line 892 "IB01MY.f"
			dwork[ipg + j - 1 + i__] = (dwork[ipg + j - 1 + i__] 
				- sn * dwork[ing + i__]) / cs;
#line 894 "IB01MY.f"
			dwork[ing + i__] = -sn * dwork[ipg + j - 1 + i__] + 
				cs * dwork[ing + i__];
#line 896 "IB01MY.f"
/* L240: */
#line 896 "IB01MY.f"
		    }

#line 898 "IB01MY.f"
		    ingc = ingp;
#line 899 "IB01MY.f"
/* L250: */
#line 899 "IB01MY.f"
		}

/*              Save one block row of R, and shift the generators for the */
/*              calculation of the following row. */

#line 904 "IB01MY.f"
		i__3 = nr - icol;
#line 904 "IB01MY.f"
		dlacpy_("Upper", m, &i__3, &dwork[ipg + jds], &nrg, &r__[icol 
			+ 1 + (icol + 1) * r_dim1], ldr, (ftnlen)5);
#line 906 "IB01MY.f"
		icol += *m;

#line 908 "IB01MY.f"
		i__3 = icol * nrg;
#line 908 "IB01MY.f"
		i__4 = -mnrg;
#line 908 "IB01MY.f"
		for (i__ = (mmnobr - *m) * nrg; i__4 < 0 ? i__ >= i__3 : i__ 
			<= i__3; i__ += i__4) {
#line 909 "IB01MY.f"
		    dlacpy_("Full", m, m, &dwork[ipg + i__ - mnrg], &nrg, &
			    dwork[ipg + i__], &nrg, (ftnlen)4);
#line 911 "IB01MY.f"
/* L260: */
#line 911 "IB01MY.f"
		}

#line 913 "IB01MY.f"
		i__4 = (mmnobr + *l) * nrg;
#line 913 "IB01MY.f"
		i__3 = -lnrg;
#line 913 "IB01MY.f"
		for (i__ = (nr - *l) * nrg; i__3 < 0 ? i__ >= i__4 : i__ <= 
			i__4; i__ += i__3) {
#line 914 "IB01MY.f"
		    dlacpy_("Full", m, l, &dwork[ipg + i__ - lnrg], &nrg, &
			    dwork[ipg + i__], &nrg, (ftnlen)4);
#line 916 "IB01MY.f"
/* L270: */
#line 916 "IB01MY.f"
		}

#line 918 "IB01MY.f"
		dlaset_("Full", m, l, &c_b111, &c_b111, &dwork[ipgc], &nrg, (
			ftnlen)4);
#line 919 "IB01MY.f"
		jds += mnrg;
#line 920 "IB01MY.f"
/* L280: */
#line 920 "IB01MY.f"
	    }

#line 922 "IB01MY.f"
	}

/*        Process the output part of the generators. */

#line 926 "IB01MY.f"
	jwork = itau + *l;

/*        Reduce the first L columns of the submatrix */
/*        G1(1:M+L+1,2*M*NOBR+1:2*(M+L)*NOBR) to upper triangular form. */
/*        Workspace: need   (M+L)*4*NOBR*(M+L+1)+2*L; */
/*                   prefer (M+L)*4*NOBR*(M+L+1)+L+L*NB. */

#line 933 "IB01MY.f"
	ingc = ing + mmnobr * nrg;
#line 934 "IB01MY.f"
	i__1 = *ldwork - jwork + 1;
#line 934 "IB01MY.f"
	dgeqrf_(&nrg, l, &dwork[ipgc], &nrg, &dwork[itau], &dwork[jwork], &
		i__1, &ierr);
/* Computing MAX */
#line 936 "IB01MY.f"
	i__1 = maxwrk, i__3 = (integer) dwork[jwork] + jwork - 1;
#line 936 "IB01MY.f"
	maxwrk = max(i__1,i__3);

/*        Workspace: need   (M+L)*4*NOBR*(M+L+1)+L*2*NOBR; */
/*                   prefer (M+L)*4*NOBR*(M+L+1)+L+(L*2*NOBR-L)*NB. */

#line 941 "IB01MY.f"
	i__1 = llnobr - *l;
#line 941 "IB01MY.f"
	i__3 = *ldwork - jwork + 1;
#line 941 "IB01MY.f"
	dormqr_("Left", "Transpose", &nrg, &i__1, l, &dwork[ipgc], &nrg, &
		dwork[itau], &dwork[ipgc + lnrg], &nrg, &dwork[jwork], &i__3, 
		&ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 944 "IB01MY.f"
	i__1 = maxwrk, i__3 = (integer) dwork[jwork] + jwork - 1;
#line 944 "IB01MY.f"
	maxwrk = max(i__1,i__3);

/*        Annihilate, column by column, the first L columns of the */
/*        output part of the matrix G2 of negative generators, using */
/*        Householder transformations and modified hyperbolic rotations. */

#line 950 "IB01MY.f"
	i__1 = *l;
#line 950 "IB01MY.f"
	for (j = 1; j <= i__1; ++j) {
#line 951 "IB01MY.f"
	    dlarfg_(&nrg, &dwork[ingc], &dwork[ingc + 1], &c__1, &tau);
#line 952 "IB01MY.f"
	    beta = dwork[ingc];
#line 953 "IB01MY.f"
	    dwork[ingc] = 1.;
#line 954 "IB01MY.f"
	    ingp = ingc + nrg;
#line 955 "IB01MY.f"
	    i__3 = llnobr - j;
#line 955 "IB01MY.f"
	    dlarf_("Left", &nrg, &i__3, &dwork[ingc], &c__1, &tau, &dwork[
		    ingp], &nrg, &dwork[itau], (ftnlen)4);
#line 957 "IB01MY.f"
	    dwork[ingc] = beta;

/*           Compute the coefficients of the modified hyperbolic */
/*           rotation. */

#line 962 "IB01MY.f"
	    ma02fd_(&dwork[ipgc + (j - 1) * (nrg + 1)], &dwork[ingc], &cs, &
		    sn, &ierr);
#line 964 "IB01MY.f"
	    if (ierr != 0) {

/*              Error return: the matrix H'*H is not (numerically) */
/*              positive definite. */

#line 969 "IB01MY.f"
		*info = 1;
#line 970 "IB01MY.f"
		return 0;
#line 971 "IB01MY.f"
	    }

#line 973 "IB01MY.f"
	    i__3 = (nr - 1) * nrg;
#line 973 "IB01MY.f"
	    i__4 = nrg;
#line 973 "IB01MY.f"
	    for (i__ = (j + mmnobr) * nrg; i__4 < 0 ? i__ >= i__3 : i__ <= 
		    i__3; i__ += i__4) {
#line 974 "IB01MY.f"
		dwork[ipg + j - 1 + i__] = (dwork[ipg + j - 1 + i__] - sn * 
			dwork[ing + i__]) / cs;
#line 976 "IB01MY.f"
		dwork[ing + i__] = -sn * dwork[ipg + j - 1 + i__] + cs * 
			dwork[ing + i__];
#line 978 "IB01MY.f"
/* L290: */
#line 978 "IB01MY.f"
	    }

#line 980 "IB01MY.f"
	    ingc = ingp;
#line 981 "IB01MY.f"
/* L300: */
#line 981 "IB01MY.f"
	}

/*        Save one block row of R, and shift the generators for the */
/*        calculation of the following row. */

#line 986 "IB01MY.f"
	dlacpy_("Upper", l, &llnobr, &dwork[ipgc], &nrg, &r__[mmnobr + 1 + (
		mmnobr + 1) * r_dim1], ldr, (ftnlen)5);

#line 989 "IB01MY.f"
	i__1 = (mmnobr + *l) * nrg;
#line 989 "IB01MY.f"
	i__4 = -lnrg;
#line 989 "IB01MY.f"
	for (i__ = (nr - *l) * nrg; i__4 < 0 ? i__ >= i__1 : i__ <= i__1; i__ 
		+= i__4) {
#line 990 "IB01MY.f"
	    dlacpy_("Full", l, l, &dwork[ipg + i__ - lnrg], &nrg, &dwork[ipg 
		    + i__], &nrg, (ftnlen)4);
#line 992 "IB01MY.f"
/* L310: */
#line 992 "IB01MY.f"
	}

/*        Update the output part of generators using the Schur algorithm. */
/*        Workspace: need   (M+L)*4*NOBR*(M+L+1)+2*NOBR*L-L. */

#line 997 "IB01MY.f"
	jds = lnrg;
#line 998 "IB01MY.f"
	icol = *l;

#line 1000 "IB01MY.f"
	i__4 = nobr2;
#line 1000 "IB01MY.f"
	for (k = 2; k <= i__4; ++k) {
#line 1001 "IB01MY.f"
	    i__1 = llnobr - icol - *l;
#line 1001 "IB01MY.f"
	    i__3 = *m + 1;
#line 1001 "IB01MY.f"
	    mb04od_("Full", l, &i__1, &i__3, &dwork[ipgc + jds], &nrg, &dwork[
		    ipgc + *l + jds], &nrg, &dwork[ipgc + jds + lnrg], &nrg, &
		    dwork[ipgc + *l + jds + lnrg], &nrg, &dwork[itau], &dwork[
		    jwork], (ftnlen)4);

#line 1007 "IB01MY.f"
	    i__1 = *l;
#line 1007 "IB01MY.f"
	    for (j = 1; j <= i__1; ++j) {
#line 1008 "IB01MY.f"
		icj = icol + j;
#line 1009 "IB01MY.f"
		dlarfg_(&nrg, &dwork[ingc], &dwork[ingc + 1], &c__1, &tau);
#line 1010 "IB01MY.f"
		beta = dwork[ingc];
#line 1011 "IB01MY.f"
		dwork[ingc] = 1.;
#line 1012 "IB01MY.f"
		ingp = ingc + nrg;
#line 1013 "IB01MY.f"
		i__3 = llnobr - icj;
#line 1013 "IB01MY.f"
		dlarf_("Left", &nrg, &i__3, &dwork[ingc], &c__1, &tau, &dwork[
			ingp], &nrg, &dwork[itau], (ftnlen)4);
#line 1015 "IB01MY.f"
		dwork[ingc] = beta;

/*              Compute the coefficients of the modified hyperbolic */
/*              rotation. */

#line 1020 "IB01MY.f"
		ma02fd_(&dwork[ipgc + j - 1 + (icj - 1) * nrg], &dwork[ingc], 
			&cs, &sn, &ierr);
#line 1022 "IB01MY.f"
		if (ierr != 0) {

/*                 Error return: the matrix H'*H is not (numerically) */
/*                 positive definite. */

#line 1027 "IB01MY.f"
		    *info = 1;
#line 1028 "IB01MY.f"
		    return 0;
#line 1029 "IB01MY.f"
		}

#line 1031 "IB01MY.f"
		i__3 = (nr - 1) * nrg;
#line 1031 "IB01MY.f"
		i__2 = nrg;
#line 1031 "IB01MY.f"
		for (i__ = (icj + mmnobr) * nrg; i__2 < 0 ? i__ >= i__3 : i__ 
			<= i__3; i__ += i__2) {
#line 1032 "IB01MY.f"
		    dwork[ipg + j - 1 + i__] = (dwork[ipg + j - 1 + i__] - sn 
			    * dwork[ing + i__]) / cs;
#line 1034 "IB01MY.f"
		    dwork[ing + i__] = -sn * dwork[ipg + j - 1 + i__] + cs * 
			    dwork[ing + i__];
#line 1036 "IB01MY.f"
/* L320: */
#line 1036 "IB01MY.f"
		}

#line 1038 "IB01MY.f"
		ingc = ingp;
#line 1039 "IB01MY.f"
/* L330: */
#line 1039 "IB01MY.f"
	    }

/*           Save one block row of R, and shift the generators for the */
/*           calculation of the following row. */

#line 1044 "IB01MY.f"
	    i__1 = llnobr - icol;
#line 1044 "IB01MY.f"
	    dlacpy_("Upper", l, &i__1, &dwork[ipgc + jds], &nrg, &r__[mmnobr 
		    + icol + 1 + (mmnobr + icol + 1) * r_dim1], ldr, (ftnlen)
		    5);

#line 1047 "IB01MY.f"
	    i__1 = (mmnobr + icol) * nrg;
#line 1047 "IB01MY.f"
	    i__2 = -lnrg;
#line 1047 "IB01MY.f"
	    for (i__ = (nr - *l) * nrg; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; 
		    i__ += i__2) {
#line 1048 "IB01MY.f"
		dlacpy_("Full", l, l, &dwork[ipg + i__ - lnrg], &nrg, &dwork[
			ipg + i__], &nrg, (ftnlen)4);
#line 1050 "IB01MY.f"
/* L340: */
#line 1050 "IB01MY.f"
	    }

#line 1052 "IB01MY.f"
	    icol += *l;
#line 1053 "IB01MY.f"
	    jds += lnrg;
#line 1054 "IB01MY.f"
/* L350: */
#line 1054 "IB01MY.f"
	}

#line 1056 "IB01MY.f"
	if (moesp && *m > 0) {

/*           For the MOESP algorithm, interchange the past and future */
/*           input parts of the R factor, and compute the new R factor */
/*           using a specialized QR factorization.  A tailored fast */
/*           QR factorization for the MOESP algorithm could be slightly */
/*           more efficient. */

#line 1064 "IB01MY.f"
	    i__4 = mnobr;
#line 1064 "IB01MY.f"
	    for (j = 1; j <= i__4; ++j) {
#line 1065 "IB01MY.f"
		dswap_(&j, &r__[j * r_dim1 + 1], &c__1, &r__[(mnobr + j) * 
			r_dim1 + 1], &c__1);
#line 1066 "IB01MY.f"
		dcopy_(&mnobr, &r__[j + 1 + (mnobr + j) * r_dim1], &c__1, &
			r__[j + 1 + j * r_dim1], &c__1);
#line 1067 "IB01MY.f"
		i__2 = mmnobr - j;
#line 1067 "IB01MY.f"
		dcopy_(&i__2, dum, &c__0, &r__[j + 1 + (mnobr + j) * r_dim1], 
			&c__1);
#line 1068 "IB01MY.f"
/* L360: */
#line 1068 "IB01MY.f"
	    }

/*           Triangularize the first two block columns (using structure), */
/*           and apply the transformation to the corresponding part of */
/*           the remaining block columns. */
/*           Workspace: need 2*(M+L)*NOBR. */

#line 1075 "IB01MY.f"
	    itau = 1;
#line 1076 "IB01MY.f"
	    jwork = itau + mmnobr;
#line 1077 "IB01MY.f"
	    i__4 = mnobr - 1;
#line 1077 "IB01MY.f"
	    i__2 = *ldwork - jwork + 1;
#line 1077 "IB01MY.f"
	    mb04id_(&mmnobr, &mmnobr, &i__4, &llnobr, &r__[r_offset], ldr, &
		    r__[(mmnobr + 1) * r_dim1 + 1], ldr, &dwork[itau], &dwork[
		    jwork], &i__2, &ierr);
/* Computing MAX */
#line 1080 "IB01MY.f"
	    i__4 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 1080 "IB01MY.f"
	    maxwrk = max(i__4,i__2);
#line 1081 "IB01MY.f"
	}
#line 1082 "IB01MY.f"
    }

#line 1084 "IB01MY.f"
    nsmpsm = 0;
#line 1085 "IB01MY.f"
    icycle = 1;

/*     Return optimal workspace in  DWORK(1). */

#line 1089 "IB01MY.f"
    dwork[1] = (doublereal) maxwrk;
#line 1090 "IB01MY.f"
    maxwrk = 1;
#line 1091 "IB01MY.f"
    return 0;

/* *** Last line of IB01MY *** */
} /* ib01my_ */

