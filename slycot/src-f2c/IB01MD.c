#line 1 "IB01MD.f"
/* IB01MD.f -- translated by f2c (version 20100827).
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

#line 1 "IB01MD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b24 = -1.;
static doublereal c_b29 = 1.;
static doublereal c_b207 = 0.;
static integer c__0 = 0;

/* Subroutine */ int ib01md_(char *meth, char *alg, char *batch, char *conct, 
	integer *nobr, integer *m, integer *l, integer *nsmp, doublereal *u, 
	integer *ldu, doublereal *y, integer *ldy, doublereal *r__, integer *
	ldr, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	iwarn, integer *info, ftnlen meth_len, ftnlen alg_len, ftnlen 
	batch_len, ftnlen conct_len)
{
    /* System generated locals */
    integer r_dim1, r_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, id, jd, ii, nr, ns;
    static doublereal dum[1];
    static integer nsf;
    static doublereal upd;
    static integer inu, nsl, iny;
    extern /* Subroutine */ int dger_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer icol, ierr, itau, init;
    static logical last;
    static doublereal temp;
    static integer irev;
    static logical linr, n4sid;
    static integer nobr2;
    static logical chalg;
    extern /* Subroutine */ int mb04od_(char *, integer *, integer *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    ftnlen), dgemm_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ib01my_(char *, char *, char *, integer *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    static integer nobr21;
    static logical qralg;
    static integer initi, lnobr, mnobr;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical moesp;
    static integer lldrw, mldrw;
    static logical first;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    static integer nobrm1, ishft2;
    static logical onebch, connec;
    static integer icycle;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static logical fqralg;
    static integer ncycle, inicyc;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static integer nicycl;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    static integer lmnobr, mmnobr;
    static logical interm;
    extern /* Subroutine */ int dpotrf_(char *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static integer ishftu, nslast, ldrwrk, ishfty, minwrk, maxwrk, ldrwmx, 
	    nsmpsm;


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
/*     block Hankel matrices using input-output data.  The input-output */
/*     data can, optionally, be processed sequentially. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     METH    CHARACTER*1 */
/*             Specifies the subspace identification method to be used, */
/*             as follows: */
/*             = 'M':  MOESP  algorithm with past inputs and outputs; */
/*             = 'N':  N4SID  algorithm. */

/*     ALG     CHARACTER*1 */
/*             Specifies the algorithm for computing the triangular */
/*             factor R, as follows: */
/*             = 'C':  Cholesky algorithm applied to the correlation */
/*                     matrix of the input-output data; */
/*             = 'F':  Fast QR algorithm; */
/*             = 'Q':  QR algorithm applied to the concatenated block */
/*                     Hankel matrices. */

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
/*             (In the MOESP theory,  NOBR  should be larger than  n, */
/*             the estimated dimension of state vector.) */

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

/*     R       (output or input/output) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             On exit, if INFO = 0 and ALG = 'Q', or (ALG = 'C' or 'F', */
/*             and BATCH = 'L' or 'O'), the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of */
/*             this array contains the (current) upper triangular factor */
/*             R from the QR factorization of the concatenated block */
/*             Hankel matrices. The diagonal elements of R are positive */
/*             when the Cholesky algorithm was successfully used. */
/*             On exit, if ALG = 'C' and BATCH = 'F' or 'I', the leading */
/*             2*(M+L)*NOBR-by-2*(M+L)*NOBR upper triangular part of this */
/*             array contains the current upper triangular part of the */
/*             correlation matrix in sequential data processing. */
/*             If ALG = 'F' and BATCH = 'F' or 'I', the array R is not */
/*             referenced. */
/*             On entry, if ALG = 'C', or ALG = 'Q', and BATCH = 'I' or */
/*             'L', the leading  2*(M+L)*NOBR-by-2*(M+L)*NOBR  upper */
/*             triangular part of this array must contain the upper */
/*             triangular matrix R computed at the previous call of this */
/*             routine in sequential data processing. The array R need */
/*             not be set on entry if ALG = 'F' or if BATCH = 'F' or 'O'. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= 2*(M+L)*NOBR. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= M+L, if ALG = 'F'; */
/*             LIWORK >= 0,   if ALG = 'C' or 'Q'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -17,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */
/*             Let */
/*             k = 0,               if CONCT = 'N' and ALG = 'C' or 'Q'; */
/*             k = 2*NOBR-1,        if CONCT = 'C' and ALG = 'C' or 'Q'; */
/*             k = 2*NOBR*(M+L+1),  if CONCT = 'N' and ALG = 'F'; */
/*             k = 2*NOBR*(M+L+2),  if CONCT = 'C' and ALG = 'F'. */
/*             The first (M+L)*k elements of  DWORK  should be preserved */
/*             during successive calls of the routine with  BATCH = 'F' */
/*             or  'I',  till the final call with  BATCH = 'L'. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= (4*NOBR-2)*(M+L), if ALG = 'C', BATCH <> 'O' and */
/*                                     CONCT = 'C'; */
/*             LDWORK >= 1,            if ALG = 'C', BATCH = 'O' or */
/*                                     CONCT = 'N'; */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+3), if ALG = 'F', */
/*                                     BATCH <> 'O' and CONCT = 'C'; */
/*             LDWORK >= (M+L)*2*NOBR*(M+L+1), if ALG = 'F', */
/*                                     BATCH = 'F', 'I' and CONCT = 'N'; */
/*             LDWORK >= (M+L)*4*NOBR*(M+L+1)+(M+L)*2*NOBR, if ALG = 'F', */
/*                                     BATCH = 'L' and CONCT = 'N', or */
/*                                     BATCH = 'O'; */
/*             LDWORK >= 4*(M+L)*NOBR, if ALG = 'Q', BATCH = 'F' or 'O', */
/*                                     and LDR >= NS = NSMP - 2*NOBR + 1; */
/*             LDWORK >= 6*(M+L)*NOBR, if ALG = 'Q', BATCH = 'F' or 'O', */
/*                                     and LDR < NS, or BATCH = 'I' or */
/*                                     'L' and CONCT = 'N'; */
/*             LDWORK >= 4*(NOBR+1)*(M+L)*NOBR, if ALG = 'Q', BATCH = 'I' */
/*                                     or 'L' and CONCT = 'C'. */
/*             The workspace used for ALG = 'Q' is */
/*                       LDRWRK*2*(M+L)*NOBR + 4*(M+L)*NOBR, */
/*             where LDRWRK = LDWORK/(2*(M+L)*NOBR) - 2; recommended */
/*             value LDRWRK = NS, assuming a large enough cache size. */
/*             For good performance,  LDWORK  should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  the number of 100 cycles in sequential data */
/*                   processing has been exhausted without signaling */
/*                   that the last block of data was get; the cycle */
/*                   counter was reinitialized; */
/*             = 2:  a fast algorithm was requested (ALG = 'C' or 'F'), */
/*                   but it failed, and the QR algorithm was then used */
/*                   (non-sequential data processing). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  a fast algorithm was requested (ALG = 'C', or 'F') */
/*                   in sequential data processing, but it failed. The */
/*                   routine can be repeatedly called again using the */
/*                   standard QR algorithm. */

/*     METHOD */

/*     1) For non-sequential data processing using QR algorithm, a */
/*     t x 2(m+l)s  matrix H is constructed, where */

/*          H = [ Uf'         Up'      Y'      ],  for METH = 'M', */
/*                  s+1,2s,t    1,s,t   1,2s,t */

/*          H = [ U'       Y'      ],              for METH = 'N', */
/*                 1,2s,t   1,2s,t */

/*     and  Up     , Uf        , U      , and  Y        are block Hankel */
/*            1,s,t    s+1,2s,t   1,2s,t        1,2s,t */
/*     matrices defined in terms of the input and output data [3]. */
/*     A QR factorization is used to compress the data. */
/*     The fast QR algorithm uses a QR factorization which exploits */
/*     the block-Hankel structure. Actually, the Cholesky factor of H'*H */
/*     is computed. */

/*     2) For sequential data processing using QR algorithm, the QR */
/*     decomposition is done sequentially, by updating the upper */
/*     triangular factor  R.  This is also performed internally if the */
/*     workspace is not large enough to accommodate an entire batch. */

/*     3) For non-sequential or sequential data processing using */
/*     Cholesky algorithm, the correlation matrix of input-output data is */
/*     computed (sequentially, if requested), taking advantage of the */
/*     block Hankel structure [7].  Then, the Cholesky factor of the */
/*     correlation matrix is found, if possible. */

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

/*     [5] Peternell, K., Scherrer, W. and Deistler, M. */
/*         Statistical Analysis of Novel Subspace Identification Methods. */
/*         Signal Processing, 52, pp. 161-177, 1996. */

/*     [6] Sima, V. */
/*         Subspace-based Algorithms for Multivariable System */
/*         Identification. */
/*         Studies in Informatics and Control, 5, pp. 335-344, 1996. */

/*     [7] Sima, V. */
/*         Cholesky or QR Factorization for Data Compression in */
/*         Subspace-based Identification ? */
/*         Proceedings of the Second NICONET Workshop on ``Numerical */
/*         Control Software: SLICOT, a Useful Tool in Industry'', */
/*         December 3, 1999, INRIA Rocquencourt, France, pp. 75-80, 1999. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable (when QR algorithm is */
/*     used), reliable and efficient. The fast Cholesky or QR algorithms */
/*     are more efficient, but the accuracy could diminish by forming the */
/*     correlation matrix. */
/*                                        2 */
/*     The QR algorithm needs 0(t(2(m+l)s) ) floating point operations. */
/*                                           2              3 */
/*     The Cholesky algorithm needs 0(2t(m+l) s)+0((2(m+l)s) ) floating */
/*     point operations. */
/*                                          2           3 2 */
/*     The fast QR algorithm needs 0(2t(m+l) s)+0(4(m+l) s ) floating */
/*     point operations. */

/*     FURTHER COMMENTS */

/*     For ALG = 'Q', BATCH = 'O' and LDR < NS, or BATCH <> 'O', the */
/*     calculations could be rather inefficient if only minimal workspace */
/*     (see argument LDWORK) is provided. It is advisable to provide as */
/*     much workspace as possible. Almost optimal efficiency can be */
/*     obtained for  LDWORK = (NS+2)*(2*(M+L)*NOBR),  assuming that the */
/*     cache size is large enough to accommodate R, U, Y, and DWORK. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Aug. 1999. */

/*     REVISIONS */

/*     Feb. 2000, Aug. 2000, Feb. 2004. */

/*     KEYWORDS */

/*     Cholesky decomposition, Hankel matrix, identification methods, */
/*     multivariable systems, QR decomposition. */

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
/*        ICYCLE  is used to count the cycles for  BATCH = 'I'. It is */
/*                reinitialized at each MAXCYC cycles. */
/*        MAXWRK  is used to store the optimal workspace. */
/*        NSMPSM  is used to sum up the  NSMP  values for  BATCH <> 'O'. */
/*     .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 374 "IB01MD.f"
    /* Parameter adjustments */
#line 374 "IB01MD.f"
    u_dim1 = *ldu;
#line 374 "IB01MD.f"
    u_offset = 1 + u_dim1;
#line 374 "IB01MD.f"
    u -= u_offset;
#line 374 "IB01MD.f"
    y_dim1 = *ldy;
#line 374 "IB01MD.f"
    y_offset = 1 + y_dim1;
#line 374 "IB01MD.f"
    y -= y_offset;
#line 374 "IB01MD.f"
    r_dim1 = *ldr;
#line 374 "IB01MD.f"
    r_offset = 1 + r_dim1;
#line 374 "IB01MD.f"
    r__ -= r_offset;
#line 374 "IB01MD.f"
    --iwork;
#line 374 "IB01MD.f"
    --dwork;
#line 374 "IB01MD.f"

#line 374 "IB01MD.f"
    /* Function Body */
#line 374 "IB01MD.f"
    moesp = lsame_(meth, "M", (ftnlen)1, (ftnlen)1);
#line 375 "IB01MD.f"
    n4sid = lsame_(meth, "N", (ftnlen)1, (ftnlen)1);
#line 376 "IB01MD.f"
    fqralg = lsame_(alg, "F", (ftnlen)1, (ftnlen)1);
#line 377 "IB01MD.f"
    qralg = lsame_(alg, "Q", (ftnlen)1, (ftnlen)1);
#line 378 "IB01MD.f"
    chalg = lsame_(alg, "C", (ftnlen)1, (ftnlen)1);
#line 379 "IB01MD.f"
    onebch = lsame_(batch, "O", (ftnlen)1, (ftnlen)1);
#line 380 "IB01MD.f"
    first = lsame_(batch, "F", (ftnlen)1, (ftnlen)1) || onebch;
#line 381 "IB01MD.f"
    interm = lsame_(batch, "I", (ftnlen)1, (ftnlen)1);
#line 382 "IB01MD.f"
    last = lsame_(batch, "L", (ftnlen)1, (ftnlen)1) || onebch;
#line 383 "IB01MD.f"
    if (! onebch) {
#line 384 "IB01MD.f"
	connec = lsame_(conct, "C", (ftnlen)1, (ftnlen)1);
#line 385 "IB01MD.f"
    } else {
#line 386 "IB01MD.f"
	connec = FALSE_;
#line 387 "IB01MD.f"
    }

#line 389 "IB01MD.f"
    mnobr = *m * *nobr;
#line 390 "IB01MD.f"
    lnobr = *l * *nobr;
#line 391 "IB01MD.f"
    lmnobr = lnobr + mnobr;
#line 392 "IB01MD.f"
    mmnobr = mnobr + mnobr;
#line 393 "IB01MD.f"
    nobrm1 = *nobr - 1;
#line 394 "IB01MD.f"
    nobr21 = *nobr + nobrm1;
#line 395 "IB01MD.f"
    nobr2 = nobr21 + 1;
#line 396 "IB01MD.f"
    *iwarn = 0;
#line 397 "IB01MD.f"
    *info = 0;
#line 398 "IB01MD.f"
    ierr = 0;
#line 399 "IB01MD.f"
    if (first) {
#line 400 "IB01MD.f"
	icycle = 1;
#line 401 "IB01MD.f"
	maxwrk = 1;
#line 402 "IB01MD.f"
	nsmpsm = 0;
#line 403 "IB01MD.f"
    }
#line 404 "IB01MD.f"
    nsmpsm += *nsmp;
#line 405 "IB01MD.f"
    nr = lmnobr + lmnobr;

/*     Check the scalar input parameters. */

#line 409 "IB01MD.f"
    if (! (moesp || n4sid)) {
#line 410 "IB01MD.f"
	*info = -1;
#line 411 "IB01MD.f"
    } else if (! (fqralg || qralg || chalg)) {
#line 412 "IB01MD.f"
	*info = -2;
#line 413 "IB01MD.f"
    } else if (! (first || interm || last)) {
#line 414 "IB01MD.f"
	*info = -3;
#line 415 "IB01MD.f"
    } else if (! onebch) {
#line 416 "IB01MD.f"
	if (! (connec || lsame_(conct, "N", (ftnlen)1, (ftnlen)1))) {
#line 416 "IB01MD.f"
	    *info = -4;
#line 416 "IB01MD.f"
	}
#line 418 "IB01MD.f"
    }
#line 419 "IB01MD.f"
    if (*info == 0) {
#line 420 "IB01MD.f"
	if (*nobr <= 0) {
#line 421 "IB01MD.f"
	    *info = -5;
#line 422 "IB01MD.f"
	} else if (*m < 0) {
#line 423 "IB01MD.f"
	    *info = -6;
#line 424 "IB01MD.f"
	} else if (*l <= 0) {
#line 425 "IB01MD.f"
	    *info = -7;
#line 426 "IB01MD.f"
	} else if (*nsmp < nobr2 || last && nsmpsm < nr + nobr21) {
#line 428 "IB01MD.f"
	    *info = -8;
#line 429 "IB01MD.f"
	} else if (*ldu < 1 || *m > 0 && *ldu < *nsmp) {
#line 430 "IB01MD.f"
	    *info = -10;
#line 431 "IB01MD.f"
	} else if (*ldy < *nsmp) {
#line 432 "IB01MD.f"
	    *info = -12;
#line 433 "IB01MD.f"
	} else if (*ldr < nr) {
#line 434 "IB01MD.f"
	    *info = -14;
#line 435 "IB01MD.f"
	} else {

/*           Compute workspace. */
/*           (Note: Comments in the code beginning "Workspace:" describe */
/*           the minimal amount of workspace needed at that point in the */
/*           code, as well as the preferred amount for good performance. */
/*           NB refers to the optimal block size for the immediately */
/*           following subroutine, as returned by ILAENV.) */

#line 444 "IB01MD.f"
	    ns = *nsmp - nobr21;
#line 445 "IB01MD.f"
	    if (chalg) {
#line 446 "IB01MD.f"
		if (! onebch && connec) {
#line 447 "IB01MD.f"
		    minwrk = nr - *m - *l << 1;
#line 448 "IB01MD.f"
		} else {
#line 449 "IB01MD.f"
		    minwrk = 1;
#line 450 "IB01MD.f"
		}
#line 451 "IB01MD.f"
	    } else if (fqralg) {
#line 452 "IB01MD.f"
		if (! onebch && connec) {
#line 453 "IB01MD.f"
		    minwrk = nr * (*m + *l + 3);
#line 454 "IB01MD.f"
		} else if (first || interm) {
#line 455 "IB01MD.f"
		    minwrk = nr * (*m + *l + 1);
#line 456 "IB01MD.f"
		} else {
#line 457 "IB01MD.f"
		    minwrk = (nr << 1) * (*m + *l + 1) + nr;
#line 458 "IB01MD.f"
		}
#line 459 "IB01MD.f"
	    } else {
#line 460 "IB01MD.f"
		minwrk = nr << 1;
#line 461 "IB01MD.f"
		maxwrk = nr + nr * ilaenv_(&c__1, "DGEQRF", " ", &ns, &nr, &
			c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 463 "IB01MD.f"
		if (first) {
#line 464 "IB01MD.f"
		    if (*ldr < ns) {
#line 465 "IB01MD.f"
			minwrk += nr;
#line 466 "IB01MD.f"
			maxwrk = ns * nr + maxwrk;
#line 467 "IB01MD.f"
		    }
#line 468 "IB01MD.f"
		} else {
#line 469 "IB01MD.f"
		    if (connec) {
#line 470 "IB01MD.f"
			minwrk *= *nobr + 1;
#line 471 "IB01MD.f"
		    } else {
#line 472 "IB01MD.f"
			minwrk += nr;
#line 473 "IB01MD.f"
		    }
#line 474 "IB01MD.f"
		    maxwrk = ns * nr + maxwrk;
#line 475 "IB01MD.f"
		}
#line 476 "IB01MD.f"
	    }
#line 477 "IB01MD.f"
	    maxwrk = max(minwrk,maxwrk);

#line 479 "IB01MD.f"
	    if (*ldwork < minwrk) {
#line 480 "IB01MD.f"
		*info = -17;
#line 481 "IB01MD.f"
		dwork[1] = (doublereal) minwrk;
#line 482 "IB01MD.f"
	    }
#line 483 "IB01MD.f"
	}
#line 484 "IB01MD.f"
    }

/*     Return if there are illegal arguments. */

#line 488 "IB01MD.f"
    if (*info != 0) {
#line 489 "IB01MD.f"
	i__1 = -(*info);
#line 489 "IB01MD.f"
	xerbla_("IB01MD", &i__1, (ftnlen)6);
#line 490 "IB01MD.f"
	return 0;
#line 491 "IB01MD.f"
    }

#line 493 "IB01MD.f"
    if (chalg) {

/*        Compute the  R  factor from a Cholesky factorization of the */
/*        input-output data correlation matrix. */

/*        Set the parameters for constructing the correlations of the */
/*        current block. */

#line 501 "IB01MD.f"
	ldrwrk = (nobr2 << 1) - 2;
#line 502 "IB01MD.f"
	if (first) {
#line 503 "IB01MD.f"
	    upd = 0.;
#line 504 "IB01MD.f"
	} else {
#line 505 "IB01MD.f"
	    upd = 1.;
#line 506 "IB01MD.f"
	}

#line 508 "IB01MD.f"
	if (! first && connec) {

/*           Restore the saved (M+L)*(2*NOBR-1) "connection" elements of */
/*           U  and  Y  into their appropriate position in sequential */
/*           processing. The process is performed column-wise, in */
/*           reverse order, first for  Y  and then for  U. */
/*           Workspace: need   (4*NOBR-2)*(M+L). */

#line 516 "IB01MD.f"
	    irev = nr - *m - *l - nobr21 + 1;
#line 517 "IB01MD.f"
	    icol = (nr - *m - *l << 1) - ldrwrk + 1;

#line 519 "IB01MD.f"
	    i__1 = *m + *l;
#line 519 "IB01MD.f"
	    for (j = 2; j <= i__1; ++j) {
#line 520 "IB01MD.f"
		for (i__ = nobr21 - 1; i__ >= 0; --i__) {
#line 521 "IB01MD.f"
		    dwork[icol + i__] = dwork[irev + i__];
#line 522 "IB01MD.f"
/* L5: */
#line 522 "IB01MD.f"
		}
#line 523 "IB01MD.f"
		irev -= nobr21;
#line 524 "IB01MD.f"
		icol -= ldrwrk;
#line 525 "IB01MD.f"
/* L10: */
#line 525 "IB01MD.f"
	    }

#line 527 "IB01MD.f"
	    if (*m > 0) {
#line 527 "IB01MD.f"
		dlacpy_("Full", &nobr21, m, &u[u_offset], ldu, &dwork[nobr2], 
			&ldrwrk, (ftnlen)4);
#line 527 "IB01MD.f"
	    }
#line 530 "IB01MD.f"
	    dlacpy_("Full", &nobr21, l, &y[y_offset], ldy, &dwork[ldrwrk * *m 
		    + nobr2], &ldrwrk, (ftnlen)4);
#line 532 "IB01MD.f"
	}

#line 534 "IB01MD.f"
	if (*m > 0) {

/*           Let  Guu(i,j) = Guu0(i,j) + u_i*u_j' + u_(i+1)*u_(j+1)' + */
/*                                 ... + u_(i+NS-1)*u_(j+NS-1)', */
/*           where  u_i'  is the i-th row of  U,  j = 1 : 2s,  i = 1 : j, */
/*           NS = NSMP - 2s + 1,  and  Guu0(i,j)  is a zero matrix for */
/*           BATCH = 'O' or 'F', and it is the matrix Guu(i,j) computed */
/*           till the current block for BATCH = 'I' or 'L'. The matrix */
/*           Guu(i,j)  is  m-by-m,  and  Guu(j,j)  is symmetric. The */
/*           upper triangle of the U-U correlations,  Guu,  is computed */
/*           (or updated) column-wise in the array  R,  that is, in the */
/*           order  Guu(1,1),  Guu(1,2),  Guu(2,2),  ...,  Guu(2s,2s). */
/*           Only the submatrices of the first block-row are fully */
/*           computed (or updated). The remaining ones are determined */
/*           exploiting the block-Hankel structure, using the updating */
/*           formula */

/*           Guu(i+1,j+1) = Guu0(i+1,j+1) - Guu0(i,j) + Guu(i,j) + */
/*                                 u_(i+NS)*u_(j+NS)' - u_i*u_j'. */

#line 554 "IB01MD.f"
	    if (! first) {

/*              Subtract the contribution of the previous block of data */
/*              in sequential processing. The columns must be processed */
/*              in backward order. */

#line 560 "IB01MD.f"
		for (i__ = nobr21 * *m; i__ >= 1; --i__) {
#line 561 "IB01MD.f"
		    daxpy_(&i__, &c_b24, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
			    m + 1 + (*m + i__) * r_dim1], &c__1);
#line 562 "IB01MD.f"
/* L20: */
#line 562 "IB01MD.f"
		}

#line 564 "IB01MD.f"
	    }

/*           Compute/update  Guu(1,1). */

#line 568 "IB01MD.f"
	    if (! first && connec) {
#line 568 "IB01MD.f"
		dsyrk_("Upper", "Transpose", m, &nobr21, &c_b29, &dwork[1], &
			ldrwrk, &upd, &r__[r_offset], ldr, (ftnlen)5, (ftnlen)
			9);
#line 568 "IB01MD.f"
	    }
#line 571 "IB01MD.f"
	    dsyrk_("Upper", "Transpose", m, &ns, &c_b29, &u[u_offset], ldu, &
		    upd, &r__[r_offset], ldr, (ftnlen)5, (ftnlen)9);

#line 574 "IB01MD.f"
	    jd = 1;

#line 576 "IB01MD.f"
	    if (first || ! connec) {

#line 578 "IB01MD.f"
		i__1 = nobr2;
#line 578 "IB01MD.f"
		for (j = 2; j <= i__1; ++j) {
#line 579 "IB01MD.f"
		    jd += *m;
#line 580 "IB01MD.f"
		    id = *m + 1;

/*                 Compute/update  Guu(1,j). */

#line 584 "IB01MD.f"
		    dgemm_("Transpose", "NoTranspose", m, m, &ns, &c_b29, &u[
			    u_offset], ldu, &u[j + u_dim1], ldu, &upd, &r__[
			    jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);

/*                 Compute/update  Guu(2:j,j), exploiting the */
/*                 block-Hankel structure. */

#line 590 "IB01MD.f"
		    if (first) {

#line 592 "IB01MD.f"
			i__2 = jd - 1;
#line 592 "IB01MD.f"
			for (i__ = jd - *m; i__ <= i__2; ++i__) {
#line 593 "IB01MD.f"
			    dcopy_(&i__, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
				    m + 1 + (*m + i__) * r_dim1], &c__1);
#line 594 "IB01MD.f"
/* L30: */
#line 594 "IB01MD.f"
			}

#line 596 "IB01MD.f"
		    } else {

#line 598 "IB01MD.f"
			i__2 = jd - 1;
#line 598 "IB01MD.f"
			for (i__ = jd - *m; i__ <= i__2; ++i__) {
#line 599 "IB01MD.f"
			    daxpy_(&i__, &c_b29, &r__[i__ * r_dim1 + 1], &
				    c__1, &r__[*m + 1 + (*m + i__) * r_dim1], 
				    &c__1);
#line 600 "IB01MD.f"
/* L40: */
#line 600 "IB01MD.f"
			}

#line 602 "IB01MD.f"
		    }

#line 604 "IB01MD.f"
		    i__2 = j - 1;
#line 604 "IB01MD.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 605 "IB01MD.f"
			dger_(m, m, &c_b29, &u[ns + i__ - 1 + u_dim1], ldu, &
				u[ns + j - 1 + u_dim1], ldu, &r__[id + jd * 
				r_dim1], ldr);
#line 607 "IB01MD.f"
			dger_(m, m, &c_b24, &u[i__ - 1 + u_dim1], ldu, &u[j - 
				1 + u_dim1], ldu, &r__[id + jd * r_dim1], ldr)
				;
#line 609 "IB01MD.f"
			id += *m;
#line 610 "IB01MD.f"
/* L50: */
#line 610 "IB01MD.f"
		    }

#line 612 "IB01MD.f"
		    i__2 = *m;
#line 612 "IB01MD.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 613 "IB01MD.f"
			daxpy_(&i__, &u[ns + j - 1 + i__ * u_dim1], &u[ns + j 
				- 1 + u_dim1], ldu, &r__[jd + (jd + i__ - 1) *
				 r_dim1], &c__1);
#line 615 "IB01MD.f"
			d__1 = -u[j - 1 + i__ * u_dim1];
#line 615 "IB01MD.f"
			daxpy_(&i__, &d__1, &u[j - 1 + u_dim1], ldu, &r__[jd 
				+ (jd + i__ - 1) * r_dim1], &c__1);
#line 617 "IB01MD.f"
/* L60: */
#line 617 "IB01MD.f"
		    }

#line 619 "IB01MD.f"
/* L70: */
#line 619 "IB01MD.f"
		}

#line 621 "IB01MD.f"
	    } else {

#line 623 "IB01MD.f"
		i__1 = nobr2;
#line 623 "IB01MD.f"
		for (j = 2; j <= i__1; ++j) {
#line 624 "IB01MD.f"
		    jd += *m;
#line 625 "IB01MD.f"
		    id = *m + 1;

/*                 Compute/update  Guu(1,j)  for sequential processing */
/*                 with connected blocks. */

#line 630 "IB01MD.f"
		    dgemm_("Transpose", "NoTranspose", m, m, &nobr21, &c_b29, 
			    &dwork[1], &ldrwrk, &dwork[j], &ldrwrk, &upd, &
			    r__[jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);
#line 633 "IB01MD.f"
		    dgemm_("Transpose", "NoTranspose", m, m, &ns, &c_b29, &u[
			    u_offset], ldu, &u[j + u_dim1], ldu, &c_b29, &r__[
			    jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);

/*                 Compute/update  Guu(2:j,j)  for sequential processing */
/*                 with connected blocks, exploiting the block-Hankel */
/*                 structure. */

#line 640 "IB01MD.f"
		    if (first) {

#line 642 "IB01MD.f"
			i__2 = jd - 1;
#line 642 "IB01MD.f"
			for (i__ = jd - *m; i__ <= i__2; ++i__) {
#line 643 "IB01MD.f"
			    dcopy_(&i__, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
				    m + 1 + (*m + i__) * r_dim1], &c__1);
#line 644 "IB01MD.f"
/* L80: */
#line 644 "IB01MD.f"
			}

#line 646 "IB01MD.f"
		    } else {

#line 648 "IB01MD.f"
			i__2 = jd - 1;
#line 648 "IB01MD.f"
			for (i__ = jd - *m; i__ <= i__2; ++i__) {
#line 649 "IB01MD.f"
			    daxpy_(&i__, &c_b29, &r__[i__ * r_dim1 + 1], &
				    c__1, &r__[*m + 1 + (*m + i__) * r_dim1], 
				    &c__1);
#line 650 "IB01MD.f"
/* L90: */
#line 650 "IB01MD.f"
			}

#line 652 "IB01MD.f"
		    }

#line 654 "IB01MD.f"
		    i__2 = j - 1;
#line 654 "IB01MD.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 655 "IB01MD.f"
			dger_(m, m, &c_b29, &u[ns + i__ - 1 + u_dim1], ldu, &
				u[ns + j - 1 + u_dim1], ldu, &r__[id + jd * 
				r_dim1], ldr);
#line 657 "IB01MD.f"
			dger_(m, m, &c_b24, &dwork[i__ - 1], &ldrwrk, &dwork[
				j - 1], &ldrwrk, &r__[id + jd * r_dim1], ldr);
#line 659 "IB01MD.f"
			id += *m;
#line 660 "IB01MD.f"
/* L100: */
#line 660 "IB01MD.f"
		    }

#line 662 "IB01MD.f"
		    i__2 = *m;
#line 662 "IB01MD.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 663 "IB01MD.f"
			daxpy_(&i__, &u[ns + j - 1 + i__ * u_dim1], &u[ns + j 
				- 1 + u_dim1], ldu, &r__[jd + (jd + i__ - 1) *
				 r_dim1], &c__1);
#line 665 "IB01MD.f"
			d__1 = -dwork[(i__ - 1) * ldrwrk + j - 1];
#line 665 "IB01MD.f"
			daxpy_(&i__, &d__1, &dwork[j - 1], &ldrwrk, &r__[jd + 
				(jd + i__ - 1) * r_dim1], &c__1);
#line 667 "IB01MD.f"
/* L110: */
#line 667 "IB01MD.f"
		    }

#line 669 "IB01MD.f"
/* L120: */
#line 669 "IB01MD.f"
		}

#line 671 "IB01MD.f"
	    }

#line 673 "IB01MD.f"
	    if (last && moesp) {

/*              Interchange past and future parts for MOESP algorithm. */
/*              (Only the upper triangular parts are interchanged, and */
/*              the (1,2) part is transposed in-situ.) */

#line 679 "IB01MD.f"
		temp = r__[r_dim1 + 1];
#line 680 "IB01MD.f"
		r__[r_dim1 + 1] = r__[mnobr + 1 + (mnobr + 1) * r_dim1];
#line 681 "IB01MD.f"
		r__[mnobr + 1 + (mnobr + 1) * r_dim1] = temp;

#line 683 "IB01MD.f"
		i__1 = mnobr;
#line 683 "IB01MD.f"
		for (j = 2; j <= i__1; ++j) {
#line 684 "IB01MD.f"
		    dswap_(&j, &r__[j * r_dim1 + 1], &c__1, &r__[mnobr + 1 + (
			    mnobr + j) * r_dim1], &c__1);
#line 685 "IB01MD.f"
		    i__2 = j - 1;
#line 685 "IB01MD.f"
		    dswap_(&i__2, &r__[(mnobr + j) * r_dim1 + 1], &c__1, &r__[
			    j + (mnobr + 1) * r_dim1], ldr);
#line 686 "IB01MD.f"
/* L130: */
#line 686 "IB01MD.f"
		}

#line 688 "IB01MD.f"
	    }

/*           Let  Guy(i,j) = Guy0(i,j) + u_i*y_j' + u_(i+1)*y_(j+1)' + */
/*                                 ... + u_(i+NS-1)*y_(j+NS-1)', */
/*           where  u_i'  is the i-th row of  U,  y_j'  is the j-th row */
/*           of  Y,  j = 1 : 2s,  i = 1 : 2s,  NS = NSMP - 2s + 1,  and */
/*           Guy0(i,j)  is a zero matrix for  BATCH = 'O' or 'F', and it */
/*           is the matrix Guy(i,j) computed till the current block for */
/*           BATCH = 'I' or 'L'.  Guy(i,j) is m-by-L. The U-Y */
/*           correlations,  Guy,  are computed (or updated) column-wise */
/*           in the array  R. Only the submatrices of the first block- */
/*           column and block-row are fully computed (or updated). The */
/*           remaining ones are determined exploiting the block-Hankel */
/*           structure, using the updating formula */

/*           Guy(i+1,j+1) = Guy0(i+1,j+1) - Guy0(i,j) + Guy(i,j) + */
/*                                 u_(i+NS)*y(j+NS)' - u_i*y_j'. */

#line 706 "IB01MD.f"
	    ii = mmnobr - *m;
#line 707 "IB01MD.f"
	    if (! first) {

/*              Subtract the contribution of the previous block of data */
/*              in sequential processing. The columns must be processed */
/*              in backward order. */

#line 713 "IB01MD.f"
		i__1 = mmnobr + 1;
#line 713 "IB01MD.f"
		for (i__ = nr - *l; i__ >= i__1; --i__) {
#line 714 "IB01MD.f"
		    daxpy_(&ii, &c_b24, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
			    m + 1 + (*l + i__) * r_dim1], &c__1);
#line 715 "IB01MD.f"
/* L140: */
#line 715 "IB01MD.f"
		}

#line 717 "IB01MD.f"
	    }

/*           Compute/update the first block-column of  Guy,  Guy(i,1). */

#line 721 "IB01MD.f"
	    if (first || ! connec) {

#line 723 "IB01MD.f"
		i__1 = nobr2;
#line 723 "IB01MD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 724 "IB01MD.f"
		    dgemm_("Transpose", "NoTranspose", m, l, &ns, &c_b29, &u[
			    i__ + u_dim1], ldu, &y[y_offset], ldy, &upd, &r__[
			    (i__ - 1) * *m + 1 + (mmnobr + 1) * r_dim1], ldr, 
			    (ftnlen)9, (ftnlen)11);
#line 727 "IB01MD.f"
/* L150: */
#line 727 "IB01MD.f"
		}

#line 729 "IB01MD.f"
	    } else {

#line 731 "IB01MD.f"
		i__1 = nobr2;
#line 731 "IB01MD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 732 "IB01MD.f"
		    dgemm_("Transpose", "NoTranspose", m, l, &nobr21, &c_b29, 
			    &dwork[i__], &ldrwrk, &dwork[ldrwrk * *m + 1], &
			    ldrwrk, &upd, &r__[(i__ - 1) * *m + 1 + (mmnobr + 
			    1) * r_dim1], ldr, (ftnlen)9, (ftnlen)11);
#line 735 "IB01MD.f"
		    dgemm_("Transpose", "NoTranspose", m, l, &ns, &c_b29, &u[
			    i__ + u_dim1], ldu, &y[y_offset], ldy, &c_b29, &
			    r__[(i__ - 1) * *m + 1 + (mmnobr + 1) * r_dim1], 
			    ldr, (ftnlen)9, (ftnlen)11);
#line 738 "IB01MD.f"
/* L160: */
#line 738 "IB01MD.f"
		}

#line 740 "IB01MD.f"
	    }

#line 742 "IB01MD.f"
	    jd = mmnobr + 1;

#line 744 "IB01MD.f"
	    if (first || ! connec) {

#line 746 "IB01MD.f"
		i__1 = nobr2;
#line 746 "IB01MD.f"
		for (j = 2; j <= i__1; ++j) {
#line 747 "IB01MD.f"
		    jd += *l;
#line 748 "IB01MD.f"
		    id = *m + 1;

/*                 Compute/update  Guy(1,j). */

#line 752 "IB01MD.f"
		    dgemm_("Transpose", "NoTranspose", m, l, &ns, &c_b29, &u[
			    u_offset], ldu, &y[j + y_dim1], ldy, &upd, &r__[
			    jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);

/*                 Compute/update  Guy(2:2*s,j), exploiting the */
/*                 block-Hankel structure. */

#line 758 "IB01MD.f"
		    if (first) {

#line 760 "IB01MD.f"
			i__2 = jd - 1;
#line 760 "IB01MD.f"
			for (i__ = jd - *l; i__ <= i__2; ++i__) {
#line 761 "IB01MD.f"
			    dcopy_(&ii, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
				    m + 1 + (*l + i__) * r_dim1], &c__1);
#line 762 "IB01MD.f"
/* L170: */
#line 762 "IB01MD.f"
			}

#line 764 "IB01MD.f"
		    } else {

#line 766 "IB01MD.f"
			i__2 = jd - 1;
#line 766 "IB01MD.f"
			for (i__ = jd - *l; i__ <= i__2; ++i__) {
#line 767 "IB01MD.f"
			    daxpy_(&ii, &c_b29, &r__[i__ * r_dim1 + 1], &c__1,
				     &r__[*m + 1 + (*l + i__) * r_dim1], &
				    c__1);
#line 768 "IB01MD.f"
/* L180: */
#line 768 "IB01MD.f"
			}

#line 770 "IB01MD.f"
		    }

#line 772 "IB01MD.f"
		    i__2 = nobr2;
#line 772 "IB01MD.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 773 "IB01MD.f"
			dger_(m, l, &c_b29, &u[ns + i__ - 1 + u_dim1], ldu, &
				y[ns + j - 1 + y_dim1], ldy, &r__[id + jd * 
				r_dim1], ldr);
#line 775 "IB01MD.f"
			dger_(m, l, &c_b24, &u[i__ - 1 + u_dim1], ldu, &y[j - 
				1 + y_dim1], ldy, &r__[id + jd * r_dim1], ldr)
				;
#line 777 "IB01MD.f"
			id += *m;
#line 778 "IB01MD.f"
/* L190: */
#line 778 "IB01MD.f"
		    }

#line 780 "IB01MD.f"
/* L200: */
#line 780 "IB01MD.f"
		}

#line 782 "IB01MD.f"
	    } else {

#line 784 "IB01MD.f"
		i__1 = nobr2;
#line 784 "IB01MD.f"
		for (j = 2; j <= i__1; ++j) {
#line 785 "IB01MD.f"
		    jd += *l;
#line 786 "IB01MD.f"
		    id = *m + 1;

/*                 Compute/update  Guy(1,j)  for sequential processing */
/*                 with connected blocks. */

#line 791 "IB01MD.f"
		    dgemm_("Transpose", "NoTranspose", m, l, &nobr21, &c_b29, 
			    &dwork[1], &ldrwrk, &dwork[ldrwrk * *m + j], &
			    ldrwrk, &upd, &r__[jd * r_dim1 + 1], ldr, (ftnlen)
			    9, (ftnlen)11);
#line 794 "IB01MD.f"
		    dgemm_("Transpose", "NoTranspose", m, l, &ns, &c_b29, &u[
			    u_offset], ldu, &y[j + y_dim1], ldy, &c_b29, &r__[
			    jd * r_dim1 + 1], ldr, (ftnlen)9, (ftnlen)11);

/*                 Compute/update  Guy(2:2*s,j)  for sequential */
/*                 processing with connected blocks, exploiting the */
/*                 block-Hankel structure. */

#line 801 "IB01MD.f"
		    if (first) {

#line 803 "IB01MD.f"
			i__2 = jd - 1;
#line 803 "IB01MD.f"
			for (i__ = jd - *l; i__ <= i__2; ++i__) {
#line 804 "IB01MD.f"
			    dcopy_(&ii, &r__[i__ * r_dim1 + 1], &c__1, &r__[*
				    m + 1 + (*l + i__) * r_dim1], &c__1);
#line 805 "IB01MD.f"
/* L210: */
#line 805 "IB01MD.f"
			}

#line 807 "IB01MD.f"
		    } else {

#line 809 "IB01MD.f"
			i__2 = jd - 1;
#line 809 "IB01MD.f"
			for (i__ = jd - *l; i__ <= i__2; ++i__) {
#line 810 "IB01MD.f"
			    daxpy_(&ii, &c_b29, &r__[i__ * r_dim1 + 1], &c__1,
				     &r__[*m + 1 + (*l + i__) * r_dim1], &
				    c__1);
#line 811 "IB01MD.f"
/* L220: */
#line 811 "IB01MD.f"
			}

#line 813 "IB01MD.f"
		    }

#line 815 "IB01MD.f"
		    i__2 = nobr2;
#line 815 "IB01MD.f"
		    for (i__ = 2; i__ <= i__2; ++i__) {
#line 816 "IB01MD.f"
			dger_(m, l, &c_b29, &u[ns + i__ - 1 + u_dim1], ldu, &
				y[ns + j - 1 + y_dim1], ldy, &r__[id + jd * 
				r_dim1], ldr);
#line 818 "IB01MD.f"
			dger_(m, l, &c_b24, &dwork[i__ - 1], &ldrwrk, &dwork[
				ldrwrk * *m + j - 1], &ldrwrk, &r__[id + jd * 
				r_dim1], ldr);
#line 821 "IB01MD.f"
			id += *m;
#line 822 "IB01MD.f"
/* L230: */
#line 822 "IB01MD.f"
		    }

#line 824 "IB01MD.f"
/* L240: */
#line 824 "IB01MD.f"
		}

#line 826 "IB01MD.f"
	    }

#line 828 "IB01MD.f"
	    if (last && moesp) {

/*              Interchange past and future parts of U-Y correlations */
/*              for MOESP algorithm. */

#line 833 "IB01MD.f"
		i__1 = nr;
#line 833 "IB01MD.f"
		for (j = mmnobr + 1; j <= i__1; ++j) {
#line 834 "IB01MD.f"
		    dswap_(&mnobr, &r__[j * r_dim1 + 1], &c__1, &r__[mnobr + 
			    1 + j * r_dim1], &c__1);
#line 835 "IB01MD.f"
/* L250: */
#line 835 "IB01MD.f"
		}

#line 837 "IB01MD.f"
	    }
#line 838 "IB01MD.f"
	}

/*        Let  Gyy(i,j) = Gyy0(i,j) + y_i*y_i' + y_(i+1)*y_(i+1)' + ... + */
/*                                    y_(i+NS-1)*y_(i+NS-1)', */
/*        where  y_i'  is the i-th row of  Y,  j = 1 : 2s,  i = 1 : j, */
/*        NS = NSMP - 2s + 1,  and  Gyy0(i,j)  is a zero matrix for */
/*        BATCH = 'O' or 'F', and it is the matrix Gyy(i,j) computed till */
/*        the current block for BATCH = 'I' or 'L'.  Gyy(i,j) is L-by-L, */
/*        and  Gyy(j,j)  is symmetric. The upper triangle of the Y-Y */
/*        correlations,  Gyy,  is computed (or updated) column-wise in */
/*        the corresponding part of the array  R,  that is, in the order */
/*        Gyy(1,1),  Gyy(1,2),  Gyy(2,2),  ...,  Gyy(2s,2s).  Only the */
/*        submatrices of the first block-row are fully computed (or */
/*        updated). The remaining ones are determined exploiting the */
/*        block-Hankel structure, using the updating formula */

/*        Gyy(i+1,j+1) = Gyy0(i+1,j+1) - Gyy0(i,j) + Gyy(i,j) + */
/*                              y_(i+NS)*y_(j+NS)' - y_i*y_j'. */

#line 857 "IB01MD.f"
	jd = mmnobr + 1;

#line 859 "IB01MD.f"
	if (! first) {

/*           Subtract the contribution of the previous block of data */
/*           in sequential processing. The columns must be processed in */
/*           backward order. */

#line 865 "IB01MD.f"
	    i__1 = mmnobr + 1;
#line 865 "IB01MD.f"
	    for (i__ = nr - *l; i__ >= i__1; --i__) {
#line 866 "IB01MD.f"
		i__2 = i__ - mmnobr;
#line 866 "IB01MD.f"
		daxpy_(&i__2, &c_b24, &r__[jd + i__ * r_dim1], &c__1, &r__[jd 
			+ *l + (*l + i__) * r_dim1], &c__1);
#line 867 "IB01MD.f"
/* L260: */
#line 867 "IB01MD.f"
	    }

#line 869 "IB01MD.f"
	}

/*        Compute/update  Gyy(1,1). */

#line 873 "IB01MD.f"
	if (! first && connec) {
#line 873 "IB01MD.f"
	    dsyrk_("Upper", "Transpose", l, &nobr21, &c_b29, &dwork[ldrwrk * *
		    m + 1], &ldrwrk, &upd, &r__[jd + jd * r_dim1], ldr, (
		    ftnlen)5, (ftnlen)9);
#line 873 "IB01MD.f"
	}
#line 876 "IB01MD.f"
	dsyrk_("Upper", "Transpose", l, &ns, &c_b29, &y[y_offset], ldy, &upd, 
		&r__[jd + jd * r_dim1], ldr, (ftnlen)5, (ftnlen)9);

#line 879 "IB01MD.f"
	if (first || ! connec) {

#line 881 "IB01MD.f"
	    i__1 = nobr2;
#line 881 "IB01MD.f"
	    for (j = 2; j <= i__1; ++j) {
#line 882 "IB01MD.f"
		jd += *l;
#line 883 "IB01MD.f"
		id = mmnobr + *l + 1;

/*              Compute/update  Gyy(1,j). */

#line 887 "IB01MD.f"
		dgemm_("Transpose", "NoTranspose", l, l, &ns, &c_b29, &y[
			y_offset], ldy, &y[j + y_dim1], ldy, &upd, &r__[
			mmnobr + 1 + jd * r_dim1], ldr, (ftnlen)9, (ftnlen)11)
			;

/*              Compute/update  Gyy(2:j,j), exploiting the block-Hankel */
/*              structure. */

#line 893 "IB01MD.f"
		if (first) {

#line 895 "IB01MD.f"
		    i__2 = jd - 1;
#line 895 "IB01MD.f"
		    for (i__ = jd - *l; i__ <= i__2; ++i__) {
#line 896 "IB01MD.f"
			i__3 = i__ - mmnobr;
#line 896 "IB01MD.f"
			dcopy_(&i__3, &r__[mmnobr + 1 + i__ * r_dim1], &c__1, 
				&r__[mmnobr + *l + 1 + (*l + i__) * r_dim1], &
				c__1);
#line 898 "IB01MD.f"
/* L270: */
#line 898 "IB01MD.f"
		    }

#line 900 "IB01MD.f"
		} else {

#line 902 "IB01MD.f"
		    i__2 = jd - 1;
#line 902 "IB01MD.f"
		    for (i__ = jd - *l; i__ <= i__2; ++i__) {
#line 903 "IB01MD.f"
			i__3 = i__ - mmnobr;
#line 903 "IB01MD.f"
			daxpy_(&i__3, &c_b29, &r__[mmnobr + 1 + i__ * r_dim1],
				 &c__1, &r__[mmnobr + *l + 1 + (*l + i__) * 
				r_dim1], &c__1);
#line 905 "IB01MD.f"
/* L280: */
#line 905 "IB01MD.f"
		    }

#line 907 "IB01MD.f"
		}

#line 909 "IB01MD.f"
		i__2 = j - 1;
#line 909 "IB01MD.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 910 "IB01MD.f"
		    dger_(l, l, &c_b29, &y[ns + i__ - 1 + y_dim1], ldy, &y[ns 
			    + j - 1 + y_dim1], ldy, &r__[id + jd * r_dim1], 
			    ldr);
#line 912 "IB01MD.f"
		    dger_(l, l, &c_b24, &y[i__ - 1 + y_dim1], ldy, &y[j - 1 + 
			    y_dim1], ldy, &r__[id + jd * r_dim1], ldr);
#line 914 "IB01MD.f"
		    id += *l;
#line 915 "IB01MD.f"
/* L290: */
#line 915 "IB01MD.f"
		}

#line 917 "IB01MD.f"
		i__2 = *l;
#line 917 "IB01MD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 918 "IB01MD.f"
		    daxpy_(&i__, &y[ns + j - 1 + i__ * y_dim1], &y[ns + j - 1 
			    + y_dim1], ldy, &r__[jd + (jd + i__ - 1) * r_dim1]
			    , &c__1);
#line 920 "IB01MD.f"
		    d__1 = -y[j - 1 + i__ * y_dim1];
#line 920 "IB01MD.f"
		    daxpy_(&i__, &d__1, &y[j - 1 + y_dim1], ldy, &r__[jd + (
			    jd + i__ - 1) * r_dim1], &c__1);
#line 922 "IB01MD.f"
/* L300: */
#line 922 "IB01MD.f"
		}

#line 924 "IB01MD.f"
/* L310: */
#line 924 "IB01MD.f"
	    }

#line 926 "IB01MD.f"
	} else {

#line 928 "IB01MD.f"
	    i__1 = nobr2;
#line 928 "IB01MD.f"
	    for (j = 2; j <= i__1; ++j) {
#line 929 "IB01MD.f"
		jd += *l;
#line 930 "IB01MD.f"
		id = mmnobr + *l + 1;

/*              Compute/update  Gyy(1,j)  for sequential processing with */
/*              connected blocks. */

#line 935 "IB01MD.f"
		dgemm_("Transpose", "NoTranspose", l, l, &nobr21, &c_b29, &
			dwork[ldrwrk * *m + 1], &ldrwrk, &dwork[ldrwrk * *m + 
			j], &ldrwrk, &upd, &r__[mmnobr + 1 + jd * r_dim1], 
			ldr, (ftnlen)9, (ftnlen)11);
#line 939 "IB01MD.f"
		dgemm_("Transpose", "NoTranspose", l, l, &ns, &c_b29, &y[
			y_offset], ldy, &y[j + y_dim1], ldy, &c_b29, &r__[
			mmnobr + 1 + jd * r_dim1], ldr, (ftnlen)9, (ftnlen)11)
			;

/*              Compute/update  Gyy(2:j,j)  for sequential processing */
/*              with connected blocks, exploiting the block-Hankel */
/*              structure. */

#line 946 "IB01MD.f"
		if (first) {

#line 948 "IB01MD.f"
		    i__2 = jd - 1;
#line 948 "IB01MD.f"
		    for (i__ = jd - *l; i__ <= i__2; ++i__) {
#line 949 "IB01MD.f"
			i__3 = i__ - mmnobr;
#line 949 "IB01MD.f"
			dcopy_(&i__3, &r__[mmnobr + 1 + i__ * r_dim1], &c__1, 
				&r__[mmnobr + *l + 1 + (*l + i__) * r_dim1], &
				c__1);
#line 951 "IB01MD.f"
/* L320: */
#line 951 "IB01MD.f"
		    }

#line 953 "IB01MD.f"
		} else {

#line 955 "IB01MD.f"
		    i__2 = jd - 1;
#line 955 "IB01MD.f"
		    for (i__ = jd - *l; i__ <= i__2; ++i__) {
#line 956 "IB01MD.f"
			i__3 = i__ - mmnobr;
#line 956 "IB01MD.f"
			daxpy_(&i__3, &c_b29, &r__[mmnobr + 1 + i__ * r_dim1],
				 &c__1, &r__[mmnobr + *l + 1 + (*l + i__) * 
				r_dim1], &c__1);
#line 958 "IB01MD.f"
/* L330: */
#line 958 "IB01MD.f"
		    }

#line 960 "IB01MD.f"
		}

#line 962 "IB01MD.f"
		i__2 = j - 1;
#line 962 "IB01MD.f"
		for (i__ = 2; i__ <= i__2; ++i__) {
#line 963 "IB01MD.f"
		    dger_(l, l, &c_b29, &y[ns + i__ - 1 + y_dim1], ldy, &y[ns 
			    + j - 1 + y_dim1], ldy, &r__[id + jd * r_dim1], 
			    ldr);
#line 965 "IB01MD.f"
		    dger_(l, l, &c_b24, &dwork[ldrwrk * *m + i__ - 1], &
			    ldrwrk, &dwork[ldrwrk * *m + j - 1], &ldrwrk, &
			    r__[id + jd * r_dim1], ldr);
#line 968 "IB01MD.f"
		    id += *l;
#line 969 "IB01MD.f"
/* L340: */
#line 969 "IB01MD.f"
		}

#line 971 "IB01MD.f"
		i__2 = *l;
#line 971 "IB01MD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 972 "IB01MD.f"
		    daxpy_(&i__, &y[ns + j - 1 + i__ * y_dim1], &y[ns + j - 1 
			    + y_dim1], ldy, &r__[jd + (jd + i__ - 1) * r_dim1]
			    , &c__1);
#line 974 "IB01MD.f"
		    d__1 = -dwork[ldrwrk * (*m + i__ - 1) + j - 1];
#line 974 "IB01MD.f"
		    daxpy_(&i__, &d__1, &dwork[ldrwrk * *m + j - 1], &ldrwrk, 
			    &r__[jd + (jd + i__ - 1) * r_dim1], &c__1);
#line 977 "IB01MD.f"
/* L350: */
#line 977 "IB01MD.f"
		}

#line 979 "IB01MD.f"
/* L360: */
#line 979 "IB01MD.f"
	    }

#line 981 "IB01MD.f"
	}

#line 983 "IB01MD.f"
	if (! last) {
#line 984 "IB01MD.f"
	    if (connec) {

/*              For sequential processing with connected data blocks, */
/*              save the remaining ("connection") elements of  U  and  Y */
/*              in the first  (M+L)*(2*NOBR-1)  locations of  DWORK. */

#line 990 "IB01MD.f"
		if (*m > 0) {
#line 990 "IB01MD.f"
		    dlacpy_("Full", &nobr21, m, &u[ns + 1 + u_dim1], ldu, &
			    dwork[1], &nobr21, (ftnlen)4);
#line 990 "IB01MD.f"
		}
#line 993 "IB01MD.f"
		dlacpy_("Full", &nobr21, l, &y[ns + 1 + y_dim1], ldy, &dwork[
			mmnobr - *m + 1], &nobr21, (ftnlen)4);
#line 995 "IB01MD.f"
	    }

/*           Return to get new data. */

#line 999 "IB01MD.f"
	    ++icycle;
#line 1000 "IB01MD.f"
	    if (icycle > 100) {
#line 1000 "IB01MD.f"
		*iwarn = 1;
#line 1000 "IB01MD.f"
	    }
#line 1002 "IB01MD.f"
	    return 0;

#line 1004 "IB01MD.f"
	} else {

/*           Try to compute the Cholesky factor of the correlation */
/*           matrix. */

#line 1009 "IB01MD.f"
	    dpotrf_("Upper", &nr, &r__[r_offset], ldr, &ierr, (ftnlen)5);
#line 1010 "IB01MD.f"
	    goto L370;
#line 1011 "IB01MD.f"
	}
#line 1012 "IB01MD.f"
    } else if (fqralg) {

/*        Compute the  R  factor from a fast QR factorization of the */
/*        input-output data correlation matrix. */

#line 1017 "IB01MD.f"
	ib01my_(meth, batch, conct, nobr, m, l, nsmp, &u[u_offset], ldu, &y[
		y_offset], ldy, &r__[r_offset], ldr, &iwork[1], &dwork[1], 
		ldwork, iwarn, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
#line 1020 "IB01MD.f"
	if (! last) {
#line 1020 "IB01MD.f"
	    return 0;
#line 1020 "IB01MD.f"
	}
/* Computing MAX */
#line 1022 "IB01MD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[1];
#line 1022 "IB01MD.f"
	maxwrk = max(i__1,i__2);
#line 1023 "IB01MD.f"
    }

#line 1025 "IB01MD.f"
L370:

#line 1027 "IB01MD.f"
    if (ierr != 0) {

/*        Error return from a fast factorization algorithm of the */
/*        input-output data correlation matrix. */

#line 1032 "IB01MD.f"
	if (onebch) {
#line 1033 "IB01MD.f"
	    qralg = TRUE_;
#line 1034 "IB01MD.f"
	    *iwarn = 2;
#line 1035 "IB01MD.f"
	    minwrk = nr << 1;
#line 1036 "IB01MD.f"
	    maxwrk = nr + nr * ilaenv_(&c__1, "DGEQRF", " ", &ns, &nr, &c_n1, 
		    &c_n1, (ftnlen)6, (ftnlen)1);
#line 1038 "IB01MD.f"
	    if (*ldr < ns) {
#line 1039 "IB01MD.f"
		minwrk += nr;
#line 1040 "IB01MD.f"
		maxwrk = ns * nr + maxwrk;
#line 1041 "IB01MD.f"
	    }
#line 1042 "IB01MD.f"
	    maxwrk = max(minwrk,maxwrk);

#line 1044 "IB01MD.f"
	    if (*ldwork < minwrk) {
#line 1045 "IB01MD.f"
		*info = -17;

/*              Return: Not enough workspace. */

#line 1049 "IB01MD.f"
		dwork[1] = (doublereal) minwrk;
#line 1050 "IB01MD.f"
		i__1 = -(*info);
#line 1050 "IB01MD.f"
		xerbla_("IB01MD", &i__1, (ftnlen)6);
#line 1051 "IB01MD.f"
		return 0;
#line 1052 "IB01MD.f"
	    }
#line 1053 "IB01MD.f"
	} else {
#line 1054 "IB01MD.f"
	    *info = 1;
#line 1055 "IB01MD.f"
	    return 0;
#line 1056 "IB01MD.f"
	}
#line 1057 "IB01MD.f"
    }

#line 1059 "IB01MD.f"
    if (qralg) {

/*        Compute the  R  factor from a QR factorization of the matrix  H */
/*        of concatenated block Hankel matrices. */

/*        Construct the matrix  H. */

/*        Set the parameters for constructing the current segment of the */
/*        Hankel matrix, taking the available memory space into account. */
/*        INITI+1 points to the beginning rows of  U  and  Y  from which */
/*                data are taken when NCYCLE > 1 inner cycles are needed, */
/*                or for sequential processing with connected blocks. */
/*        LDRWMX is the number of rows that can fit in the working space. */
/*        LDRWRK is the actual number of rows processed in this space. */
/*        NSLAST is the number of samples to be processed at the last */
/*               inner cycle. */

#line 1076 "IB01MD.f"
	initi = 0;
#line 1077 "IB01MD.f"
	ldrwmx = *ldwork / nr - 2;
#line 1078 "IB01MD.f"
	ncycle = 1;
#line 1079 "IB01MD.f"
	nslast = *nsmp;
#line 1080 "IB01MD.f"
	linr = FALSE_;
#line 1081 "IB01MD.f"
	if (first) {
#line 1082 "IB01MD.f"
	    linr = *ldr >= ns;
#line 1083 "IB01MD.f"
	    ldrwrk = ns;
#line 1084 "IB01MD.f"
	} else if (connec) {
#line 1085 "IB01MD.f"
	    ldrwrk = *nsmp;
#line 1086 "IB01MD.f"
	} else {
#line 1087 "IB01MD.f"
	    ldrwrk = ns;
#line 1088 "IB01MD.f"
	}
#line 1089 "IB01MD.f"
	inicyc = 1;

#line 1091 "IB01MD.f"
	if (! linr) {
#line 1092 "IB01MD.f"
	    if (ldrwmx < ldrwrk) {

/*              Not enough working space for doing a single inner cycle. */
/*              NCYCLE inner cycles are to be performed for the current */
/*              data block using the working space. */

#line 1098 "IB01MD.f"
		ncycle = ldrwrk / ldrwmx;
#line 1099 "IB01MD.f"
		nslast = ldrwrk % ldrwmx;
#line 1100 "IB01MD.f"
		if (nslast != 0) {
#line 1101 "IB01MD.f"
		    ++ncycle;
#line 1102 "IB01MD.f"
		} else {
#line 1103 "IB01MD.f"
		    nslast = ldrwmx;
#line 1104 "IB01MD.f"
		}
#line 1105 "IB01MD.f"
		ldrwrk = ldrwmx;
#line 1106 "IB01MD.f"
		ns = ldrwrk;
#line 1107 "IB01MD.f"
		if (first) {
#line 1107 "IB01MD.f"
		    inicyc = 2;
#line 1107 "IB01MD.f"
		}
#line 1108 "IB01MD.f"
	    }
#line 1109 "IB01MD.f"
	    mldrw = *m * ldrwrk;
#line 1110 "IB01MD.f"
	    lldrw = *l * ldrwrk;
#line 1111 "IB01MD.f"
	    inu = mldrw * *nobr + 1;
#line 1112 "IB01MD.f"
	    iny = mldrw * nobr2 + 1;
#line 1113 "IB01MD.f"
	}

/*        Process the data given at the current call. */

#line 1117 "IB01MD.f"
	if (! first && connec) {

/*           Restore the saved (M+L)*(2*NOBR-1) "connection" elements of */
/*           U  and  Y  into their appropriate position in sequential */
/*           processing. The process is performed column-wise, in */
/*           reverse order, first for  Y  and then for  U. */

#line 1124 "IB01MD.f"
	    irev = nr - *m - *l - nobr21 + 1;
#line 1125 "IB01MD.f"
	    icol = iny + lldrw - ldrwrk;

#line 1127 "IB01MD.f"
	    i__1 = *l;
#line 1127 "IB01MD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 1128 "IB01MD.f"
		for (i__ = nobr21 - 1; i__ >= 0; --i__) {
#line 1129 "IB01MD.f"
		    dwork[icol + i__] = dwork[irev + i__];
#line 1130 "IB01MD.f"
/* L375: */
#line 1130 "IB01MD.f"
		}
#line 1131 "IB01MD.f"
		irev -= nobr21;
#line 1132 "IB01MD.f"
		icol -= ldrwrk;
#line 1133 "IB01MD.f"
/* L380: */
#line 1133 "IB01MD.f"
	    }

#line 1135 "IB01MD.f"
	    if (moesp) {
#line 1136 "IB01MD.f"
		icol = inu + mldrw - ldrwrk;
#line 1137 "IB01MD.f"
	    } else {
#line 1138 "IB01MD.f"
		icol = mldrw - ldrwrk + 1;
#line 1139 "IB01MD.f"
	    }

#line 1141 "IB01MD.f"
	    i__1 = *m;
#line 1141 "IB01MD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 1142 "IB01MD.f"
		for (i__ = nobr21 - 1; i__ >= 0; --i__) {
#line 1143 "IB01MD.f"
		    dwork[icol + i__] = dwork[irev + i__];
#line 1144 "IB01MD.f"
/* L385: */
#line 1144 "IB01MD.f"
		}
#line 1145 "IB01MD.f"
		irev -= nobr21;
#line 1146 "IB01MD.f"
		icol -= ldrwrk;
#line 1147 "IB01MD.f"
/* L390: */
#line 1147 "IB01MD.f"
	    }

#line 1149 "IB01MD.f"
	    if (moesp) {
#line 1149 "IB01MD.f"
		dlacpy_("Full", &nobrm1, m, &dwork[inu + *nobr], &ldrwrk, &
			dwork[1], &ldrwrk, (ftnlen)4);
#line 1149 "IB01MD.f"
	    }
#line 1152 "IB01MD.f"
	}

/*        Data compression using QR factorization. */

#line 1156 "IB01MD.f"
	if (first) {

/*           Non-sequential data processing or first block in */
/*           sequential data processing: */
/*           Use the general QR factorization algorithm. */

#line 1162 "IB01MD.f"
	    if (linr) {

/*              Put the input-output data in the array  R. */

#line 1166 "IB01MD.f"
		if (*m > 0) {
#line 1167 "IB01MD.f"
		    if (moesp) {

#line 1169 "IB01MD.f"
			i__1 = *nobr;
#line 1169 "IB01MD.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 1170 "IB01MD.f"
			    dlacpy_("Full", &ns, m, &u[*nobr + i__ + u_dim1], 
				    ldu, &r__[(*m * (i__ - 1) + 1) * r_dim1 + 
				    1], ldr, (ftnlen)4);
#line 1172 "IB01MD.f"
/* L400: */
#line 1172 "IB01MD.f"
			}

#line 1174 "IB01MD.f"
			i__1 = *nobr;
#line 1174 "IB01MD.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 1175 "IB01MD.f"
			    dlacpy_("Full", &ns, m, &u[i__ + u_dim1], ldu, &
				    r__[(mnobr + *m * (i__ - 1) + 1) * r_dim1 
				    + 1], ldr, (ftnlen)4);
#line 1177 "IB01MD.f"
/* L410: */
#line 1177 "IB01MD.f"
			}

#line 1179 "IB01MD.f"
		    } else {

#line 1181 "IB01MD.f"
			i__1 = nobr2;
#line 1181 "IB01MD.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 1182 "IB01MD.f"
			    dlacpy_("Full", &ns, m, &u[i__ + u_dim1], ldu, &
				    r__[(*m * (i__ - 1) + 1) * r_dim1 + 1], 
				    ldr, (ftnlen)4);
#line 1184 "IB01MD.f"
/* L420: */
#line 1184 "IB01MD.f"
			}

#line 1186 "IB01MD.f"
		    }
#line 1187 "IB01MD.f"
		}

#line 1189 "IB01MD.f"
		i__1 = nobr2;
#line 1189 "IB01MD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 1190 "IB01MD.f"
		    dlacpy_("Full", &ns, l, &y[i__ + y_dim1], ldy, &r__[(
			    mmnobr + *l * (i__ - 1) + 1) * r_dim1 + 1], ldr, (
			    ftnlen)4);
#line 1192 "IB01MD.f"
/* L430: */
#line 1192 "IB01MD.f"
		}

/*              Workspace: need   4*(M+L)*NOBR, */
/*                         prefer 2*(M+L)*NOBR+2*(M+L)*NOBR*NB. */

#line 1197 "IB01MD.f"
		itau = 1;
#line 1198 "IB01MD.f"
		jwork = itau + nr;
#line 1199 "IB01MD.f"
		i__1 = *ldwork - jwork + 1;
#line 1199 "IB01MD.f"
		dgeqrf_(&ns, &nr, &r__[r_offset], ldr, &dwork[itau], &dwork[
			jwork], &i__1, &ierr);
#line 1201 "IB01MD.f"
	    } else {

/*              Put the input-output data in the array  DWORK. */

#line 1205 "IB01MD.f"
		if (*m > 0) {
#line 1206 "IB01MD.f"
		    ishftu = 1;
#line 1207 "IB01MD.f"
		    if (moesp) {
#line 1208 "IB01MD.f"
			ishft2 = inu;

#line 1210 "IB01MD.f"
			i__1 = *nobr;
#line 1210 "IB01MD.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 1211 "IB01MD.f"
			    dlacpy_("Full", &ns, m, &u[*nobr + i__ + u_dim1], 
				    ldu, &dwork[ishftu], &ldrwrk, (ftnlen)4);
#line 1213 "IB01MD.f"
			    ishftu += mldrw;
#line 1214 "IB01MD.f"
/* L440: */
#line 1214 "IB01MD.f"
			}

#line 1216 "IB01MD.f"
			i__1 = *nobr;
#line 1216 "IB01MD.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 1217 "IB01MD.f"
			    dlacpy_("Full", &ns, m, &u[i__ + u_dim1], ldu, &
				    dwork[ishft2], &ldrwrk, (ftnlen)4);
#line 1219 "IB01MD.f"
			    ishft2 += mldrw;
#line 1220 "IB01MD.f"
/* L450: */
#line 1220 "IB01MD.f"
			}

#line 1222 "IB01MD.f"
		    } else {

#line 1224 "IB01MD.f"
			i__1 = nobr2;
#line 1224 "IB01MD.f"
			for (i__ = 1; i__ <= i__1; ++i__) {
#line 1225 "IB01MD.f"
			    dlacpy_("Full", &ns, m, &u[i__ + u_dim1], ldu, &
				    dwork[ishftu], &ldrwrk, (ftnlen)4);
#line 1227 "IB01MD.f"
			    ishftu += mldrw;
#line 1228 "IB01MD.f"
/* L460: */
#line 1228 "IB01MD.f"
			}

#line 1230 "IB01MD.f"
		    }
#line 1231 "IB01MD.f"
		}

#line 1233 "IB01MD.f"
		ishfty = iny;

#line 1235 "IB01MD.f"
		i__1 = nobr2;
#line 1235 "IB01MD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 1236 "IB01MD.f"
		    dlacpy_("Full", &ns, l, &y[i__ + y_dim1], ldy, &dwork[
			    ishfty], &ldrwrk, (ftnlen)4);
#line 1238 "IB01MD.f"
		    ishfty += lldrw;
#line 1239 "IB01MD.f"
/* L470: */
#line 1239 "IB01MD.f"
		}

/*              Workspace: need   2*(M+L)*NOBR + 4*(M+L)*NOBR, */
/*                         prefer NS*2*(M+L)*NOBR + 2*(M+L)*NOBR */
/*                                                + 2*(M+L)*NOBR*NB, */
/*                         used LDRWRK*2*(M+L)*NOBR + 4*(M+L)*NOBR, */
/*                         where  NS = NSMP - 2*NOBR + 1, */
/*                            LDRWRK = min(NS, LDWORK/(2*(M+L)*NOBR)-2). */

#line 1248 "IB01MD.f"
		itau = ldrwrk * nr + 1;
#line 1249 "IB01MD.f"
		jwork = itau + nr;
#line 1250 "IB01MD.f"
		i__1 = *ldwork - jwork + 1;
#line 1250 "IB01MD.f"
		dgeqrf_(&ns, &nr, &dwork[1], &ldrwrk, &dwork[itau], &dwork[
			jwork], &i__1, &ierr);
#line 1252 "IB01MD.f"
		i__1 = min(ns,nr);
#line 1252 "IB01MD.f"
		dlacpy_("Upper ", &i__1, &nr, &dwork[1], &ldrwrk, &r__[
			r_offset], ldr, (ftnlen)6);
#line 1254 "IB01MD.f"
	    }

#line 1256 "IB01MD.f"
	    if (ns < nr) {
#line 1256 "IB01MD.f"
		i__1 = nr - ns;
#line 1256 "IB01MD.f"
		i__2 = nr - ns;
#line 1256 "IB01MD.f"
		dlaset_("Upper ", &i__1, &i__2, &c_b207, &c_b207, &r__[ns + 1 
			+ (ns + 1) * r_dim1], ldr, (ftnlen)6);
#line 1256 "IB01MD.f"
	    }
#line 1259 "IB01MD.f"
	    initi += ns;
#line 1260 "IB01MD.f"
	}

#line 1262 "IB01MD.f"
	if (ncycle > 1 || ! first) {

/*           Remaining segments of the first data block or */
/*           remaining segments/blocks in sequential data processing: */
/*           Use a structure-exploiting QR factorization algorithm. */

#line 1268 "IB01MD.f"
	    nsl = ldrwrk;
#line 1269 "IB01MD.f"
	    if (! connec) {
#line 1269 "IB01MD.f"
		nsl = ns;
#line 1269 "IB01MD.f"
	    }
#line 1270 "IB01MD.f"
	    itau = ldrwrk * nr + 1;
#line 1271 "IB01MD.f"
	    jwork = itau + nr;

#line 1273 "IB01MD.f"
	    i__1 = ncycle;
#line 1273 "IB01MD.f"
	    for (nicycl = inicyc; nicycl <= i__1; ++nicycl) {

/*              INIT  denotes the beginning row where new data are put. */

#line 1277 "IB01MD.f"
		if (connec && nicycl == 1) {
#line 1278 "IB01MD.f"
		    init = nobr2;
#line 1279 "IB01MD.f"
		} else {
#line 1280 "IB01MD.f"
		    init = 1;
#line 1281 "IB01MD.f"
		}
#line 1282 "IB01MD.f"
		if (ncycle > 1 && nicycl == ncycle) {

/*                 Last samples in the last data segment of a block. */

#line 1286 "IB01MD.f"
		    ns = nslast;
#line 1287 "IB01MD.f"
		    nsl = nslast;
#line 1288 "IB01MD.f"
		}

/*              Put the input-output data in the array  DWORK. */

#line 1292 "IB01MD.f"
		nsf = ns;
#line 1293 "IB01MD.f"
		if (init > 1 && ncycle > 1) {
#line 1293 "IB01MD.f"
		    nsf -= nobr21;
#line 1293 "IB01MD.f"
		}
#line 1294 "IB01MD.f"
		if (*m > 0) {
#line 1295 "IB01MD.f"
		    ishftu = init;

#line 1297 "IB01MD.f"
		    if (moesp) {
#line 1298 "IB01MD.f"
			ishft2 = init + inu - 1;

#line 1300 "IB01MD.f"
			i__2 = *nobr;
#line 1300 "IB01MD.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 1301 "IB01MD.f"
			    dlacpy_("Full", &nsf, m, &u[initi + *nobr + i__ + 
				    u_dim1], ldu, &dwork[ishftu], &ldrwrk, (
				    ftnlen)4);
#line 1303 "IB01MD.f"
			    ishftu += mldrw;
#line 1304 "IB01MD.f"
/* L480: */
#line 1304 "IB01MD.f"
			}

#line 1306 "IB01MD.f"
			i__2 = *nobr;
#line 1306 "IB01MD.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 1307 "IB01MD.f"
			    dlacpy_("Full", &nsf, m, &u[initi + i__ + u_dim1],
				     ldu, &dwork[ishft2], &ldrwrk, (ftnlen)4);
#line 1309 "IB01MD.f"
			    ishft2 += mldrw;
#line 1310 "IB01MD.f"
/* L490: */
#line 1310 "IB01MD.f"
			}

#line 1312 "IB01MD.f"
		    } else {

#line 1314 "IB01MD.f"
			i__2 = nobr2;
#line 1314 "IB01MD.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 1315 "IB01MD.f"
			    dlacpy_("Full", &nsf, m, &u[initi + i__ + u_dim1],
				     ldu, &dwork[ishftu], &ldrwrk, (ftnlen)4);
#line 1317 "IB01MD.f"
			    ishftu += mldrw;
#line 1318 "IB01MD.f"
/* L500: */
#line 1318 "IB01MD.f"
			}

#line 1320 "IB01MD.f"
		    }
#line 1321 "IB01MD.f"
		}

#line 1323 "IB01MD.f"
		ishfty = init + iny - 1;

#line 1325 "IB01MD.f"
		i__2 = nobr2;
#line 1325 "IB01MD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1326 "IB01MD.f"
		    dlacpy_("Full", &nsf, l, &y[initi + i__ + y_dim1], ldy, &
			    dwork[ishfty], &ldrwrk, (ftnlen)4);
#line 1328 "IB01MD.f"
		    ishfty += lldrw;
#line 1329 "IB01MD.f"
/* L510: */
#line 1329 "IB01MD.f"
		}

#line 1331 "IB01MD.f"
		if (init > 1) {

/*                 Prepare the connection to the previous block of data */
/*                 in sequential processing. */

#line 1336 "IB01MD.f"
		    if (moesp && *m > 0) {
#line 1336 "IB01MD.f"
			dlacpy_("Full", nobr, m, &u[u_offset], ldu, &dwork[*
				nobr], &ldrwrk, (ftnlen)4);
#line 1336 "IB01MD.f"
		    }

/*                 Shift the elements from the connection to the previous */
/*                 block of data in sequential processing. */

#line 1343 "IB01MD.f"
		    if (*m > 0) {
#line 1344 "IB01MD.f"
			ishftu = mldrw + 1;

#line 1346 "IB01MD.f"
			if (moesp) {
#line 1347 "IB01MD.f"
			    ishft2 = mldrw + inu;

#line 1349 "IB01MD.f"
			    i__2 = nobrm1;
#line 1349 "IB01MD.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1350 "IB01MD.f"
				dlacpy_("Full", &nobr21, m, &dwork[ishftu - 
					mldrw + 1], &ldrwrk, &dwork[ishftu], &
					ldrwrk, (ftnlen)4);
#line 1353 "IB01MD.f"
				ishftu += mldrw;
#line 1354 "IB01MD.f"
/* L520: */
#line 1354 "IB01MD.f"
			    }

#line 1356 "IB01MD.f"
			    i__2 = nobrm1;
#line 1356 "IB01MD.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1357 "IB01MD.f"
				dlacpy_("Full", &nobr21, m, &dwork[ishft2 - 
					mldrw + 1], &ldrwrk, &dwork[ishft2], &
					ldrwrk, (ftnlen)4);
#line 1360 "IB01MD.f"
				ishft2 += mldrw;
#line 1361 "IB01MD.f"
/* L530: */
#line 1361 "IB01MD.f"
			    }

#line 1363 "IB01MD.f"
			} else {

#line 1365 "IB01MD.f"
			    i__2 = nobr21;
#line 1365 "IB01MD.f"
			    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1366 "IB01MD.f"
				dlacpy_("Full", &nobr21, m, &dwork[ishftu - 
					mldrw + 1], &ldrwrk, &dwork[ishftu], &
					ldrwrk, (ftnlen)4);
#line 1369 "IB01MD.f"
				ishftu += mldrw;
#line 1370 "IB01MD.f"
/* L540: */
#line 1370 "IB01MD.f"
			    }

#line 1372 "IB01MD.f"
			}
#line 1373 "IB01MD.f"
		    }

#line 1375 "IB01MD.f"
		    ishfty = lldrw + iny;

#line 1377 "IB01MD.f"
		    i__2 = nobr21;
#line 1377 "IB01MD.f"
		    for (i__ = 1; i__ <= i__2; ++i__) {
#line 1378 "IB01MD.f"
			dlacpy_("Full", &nobr21, l, &dwork[ishfty - lldrw + 1]
				, &ldrwrk, &dwork[ishfty], &ldrwrk, (ftnlen)4)
				;
#line 1381 "IB01MD.f"
			ishfty += lldrw;
#line 1382 "IB01MD.f"
/* L550: */
#line 1382 "IB01MD.f"
		    }

#line 1384 "IB01MD.f"
		}

/*              Workspace: need LDRWRK*2*(M+L)*NOBR + 4*(M+L)*NOBR. */

#line 1388 "IB01MD.f"
		mb04od_("Full", &nr, &c__0, &nsl, &r__[r_offset], ldr, &dwork[
			1], &ldrwrk, dum, &nr, dum, &nr, &dwork[itau], &dwork[
			jwork], (ftnlen)4);
#line 1391 "IB01MD.f"
		initi += nsf;
#line 1392 "IB01MD.f"
/* L560: */
#line 1392 "IB01MD.f"
	    }

#line 1394 "IB01MD.f"
	}

#line 1396 "IB01MD.f"
	if (! last) {
#line 1397 "IB01MD.f"
	    if (connec) {

/*              For sequential processing with connected data blocks, */
/*              save the remaining ("connection") elements of  U  and  Y */
/*              in the first  (M+L)*(2*NOBR-1)  locations of  DWORK. */

#line 1403 "IB01MD.f"
		if (*m > 0) {
#line 1403 "IB01MD.f"
		    dlacpy_("Full", &nobr21, m, &u[initi + 1 + u_dim1], ldu, &
			    dwork[1], &nobr21, (ftnlen)4);
#line 1403 "IB01MD.f"
		}
#line 1406 "IB01MD.f"
		dlacpy_("Full", &nobr21, l, &y[initi + 1 + y_dim1], ldy, &
			dwork[mmnobr - *m + 1], &nobr21, (ftnlen)4);
#line 1408 "IB01MD.f"
	    }

/*           Return to get new data. */

#line 1412 "IB01MD.f"
	    ++icycle;
#line 1413 "IB01MD.f"
	    if (icycle <= 100) {
#line 1413 "IB01MD.f"
		return 0;
#line 1413 "IB01MD.f"
	    }
#line 1415 "IB01MD.f"
	    *iwarn = 1;
#line 1416 "IB01MD.f"
	    icycle = 1;

#line 1418 "IB01MD.f"
	}

#line 1420 "IB01MD.f"
    }

/*     Return optimal workspace in  DWORK(1). */

#line 1424 "IB01MD.f"
    dwork[1] = (doublereal) maxwrk;
#line 1425 "IB01MD.f"
    if (last) {
#line 1426 "IB01MD.f"
	icycle = 1;
#line 1427 "IB01MD.f"
	maxwrk = 1;
#line 1428 "IB01MD.f"
	nsmpsm = 0;
#line 1429 "IB01MD.f"
    }
#line 1430 "IB01MD.f"
    return 0;

/* *** Last line of IB01MD *** */
} /* ib01md_ */

