#line 1 "SB16AD.f"
/* SB16AD.f -- translated by f2c (version 20100827).
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

#line 1 "SB16AD.f"
/* Subroutine */ int sb16ad_(char *dico, char *jobc, char *jobo, char *jobmr, 
	char *weight, char *equil, char *ordsel, integer *n, integer *m, 
	integer *p, integer *nc, integer *ncr, doublereal *alpha, doublereal *
	a, integer *lda, doublereal *b, integer *ldb, doublereal *c__, 
	integer *ldc, doublereal *d__, integer *ldd, doublereal *ac, integer *
	ldac, doublereal *bc, integer *ldbc, doublereal *cc, integer *ldcc, 
	doublereal *dc, integer *lddc, integer *ncs, doublereal *hsvc, 
	doublereal *tol1, doublereal *tol2, integer *iwork, doublereal *dwork,
	 integer *ldwork, integer *iwarn, integer *info, ftnlen dico_len, 
	ftnlen jobc_len, ftnlen jobo_len, ftnlen jobmr_len, ftnlen weight_len,
	 ftnlen equil_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, ac_dim1, ac_offset, b_dim1, b_offset, bc_dim1, 
	    bc_offset, c_dim1, c_offset, cc_dim1, cc_offset, d_dim1, d_offset,
	     dc_dim1, dc_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer ki, kr, mp, kt, ku, kw, lw;
    static logical bal, bta;
    static integer nnc, nra;
    static logical spa;
    static integer ncu, kti, nmr, ncu1;
    static logical perf;
    static integer ierr;
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    tb01kd_(char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen), ab09ix_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical istab, discr;
    extern /* Subroutine */ int sb16ay_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, integer *, doublereal *, integer *, integer *, ftnlen,
	     ftnlen, ftnlen, ftnlen);
    static logical ostab, leftw;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scalec, scaleo;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical fixord;
    static integer iwarnl;
    static doublereal alpwrk;
    static logical frwght, rightw;
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

/*     To compute a reduced order controller (Acr,Bcr,Ccr,Dcr) for an */
/*     original state-space controller representation (Ac,Bc,Cc,Dc) by */
/*     using the frequency-weighted square-root or balancing-free */
/*     square-root Balance & Truncate (B&T) or Singular Perturbation */
/*     Approximation (SPA) model reduction methods. The algorithm tries */
/*     to minimize the norm of the frequency-weighted error */

/*           ||V*(K-Kr)*W|| */

/*     where K and Kr are the transfer-function matrices of the original */
/*     and reduced order controllers, respectively. V and W are special */
/*     frequency-weighting transfer-function matrices constructed */
/*     to enforce closed-loop stability and/or closed-loop performance. */
/*     If G is the transfer-function matrix of the open-loop system, then */
/*     the following weightings V and W can be used: */
/*                      -1 */
/*      (a)   V = (I-G*K) *G, W = I - to enforce closed-loop stability; */
/*                              -1 */
/*      (b)   V = I,  W = (I-G*K) *G - to enforce closed-loop stability; */
/*                      -1              -1 */
/*      (c)   V = (I-G*K) *G, W = (I-G*K)  - to enforce closed-loop */
/*            stability and performance. */

/*     G has the state space representation (A,B,C,D). */
/*     If K is unstable, only the ALPHA-stable part of K is reduced. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original controller as follows: */
/*             = 'C':  continuous-time controller; */
/*             = 'D':  discrete-time controller. */

/*     JOBC    CHARACTER*1 */
/*             Specifies the choice of frequency-weighted controllability */
/*             Grammian as follows: */
/*             = 'S': choice corresponding to standard Enns' method [1]; */
/*             = 'E': choice corresponding to the stability enhanced */
/*                    modified Enns' method of [2]. */

/*     JOBO    CHARACTER*1 */
/*             Specifies the choice of frequency-weighted observability */
/*             Grammian as follows: */
/*             = 'S': choice corresponding to standard Enns' method [1]; */
/*             = 'E': choice corresponding to the stability enhanced */
/*                    modified combination method of [2]. */

/*     JOBMR   CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root B&T method; */
/*             = 'F':  use the balancing-free square-root B&T method; */
/*             = 'S':  use the square-root SPA method; */
/*             = 'P':  use the balancing-free square-root SPA method. */

/*     WEIGHT  CHARACTER*1 */
/*             Specifies the type of frequency-weighting, as follows: */
/*             = 'N':  no weightings are used (V = I, W = I); */
/*             = 'O':  stability enforcing left (output) weighting */
/*                               -1 */
/*                     V = (I-G*K) *G is used (W = I); */
/*             = 'I':  stability enforcing right (input) weighting */
/*                               -1 */
/*                     W = (I-G*K) *G is used (V = I); */
/*             = 'P':  stability and performance enforcing weightings */
/*                               -1                -1 */
/*                     V = (I-G*K) *G ,  W = (I-G*K)  are used. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the triplets (A,B,C) and (Ac,Bc,Cc) as */
/*             follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     ORDSEL  CHARACTER*1 */
/*             Specifies the order selection method as follows: */
/*             = 'F':  the resulting order NCR is fixed; */
/*             = 'A':  the resulting order NCR is automatically */
/*                     determined on basis of the given tolerance TOL1. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the open-loop system state-space */
/*             representation, i.e., the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NC      (input) INTEGER */
/*             The order of the controller state-space representation, */
/*             i.e., the order of the matrix AC.  NC >= 0. */

/*     NCR     (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NCR is the desired order of */
/*             the resulting reduced order controller.  0 <= NCR <= NC. */
/*             On exit, if INFO = 0, NCR is the order of the resulting */
/*             reduced order controller. For a controller with NCU */
/*             ALPHA-unstable eigenvalues and NCS ALPHA-stable */
/*             eigenvalues (NCU+NCS = NC), NCR is set as follows: */
/*             if ORDSEL = 'F', NCR is equal to */
/*             NCU+MIN(MAX(0,NCR-NCU),NCMIN), where NCR is the desired */
/*             order on entry, NCMIN is the number of frequency-weighted */
/*             Hankel singular values greater than NCS*EPS*S1, EPS is the */
/*             machine precision (see LAPACK Library Routine DLAMCH) and */
/*             S1 is the largest Hankel singular value (computed in */
/*             HSVC(1)); NCR can be further reduced to ensure */
/*             HSVC(NCR-NCU) > HSVC(NCR+1-NCU); */
/*             if ORDSEL = 'A', NCR is the sum of NCU and the number of */
/*             Hankel singular values greater than MAX(TOL1,NCS*EPS*S1). */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             Specifies the ALPHA-stability boundary for the eigenvalues */
/*             of the state dynamics matrix AC. For a continuous-time */
/*             controller (DICO = 'C'), ALPHA <= 0 is the boundary value */
/*             for the real parts of eigenvalues; for a discrete-time */
/*             controller (DICO = 'D'), 0 <= ALPHA <= 1 represents the */
/*             boundary value for the moduli of eigenvalues. */
/*             The ALPHA-stability domain does not include the boundary. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A of the open-loop */
/*             system. */
/*             On exit, if INFO = 0 and EQUIL = 'S', the leading N-by-N */
/*             part of this array contains the scaled state dynamics */
/*             matrix of the open-loop system. */
/*             If EQUIL = 'N', this array is unchanged on exit. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the input/state matrix B of the open-loop system. */
/*             On exit, if INFO = 0 and EQUIL = 'S', the leading N-by-M */
/*             part of this array contains the scaled input/state matrix */
/*             of the open-loop system. */
/*             If EQUIL = 'N', this array is unchanged on exit. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the state/output matrix C of the open-loop system. */
/*             On exit, if INFO = 0 and EQUIL = 'S', the leading P-by-N */
/*             part of this array contains the scaled state/output matrix */
/*             of the open-loop system. */
/*             If EQUIL = 'N', this array is unchanged on exit. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             The leading P-by-M part of this array must contain the */
/*             input/output matrix D of the open-loop system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     AC      (input/output) DOUBLE PRECISION array, dimension (LDAC,NC) */
/*             On entry, the leading NC-by-NC part of this array must */
/*             contain the state dynamics matrix Ac of the original */
/*             controller. */
/*             On exit, if INFO = 0, the leading NCR-by-NCR part of this */
/*             array contains the state dynamics matrix Acr of the */
/*             reduced controller. The resulting Ac has a */
/*             block-diagonal form with two blocks. */
/*             For a system with NCU ALPHA-unstable eigenvalues and */
/*             NCS ALPHA-stable eigenvalues (NCU+NCS = NC), the leading */
/*             NCU-by-NCU block contains the unreduced part of Ac */
/*             corresponding to the ALPHA-unstable eigenvalues. */
/*             The trailing (NCR+NCS-NC)-by-(NCR+NCS-NC) block contains */
/*             the reduced part of Ac corresponding to ALPHA-stable */
/*             eigenvalues. */

/*     LDAC    INTEGER */
/*             The leading dimension of array AC.  LDAC >= MAX(1,NC). */

/*     BC      (input/output) DOUBLE PRECISION array, dimension (LDBC,P) */
/*             On entry, the leading NC-by-P part of this array must */
/*             contain the input/state matrix Bc of the original */
/*             controller. */
/*             On exit, if INFO = 0, the leading NCR-by-P part of this */
/*             array contains the input/state matrix Bcr of the reduced */
/*             controller. */

/*     LDBC    INTEGER */
/*             The leading dimension of array BC.  LDBC >= MAX(1,NC). */

/*     CC      (input/output) DOUBLE PRECISION array, dimension (LDCC,NC) */
/*             On entry, the leading M-by-NC part of this array must */
/*             contain the state/output matrix Cc of the original */
/*             controller. */
/*             On exit, if INFO = 0, the leading M-by-NCR part of this */
/*             array contains the state/output matrix Ccr of the reduced */
/*             controller. */

/*     LDCC    INTEGER */
/*             The leading dimension of array CC.  LDCC >= MAX(1,M). */

/*     DC      (input/output) DOUBLE PRECISION array, dimension (LDDC,P) */
/*             On entry, the leading M-by-P part of this array must */
/*             contain the input/output matrix Dc of the original */
/*             controller. */
/*             On exit, if INFO = 0, the leading M-by-P part of this */
/*             array contains the input/output matrix Dcr of the reduced */
/*             controller. */

/*     LDDC    INTEGER */
/*             The leading dimension of array DC.  LDDC >= MAX(1,M). */

/*     NCS     (output) INTEGER */
/*             The dimension of the ALPHA-stable part of the controller. */

/*     HSVC    (output) DOUBLE PRECISION array, dimension (NC) */
/*             If INFO = 0, the leading NCS elements of this array */
/*             contain the frequency-weighted Hankel singular values, */
/*             ordered decreasingly, of the ALPHA-stable part of the */
/*             controller. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of the reduced controller. */
/*             For model reduction, the recommended value is */
/*             TOL1 = c*S1, where c is a constant in the */
/*             interval [0.00001,0.001], and S1 is the largest */
/*             frequency-weighted Hankel singular value of the */
/*             ALPHA-stable part of the original controller */
/*             (computed in HSVC(1)). */
/*             If TOL1 <= 0 on entry, the used default value is */
/*             TOL1 = NCS*EPS*S1, where NCS is the number of */
/*             ALPHA-stable eigenvalues of Ac and EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the ALPHA-stable part of the given */
/*             controller. The recommended value is TOL2 = NCS*EPS*S1. */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension MAX(1,LIWRK1,LIWRK2) */
/*             LIWRK1 = 0,       if JOBMR  = 'B'; */
/*             LIWRK1 = NC,      if JOBMR  = 'F'; */
/*             LIWRK1 = 2*NC,    if JOBMR  = 'S' or 'P'; */
/*             LIWRK2 = 0,       if WEIGHT = 'N'; */
/*             LIWRK2 = 2*(M+P), if WEIGHT = 'O', 'I', or 'P'. */
/*             On exit, if INFO = 0, IWORK(1) contains NCMIN, the order */
/*             of the computed minimal realization of the stable part of */
/*             the controller. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 2*NC*NC + MAX( 1, LFREQ, LSQRED ), */
/*             where */
/*             LFREQ = (N+NC)*(N+NC+2*M+2*P)+ */
/*                     MAX((N+NC)*(N+NC+MAX(N+NC,M,P)+7), (M+P)*(M+P+4)) */
/*                                      if WEIGHT = 'I' or 'O' or 'P'; */
/*             LFREQ  = NC*(MAX(M,P)+5) if WEIGHT = 'N' and EQUIL = 'N'; */
/*             LFREQ  = MAX(N,NC*(MAX(M,P)+5)) if WEIGHT = 'N' and */
/*                                                EQUIL  = 'S'; */
/*             LSQRED = MAX( 1, 2*NC*NC+5*NC ); */
/*             For optimum performance LDWORK should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NCR is greater */
/*                   than NSMIN, the sum of the order of the */
/*                   ALPHA-unstable part and the order of a minimal */
/*                   realization of the ALPHA-stable part of the given */
/*                   controller; in this case, the resulting NCR is set */
/*                   equal to NSMIN; */
/*             = 2:  with ORDSEL = 'F', the selected order NCR */
/*                   corresponds to repeated singular values for the */
/*                   ALPHA-stable part of the controller, which are */
/*                   neither all included nor all excluded from the */
/*                   reduced model; in this case, the resulting NCR is */
/*                   automatically decreased to exclude all repeated */
/*                   singular values; */
/*             = 3:  with ORDSEL = 'F', the selected order NCR is less */
/*                   than the order of the ALPHA-unstable part of the */
/*                   given controller. In this case NCR is set equal to */
/*                   the order of the ALPHA-unstable part. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the closed-loop system is not well-posed; */
/*                   its feedthrough matrix is (numerically) singular; */
/*             = 2:  the computation of the real Schur form of the */
/*                   closed-loop state matrix failed; */
/*             = 3:  the closed-loop state matrix is not stable; */
/*             = 4:  the solution of a symmetric eigenproblem failed; */
/*             = 5:  the computation of the ordered real Schur form of Ac */
/*                   failed; */
/*             = 6:  the separation of the ALPHA-stable/unstable */
/*                   diagonal blocks failed because of very close */
/*                   eigenvalues; */
/*             = 7:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let K be the transfer-function matrix of the original linear */
/*     controller */

/*          d[xc(t)] = Ac*xc(t) + Bc*y(t) */
/*          u(t)     = Cc*xc(t) + Dc*y(t),                      (1) */

/*     where d[xc(t)] is dxc(t)/dt for a continuous-time system and */
/*     xc(t+1) for a discrete-time system. The subroutine SB16AD */
/*     determines the matrices of a reduced order controller */

/*          d[z(t)] = Acr*z(t) + Bcr*y(t) */
/*          u(t)    = Ccr*z(t) + Dcr*y(t),                      (2) */

/*     such that the corresponding transfer-function matrix Kr minimizes */
/*     the norm of the frequency-weighted error */

/*             V*(K-Kr)*W,                                      (3) */

/*     where V and W are special stable transfer-function matrices */
/*     chosen to enforce stability and/or performance of the closed-loop */
/*     system [3] (see description of the parameter WEIGHT). */

/*     The following procedure is used to reduce K in conjunction */
/*     with the frequency-weighted balancing approach of [2] */
/*     (see also [3]): */

/*     1) Decompose additively K, of order NC, as */

/*          K = K1 + K2, */

/*        such that K1 has only ALPHA-stable poles and K2, of order NCU, */
/*        has only ALPHA-unstable poles. */

/*     2) Compute for K1 a B&T or SPA frequency-weighted approximation */
/*        K1r of order NCR-NCU using the frequency-weighted balancing */
/*        approach of [1] in conjunction with accuracy enhancing */
/*        techniques specified by the parameter JOBMR. */

/*     3) Assemble the reduced model Kr as */

/*           Kr = K1r + K2. */

/*     For the reduction of the ALPHA-stable part, several accuracy */
/*     enhancing techniques can be employed (see [2] for details). */

/*     If JOBMR = 'B', the square-root B&T method of [1] is used. */

/*     If JOBMR = 'F', the balancing-free square-root version of the */
/*     B&T method [1] is used. */

/*     If JOBMR = 'S', the square-root version of the SPA method [2,3] */
/*     is used. */

/*     If JOBMR = 'P', the balancing-free square-root version of the */
/*     SPA method [2,3] is used. */

/*     For each of these methods, two left and right truncation matrices */
/*     are determined using the Cholesky factors of an input */
/*     frequency-weighted controllability Grammian P and an output */
/*     frequency-weighted observability Grammian Q. */
/*     P and Q are determined as the leading NC-by-NC diagonal blocks */
/*     of the controllability Grammian of K*W and of the */
/*     observability Grammian of V*K. Special techniques developed in [2] */
/*     are used to compute the Cholesky factors of P and Q directly */
/*     (see also SLICOT Library routine SB16AY). */
/*     The frequency-weighted Hankel singular values HSVC(1), ...., */
/*     HSVC(NC) are computed as the square roots of the eigenvalues */
/*     of the product P*Q. */

/*     REFERENCES */

/*     [1] Enns, D. */
/*         Model reduction with balanced realizations: An error bound */
/*         and a frequency weighted generalization. */
/*         Proc. 23-th CDC, Las Vegas, pp. 127-132, 1984. */

/*     [2] Varga, A. and Anderson, B.D.O. */
/*         Square-root balancing-free methods for frequency-weighted */
/*         balancing related model reduction. */
/*         (report in preparation) */

/*     [3] Anderson, B.D.O and Liu, Y. */
/*         Controller reduction: concepts and approaches. */
/*         IEEE Trans. Autom. Control, Vol. 34, pp. 802-812, 1989. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root */
/*     techniques. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, Sept. 2000. */
/*     D. Sima, University of Bucharest, Sept. 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Sept.2000. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000, */
/*              Sep. 2001. */

/*     KEYWORDS */

/*     Controller reduction, frequency weighting, multivariable system, */
/*     state-space model, state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 489 "SB16AD.f"
    /* Parameter adjustments */
#line 489 "SB16AD.f"
    a_dim1 = *lda;
#line 489 "SB16AD.f"
    a_offset = 1 + a_dim1;
#line 489 "SB16AD.f"
    a -= a_offset;
#line 489 "SB16AD.f"
    b_dim1 = *ldb;
#line 489 "SB16AD.f"
    b_offset = 1 + b_dim1;
#line 489 "SB16AD.f"
    b -= b_offset;
#line 489 "SB16AD.f"
    c_dim1 = *ldc;
#line 489 "SB16AD.f"
    c_offset = 1 + c_dim1;
#line 489 "SB16AD.f"
    c__ -= c_offset;
#line 489 "SB16AD.f"
    d_dim1 = *ldd;
#line 489 "SB16AD.f"
    d_offset = 1 + d_dim1;
#line 489 "SB16AD.f"
    d__ -= d_offset;
#line 489 "SB16AD.f"
    ac_dim1 = *ldac;
#line 489 "SB16AD.f"
    ac_offset = 1 + ac_dim1;
#line 489 "SB16AD.f"
    ac -= ac_offset;
#line 489 "SB16AD.f"
    bc_dim1 = *ldbc;
#line 489 "SB16AD.f"
    bc_offset = 1 + bc_dim1;
#line 489 "SB16AD.f"
    bc -= bc_offset;
#line 489 "SB16AD.f"
    cc_dim1 = *ldcc;
#line 489 "SB16AD.f"
    cc_offset = 1 + cc_dim1;
#line 489 "SB16AD.f"
    cc -= cc_offset;
#line 489 "SB16AD.f"
    dc_dim1 = *lddc;
#line 489 "SB16AD.f"
    dc_offset = 1 + dc_dim1;
#line 489 "SB16AD.f"
    dc -= dc_offset;
#line 489 "SB16AD.f"
    --hsvc;
#line 489 "SB16AD.f"
    --iwork;
#line 489 "SB16AD.f"
    --dwork;
#line 489 "SB16AD.f"

#line 489 "SB16AD.f"
    /* Function Body */
#line 489 "SB16AD.f"
    *info = 0;
#line 490 "SB16AD.f"
    *iwarn = 0;
#line 491 "SB16AD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 492 "SB16AD.f"
    bta = lsame_(jobmr, "B", (ftnlen)1, (ftnlen)1) || lsame_(jobmr, "F", (
	    ftnlen)1, (ftnlen)1);
#line 493 "SB16AD.f"
    spa = lsame_(jobmr, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobmr, "P", (
	    ftnlen)1, (ftnlen)1);
#line 494 "SB16AD.f"
    bal = lsame_(jobmr, "B", (ftnlen)1, (ftnlen)1) || lsame_(jobmr, "S", (
	    ftnlen)1, (ftnlen)1);
#line 495 "SB16AD.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
#line 496 "SB16AD.f"
    istab = lsame_(weight, "I", (ftnlen)1, (ftnlen)1);
#line 497 "SB16AD.f"
    ostab = lsame_(weight, "O", (ftnlen)1, (ftnlen)1);
#line 498 "SB16AD.f"
    perf = lsame_(weight, "P", (ftnlen)1, (ftnlen)1);
#line 499 "SB16AD.f"
    leftw = ostab || perf;
#line 500 "SB16AD.f"
    rightw = istab || perf;
#line 501 "SB16AD.f"
    frwght = leftw || rightw;

#line 503 "SB16AD.f"
    lw = 1;
#line 504 "SB16AD.f"
    nnc = *n + *nc;
#line 505 "SB16AD.f"
    mp = *m + *p;
#line 506 "SB16AD.f"
    if (frwght) {
/* Computing MAX */
/* Computing MAX */
#line 507 "SB16AD.f"
	i__3 = max(nnc,*m);
#line 507 "SB16AD.f"
	i__1 = nnc * (nnc + max(i__3,*p) + 7), i__2 = mp * (mp + 4);
#line 507 "SB16AD.f"
	lw = nnc * (nnc + (mp << 1)) + max(i__1,i__2);
#line 509 "SB16AD.f"
    } else {
#line 510 "SB16AD.f"
	lw = *nc * (max(*m,*p) + 5);
#line 511 "SB16AD.f"
	if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {
#line 511 "SB16AD.f"
	    lw = max(*n,lw);
#line 511 "SB16AD.f"
	}
#line 513 "SB16AD.f"
    }
/* Computing MAX */
#line 514 "SB16AD.f"
    i__1 = max(1,lw), i__2 = *nc * ((*nc << 1) + 5);
#line 514 "SB16AD.f"
    lw = (*nc << 1) * *nc + max(i__1,i__2);

/*     Check the input scalar arguments. */

#line 518 "SB16AD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 519 "SB16AD.f"
	*info = -1;
#line 520 "SB16AD.f"
    } else if (! (lsame_(jobc, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobc, 
	    "E", (ftnlen)1, (ftnlen)1))) {
#line 522 "SB16AD.f"
	*info = -2;
#line 523 "SB16AD.f"
    } else if (! (lsame_(jobo, "S", (ftnlen)1, (ftnlen)1) || lsame_(jobo, 
	    "E", (ftnlen)1, (ftnlen)1))) {
#line 525 "SB16AD.f"
	*info = -3;
#line 526 "SB16AD.f"
    } else if (! (bta || spa)) {
#line 527 "SB16AD.f"
	*info = -4;
#line 528 "SB16AD.f"
    } else if (! (frwght || lsame_(weight, "N", (ftnlen)1, (ftnlen)1))) {
#line 529 "SB16AD.f"
	*info = -5;
#line 530 "SB16AD.f"
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 532 "SB16AD.f"
	*info = -6;
#line 533 "SB16AD.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 534 "SB16AD.f"
	*info = -7;
#line 535 "SB16AD.f"
    } else if (*n < 0) {
#line 536 "SB16AD.f"
	*info = -8;
#line 537 "SB16AD.f"
    } else if (*m < 0) {
#line 538 "SB16AD.f"
	*info = -9;
#line 539 "SB16AD.f"
    } else if (*p < 0) {
#line 540 "SB16AD.f"
	*info = -10;
#line 541 "SB16AD.f"
    } else if (*nc < 0) {
#line 542 "SB16AD.f"
	*info = -11;
#line 543 "SB16AD.f"
    } else if (fixord && (*ncr < 0 || *ncr > *nc)) {
#line 544 "SB16AD.f"
	*info = -12;
#line 545 "SB16AD.f"
    } else if (discr && (*alpha < 0. || *alpha > 1.) || ! discr && *alpha > 
	    0.) {
#line 547 "SB16AD.f"
	*info = -13;
#line 548 "SB16AD.f"
    } else if (*lda < max(1,*n)) {
#line 549 "SB16AD.f"
	*info = -15;
#line 550 "SB16AD.f"
    } else if (*ldb < max(1,*n)) {
#line 551 "SB16AD.f"
	*info = -17;
#line 552 "SB16AD.f"
    } else if (*ldc < max(1,*p)) {
#line 553 "SB16AD.f"
	*info = -19;
#line 554 "SB16AD.f"
    } else if (*ldd < max(1,*p)) {
#line 555 "SB16AD.f"
	*info = -21;
#line 556 "SB16AD.f"
    } else if (*ldac < max(1,*nc)) {
#line 557 "SB16AD.f"
	*info = -23;
#line 558 "SB16AD.f"
    } else if (*ldbc < max(1,*nc)) {
#line 559 "SB16AD.f"
	*info = -25;
#line 560 "SB16AD.f"
    } else if (*ldcc < max(1,*m)) {
#line 561 "SB16AD.f"
	*info = -27;
#line 562 "SB16AD.f"
    } else if (*lddc < max(1,*m)) {
#line 563 "SB16AD.f"
	*info = -29;
#line 564 "SB16AD.f"
    } else if (*tol2 > 0. && ! fixord && *tol2 > *tol1) {
#line 565 "SB16AD.f"
	*info = -33;
#line 566 "SB16AD.f"
    } else if (*ldwork < lw) {
#line 567 "SB16AD.f"
	*info = -36;
#line 568 "SB16AD.f"
    }

#line 570 "SB16AD.f"
    if (*info != 0) {

/*        Error return. */

#line 574 "SB16AD.f"
	i__1 = -(*info);
#line 574 "SB16AD.f"
	xerbla_("SB16AD", &i__1, (ftnlen)6);
#line 575 "SB16AD.f"
	return 0;
#line 576 "SB16AD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 580 "SB16AD.f"
    i__1 = min(*nc,*m);
#line 580 "SB16AD.f"
    if (min(i__1,*p) == 0) {
#line 581 "SB16AD.f"
	*ncr = 0;
#line 582 "SB16AD.f"
	*ncs = 0;
#line 583 "SB16AD.f"
	iwork[1] = 0;
#line 584 "SB16AD.f"
	dwork[1] = 1.;
#line 585 "SB16AD.f"
	return 0;
#line 586 "SB16AD.f"
    }

#line 588 "SB16AD.f"
    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C and AC, BC and CC; */
/*        A <- inv(T1)*A*T1, B <- inv(T1)*B and C <- C*T1, where T1 is a */
/*        diagonal matrix; */
/*        AC <- inv(T2)*AC*T2, BC <- inv(T2)*BC and CC <- CC*T2, where T2 */
/*        is a diagonal matrix. */

/*        Real workspace: need MAX(N,NC). */

#line 598 "SB16AD.f"
	maxred = 100.;
#line 599 "SB16AD.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 601 "SB16AD.f"
	maxred = 100.;
#line 602 "SB16AD.f"
	tb01id_("All", nc, p, m, &maxred, &ac[ac_offset], ldac, &bc[bc_offset]
		, ldbc, &cc[cc_offset], ldcc, &dwork[1], info, (ftnlen)3);
#line 604 "SB16AD.f"
    }

/*     Correct the value of ALPHA to ensure stability. */

#line 608 "SB16AD.f"
    alpwrk = *alpha;
#line 609 "SB16AD.f"
    if (discr) {
#line 610 "SB16AD.f"
	if (*alpha == 1.) {
#line 610 "SB16AD.f"
	    alpwrk = 1. - sqrt(dlamch_("E", (ftnlen)1));
#line 610 "SB16AD.f"
	}
#line 611 "SB16AD.f"
    } else {
#line 612 "SB16AD.f"
	if (*alpha == 0.) {
#line 612 "SB16AD.f"
	    alpwrk = -sqrt(dlamch_("E", (ftnlen)1));
#line 612 "SB16AD.f"
	}
#line 613 "SB16AD.f"
    }

/*     Reduce Ac to a block-diagonal real Schur form, with the */
/*     ALPHA-unstable part in the leading diagonal position, using a */
/*     non-orthogonal similarity transformation, AC <- inv(T)*AC*T, and */
/*     apply the transformation to BC and CC: */
/*     BC <- inv(T)*BC and CC <- CC*T. */

/*     Workspace:  need   NC*(NC+5); */
/*                 prefer larger. */

#line 624 "SB16AD.f"
    wrkopt = 1;
#line 625 "SB16AD.f"
    ku = 1;
#line 626 "SB16AD.f"
    kr = ku + *nc * *nc;
#line 627 "SB16AD.f"
    ki = kr + *nc;
#line 628 "SB16AD.f"
    kw = ki + *nc;

#line 630 "SB16AD.f"
    i__1 = *ldwork - kw + 1;
#line 630 "SB16AD.f"
    tb01kd_(dico, "Unstable", "General", nc, p, m, &alpwrk, &ac[ac_offset], 
	    ldac, &bc[bc_offset], ldbc, &cc[cc_offset], ldcc, &ncu, &dwork[ku]
	    , nc, &dwork[kr], &dwork[ki], &dwork[kw], &i__1, &ierr, (ftnlen)1,
	     (ftnlen)8, (ftnlen)7);

#line 634 "SB16AD.f"
    if (ierr != 0) {
#line 635 "SB16AD.f"
	if (ierr != 3) {
#line 636 "SB16AD.f"
	    *info = 5;
#line 637 "SB16AD.f"
	} else {
#line 638 "SB16AD.f"
	    *info = 6;
#line 639 "SB16AD.f"
	}
#line 640 "SB16AD.f"
	return 0;
#line 641 "SB16AD.f"
    }
/* Computing MAX */
#line 642 "SB16AD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 642 "SB16AD.f"
    wrkopt = max(i__1,i__2);

#line 644 "SB16AD.f"
    iwarnl = 0;
#line 645 "SB16AD.f"
    *ncs = *nc - ncu;
#line 646 "SB16AD.f"
    if (fixord) {
/* Computing MAX */
#line 647 "SB16AD.f"
	i__1 = 0, i__2 = *ncr - ncu;
#line 647 "SB16AD.f"
	nra = max(i__1,i__2);
#line 648 "SB16AD.f"
	if (*ncr < ncu) {
#line 648 "SB16AD.f"
	    iwarnl = 3;
#line 648 "SB16AD.f"
	}
#line 650 "SB16AD.f"
    } else {
#line 651 "SB16AD.f"
	nra = 0;
#line 652 "SB16AD.f"
    }

/*     Finish if only unstable part is present. */

#line 656 "SB16AD.f"
    if (*ncs == 0) {
#line 657 "SB16AD.f"
	*ncr = ncu;
#line 658 "SB16AD.f"
	iwork[1] = 0;
#line 659 "SB16AD.f"
	dwork[1] = (doublereal) wrkopt;
#line 660 "SB16AD.f"
	return 0;
#line 661 "SB16AD.f"
    }

/*     Allocate working storage. */

#line 665 "SB16AD.f"
    kt = 1;
#line 666 "SB16AD.f"
    kti = kt + *nc * *nc;
#line 667 "SB16AD.f"
    kw = kti + *nc * *nc;

/*     Compute in DWORK(KTI) and DWORK(KT) the Cholesky factors S and R */
/*     of the frequency-weighted controllability and observability */
/*     Grammians, respectively. */

/*     Real workspace:  need  2*NC*NC + MAX( 1, LFREQ ), */
/*                      where */
/*                      LFREQ = (N+NC)*(N+NC+2*M+2*P)+ */
/*                              MAX((N+NC)*(N+NC+MAX(N+NC,M,P)+7), */
/*                                  (M+P)*(M+P+4)) */
/*                                         if WEIGHT = 'I' or 'O' or 'P'; */
/*                      LFREQ = NCS*(MAX(M,P)+5) if WEIGHT = 'N'; */
/*                      prefer larger. */
/*     Integer workspace:      2*(M+P) if WEIGHT = 'I' or 'O' or 'P'; */
/*                             0,      if WEIGHT = 'N'. */

#line 684 "SB16AD.f"
    i__1 = *ldwork - kw + 1;
#line 684 "SB16AD.f"
    sb16ay_(dico, jobc, jobo, weight, n, m, p, nc, ncs, &a[a_offset], lda, &b[
	    b_offset], ldb, &c__[c_offset], ldc, &d__[d_offset], ldd, &ac[
	    ac_offset], ldac, &bc[bc_offset], ldbc, &cc[cc_offset], ldcc, &dc[
	    dc_offset], lddc, &scalec, &scaleo, &dwork[kti], nc, &dwork[kt], 
	    nc, &iwork[1], &dwork[kw], &i__1, info, (ftnlen)1, (ftnlen)1, (
	    ftnlen)1, (ftnlen)1);
#line 689 "SB16AD.f"
    if (*info != 0) {
#line 689 "SB16AD.f"
	return 0;
#line 689 "SB16AD.f"
    }
/* Computing MAX */
#line 691 "SB16AD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 691 "SB16AD.f"
    wrkopt = max(i__1,i__2);

/*     Compute a BTA or SPA of the stable part. */
/*     Real workspace:  need   2*NC*NC + MAX( 1, 2*NC*NC+5*NC, */
/*                                               NC*MAX(M,P) ); */
/*                      prefer larger. */
/*     Integer workspace:      0,     if JOBMR = 'B'; */
/*                             NC,    if JOBMR = 'F'; */
/*                             2*NC,  if JOBMR = 'S' or 'P'. */

#line 701 "SB16AD.f"
    ncu1 = ncu + 1;
#line 702 "SB16AD.f"
    i__1 = *ldwork - kw + 1;
#line 702 "SB16AD.f"
    ab09ix_(dico, jobmr, "Schur", ordsel, ncs, p, m, &nra, &scalec, &scaleo, &
	    ac[ncu1 + ncu1 * ac_dim1], ldac, &bc[ncu1 + bc_dim1], ldbc, &cc[
	    ncu1 * cc_dim1 + 1], ldcc, &dc[dc_offset], lddc, &dwork[kti], nc, 
	    &dwork[kt], nc, &nmr, &hsvc[1], tol1, tol2, &iwork[1], &dwork[kw],
	     &i__1, iwarn, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen)5, (ftnlen)1);
#line 707 "SB16AD.f"
    *iwarn = max(*iwarn,iwarnl);
#line 708 "SB16AD.f"
    if (ierr != 0) {
#line 709 "SB16AD.f"
	*info = 7;
#line 710 "SB16AD.f"
	return 0;
#line 711 "SB16AD.f"
    }
#line 712 "SB16AD.f"
    *ncr = nra + ncu;
#line 713 "SB16AD.f"
    iwork[1] = nmr;

/* Computing MAX */
#line 715 "SB16AD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 715 "SB16AD.f"
    dwork[1] = (doublereal) max(i__1,i__2);

#line 717 "SB16AD.f"
    return 0;
/* *** Last line of SB16AD *** */
} /* sb16ad_ */

