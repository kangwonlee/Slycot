#line 1 "AB09HD.f"
/* AB09HD.f -- translated by f2c (version 20100827).
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

#line 1 "AB09HD.f"
/* Table of constant values */

static doublereal c_b19 = 0.;
static doublereal c_b24 = 1.;
static doublereal c_b29 = .66666666666666663;

/* Subroutine */ int ab09hd_(char *dico, char *job, char *equil, char *ordsel,
	 integer *n, integer *m, integer *p, integer *nr, doublereal *alpha, 
	doublereal *beta, doublereal *a, integer *lda, doublereal *b, integer 
	*ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	integer *ns, doublereal *hsv, doublereal *tol1, doublereal *tol2, 
	integer *iwork, doublereal *dwork, integer *ldwork, logical *bwork, 
	integer *iwarn, integer *info, ftnlen dico_len, ftnlen job_len, 
	ftnlen equil_len, ftnlen ordsel_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer n2, kb, kd, mb, nn, kt, ku, kw, lw, nu, nu1;
    static logical bta;
    static integer nra;
    static logical spa;
    static integer kti, kwi, nmr, kwr, lwr, ierr;
    static doublereal epsm;
    extern /* Subroutine */ int ab04md_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen),
	     tb01id_(char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen), tb01kd_(char *, char 
	    *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen), ab09hy_(integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, logical *, integer *), ab09ix_(char *, char *, char *, 
	    char *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal scalec, scaleo;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen);
    static doublereal ricond;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred;
    static logical fixord, lequil;
    static integer iwarnl, wrkopt;


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

/*     To compute a reduced order model (Ar,Br,Cr,Dr) for an original */
/*     state-space representation (A,B,C,D) by using the stochastic */
/*     balancing approach in conjunction with the square-root or */
/*     the balancing-free square-root Balance & Truncate (B&T) */
/*     or Singular Perturbation Approximation (SPA) model reduction */
/*     methods for the ALPHA-stable part of the system. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the original system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOB     CHARACTER*1 */
/*             Specifies the model reduction approach to be used */
/*             as follows: */
/*             = 'B':  use the square-root Balance & Truncate method; */
/*             = 'F':  use the balancing-free square-root */
/*                     Balance & Truncate method; */
/*             = 'S':  use the square-root Singular Perturbation */
/*                     Approximation method; */
/*             = 'P':  use the balancing-free square-root */
/*                     Singular Perturbation Approximation method. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the triplet (A,B,C) as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     ORDSEL  CHARACTER*1 */
/*             Specifies the order selection method as follows: */
/*             = 'F':  the resulting order NR is fixed; */
/*             = 'A':  the resulting order NR is automatically determined */
/*                     on basis of the given tolerance TOL1. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the original state-space representation, */
/*             i.e., the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */
/*             P <= M if BETA = 0. */

/*     NR      (input/output) INTEGER */
/*             On entry with ORDSEL = 'F', NR is the desired order of the */
/*             resulting reduced order system.  0 <= NR <= N. */
/*             On exit, if INFO = 0, NR is the order of the resulting */
/*             reduced order model. For a system with NU ALPHA-unstable */
/*             eigenvalues and NS ALPHA-stable eigenvalues (NU+NS = N), */
/*             NR is set as follows: if ORDSEL = 'F', NR is equal to */
/*             NU+MIN(MAX(0,NR-NU),NMIN), where NR is the desired order */
/*             on entry, and NMIN is the order of a minimal realization */
/*             of the ALPHA-stable part of the given system; NMIN is */
/*             determined as the number of Hankel singular values greater */
/*             than NS*EPS, where EPS is the machine precision */
/*             (see LAPACK Library Routine DLAMCH); */
/*             if ORDSEL = 'A', NR is the sum of NU and the number of */
/*             Hankel singular values greater than MAX(TOL1,NS*EPS); */
/*             NR can be further reduced to ensure that */
/*             HSV(NR-NU) > HSV(NR+1-NU). */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             Specifies the ALPHA-stability boundary for the eigenvalues */
/*             of the state dynamics matrix A. For a continuous-time */
/*             system (DICO = 'C'), ALPHA <= 0 is the boundary value for */
/*             the real parts of eigenvalues, while for a discrete-time */
/*             system (DICO = 'D'), 0 <= ALPHA <= 1 represents the */
/*             boundary value for the moduli of eigenvalues. */
/*             The ALPHA-stability domain does not include the boundary. */

/*     BETA    (input) DOUBLE PRECISION */
/*             BETA > 0 specifies the absolute/relative error weighting */
/*             parameter. A large positive value of BETA favours the */
/*             minimization of the absolute approximation error, while a */
/*             small value of BETA is appropriate for the minimization */
/*             of the relative error. */
/*             BETA = 0 means a pure relative error method and can be */
/*             used only if rank(D) = P. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, if INFO = 0, the leading NR-by-NR part of this */
/*             array contains the state dynamics matrix Ar of the reduced */
/*             order system. */
/*             The resulting A has a block-diagonal form with two blocks. */
/*             For a system with NU ALPHA-unstable eigenvalues and */
/*             NS ALPHA-stable eigenvalues (NU+NS = N), the leading */
/*             NU-by-NU block contains the unreduced part of A */
/*             corresponding to ALPHA-unstable eigenvalues in an */
/*             upper real Schur form. */
/*             The trailing (NR+NS-N)-by-(NR+NS-N) block contains */
/*             the reduced part of A corresponding to ALPHA-stable */
/*             eigenvalues. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */
/*             On exit, if INFO = 0, the leading NR-by-M part of this */
/*             array contains the input/state matrix Br of the reduced */
/*             order system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */
/*             On exit, if INFO = 0, the leading P-by-NR part of this */
/*             array contains the state/output matrix Cr of the reduced */
/*             order system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     D       (input/output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             On entry, the leading P-by-M part of this array must */
/*             contain the original input/output matrix D. */
/*             On exit, if INFO = 0, the leading P-by-M part of this */
/*             array contains the input/output matrix Dr of the reduced */
/*             order system. */

/*     LDD     INTEGER */
/*             The leading dimension of array D.  LDD >= MAX(1,P). */

/*     NS      (output) INTEGER */
/*             The dimension of the ALPHA-stable subsystem. */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, the leading NS elements of HSV contain the */
/*             Hankel singular values of the phase system corresponding */
/*             to the ALPHA-stable part of the original system. */
/*             The Hankel singular values are ordered decreasingly. */

/*     Tolerances */

/*     TOL1    DOUBLE PRECISION */
/*             If ORDSEL = 'A', TOL1 contains the tolerance for */
/*             determining the order of reduced system. */
/*             For model reduction, the recommended value of TOL1 lies */
/*             in the interval [0.00001,0.001]. */
/*             If TOL1 <= 0 on entry, the used default value is */
/*             TOL1 = NS*EPS, where NS is the number of */
/*             ALPHA-stable eigenvalues of A and EPS is the machine */
/*             precision (see LAPACK Library Routine DLAMCH). */
/*             If ORDSEL = 'F', the value of TOL1 is ignored. */
/*             TOL1 < 1. */

/*     TOL2    DOUBLE PRECISION */
/*             The tolerance for determining the order of a minimal */
/*             realization of the phase system (see METHOD) corresponding */
/*             to the ALPHA-stable part of the given system. */
/*             The recommended value is TOL2 = NS*EPS. */
/*             This value is used by default if TOL2 <= 0 on entry. */
/*             If TOL2 > 0 and ORDSEL = 'A', then TOL2 <= TOL1. */
/*             TOL2 < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension MAX(1,2*N) */
/*             On exit with INFO = 0, IWORK(1) contains the order of the */
/*             minimal realization of the system. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK and DWORK(2) contains RCOND, the reciprocal */
/*             condition number of the U11 matrix from the expression */
/*             used to compute the solution X = U21*inv(U11) of the */
/*             Riccati equation for spectral factorization. */
/*             A small value RCOND indicates possible ill-conditioning */
/*             of the respective Riccati equation. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 2*N*N + MB*(N+P) + MAX( 2, N*(MAX(N,MB,P)+5), */
/*                                    2*N*P+MAX(P*(MB+2),10*N*(N+1) ) ), */
/*             where MB = M if BETA = 0 and MB = M+P if BETA > 0. */
/*             For optimum performance LDWORK should be larger. */

/*     BWORK   LOGICAL array, dimension 2*N */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 1:  with ORDSEL = 'F', the selected order NR is greater */
/*                   than NSMIN, the sum of the order of the */
/*                   ALPHA-unstable part and the order of a minimal */
/*                   realization of the ALPHA-stable part of the given */
/*                   system; in this case, the resulting NR is set equal */
/*                   to NSMIN; */
/*             = 2:  with ORDSEL = 'F', the selected order NR corresponds */
/*                   to repeated singular values for the ALPHA-stable */
/*                   part, which are neither all included nor all */
/*                   excluded from the reduced model; in this case, the */
/*                   resulting NR is automatically decreased to exclude */
/*                   all repeated singular values; */
/*             = 3:  with ORDSEL = 'F', the selected order NR is less */
/*                   than the order of the ALPHA-unstable part of the */
/*                   given system; in this case NR is set equal to the */
/*                   order of the ALPHA-unstable part. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the computation of the ordered real Schur form of A */
/*                   failed; */
/*             = 2:  the reduction of the Hamiltonian matrix to real */
/*                   Schur form failed; */
/*             = 3:  the reordering of the real Schur form of the */
/*                   Hamiltonian matrix failed; */
/*             = 4:  the Hamiltonian matrix has less than N stable */
/*                   eigenvalues; */
/*             = 5:  the coefficient matrix U11 in the linear system */
/*                   X*U11 = U21 to determine X is singular to working */
/*                   precision; */
/*             = 6:  BETA = 0 and D has not a maximal row rank; */
/*             = 7:  the computation of Hankel singular values failed; */
/*             = 8:  the separation of the ALPHA-stable/unstable diagonal */
/*                   blocks failed because of very close eigenvalues; */
/*             = 9:  the resulting order of reduced stable part is less */
/*                   than the number of unstable zeros of the stable */
/*                   part. */
/*     METHOD */

/*     Let be the following linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t) + Du(t),                      (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system. The subroutine AB09HD determines for */
/*     the given system (1), the matrices of a reduced order system */

/*          d[z(t)] = Ar*z(t) + Br*u(t) */
/*          yr(t)   = Cr*z(t) + Dr*u(t),                  (2) */

/*     such that */

/*          INFNORM[inv(conj(W))*(G-Gr)] <= */
/*                       (1+HSV(NR+NS-N+1)) / (1-HSV(NR+NS-N+1)) + ... */
/*                       + (1+HSV(NS)) / (1-HSV(NS)) - 1, */

/*     where G and Gr are transfer-function matrices of the systems */
/*     (A,B,C,D) and (Ar,Br,Cr,Dr), respectively, W is the right, minimum */
/*     phase spectral factor satisfying */

/*         G1*conj(G1) = conj(W)* W,                      (3) */

/*     G1 is the NS-order ALPHA-stable part of G, and INFNORM(G) is the */
/*     infinity-norm of G. HSV(1), ... , HSV(NS) are the Hankel-singular */
/*     values of the stable part of the phase system (Ap,Bp,Cp) */
/*     with the transfer-function matrix */

/*          P = inv(conj(W))*G1. */

/*     If BETA > 0, then the model reduction is performed on [G BETA*I] */
/*     instead of G. This is the recommended approach to be used when D */
/*     has not a maximal row rank or when a certain balance between */
/*     relative and absolute approximation errors is desired. For */
/*     increasingly large values of BETA, the obtained reduced system */
/*     assymptotically approaches that computed by using the */
/*     Balance & Truncate or Singular Perturbation Approximation methods. */

/*     Note: conj(G)  denotes either G'(-s) for a continuous-time system */
/*           or G'(1/z) for a discrete-time system. */
/*           inv(G) is the inverse of G. */

/*     The following procedure is used to reduce a given G: */

/*     1) Decompose additively G as */

/*          G = G1 + G2, */

/*        such that G1 = (As,Bs,Cs,D) has only ALPHA-stable poles and */
/*        G2 = (Au,Bu,Cu) has only ALPHA-unstable poles. */

/*     2) Determine G1r, a reduced order approximation of the */
/*        ALPHA-stable part G1 using the balancing stochastic method */
/*        in conjunction with either the B&T [1,2] or SPA methods [3]. */

/*     3) Assemble the reduced model Gr as */

/*           Gr = G1r + G2. */

/*     Note: The employed stochastic truncation algorithm [2,3] has the */
/*     property that right half plane zeros of G1 remain as right half */
/*     plane zeros of G1r. Thus, the order can not be chosen smaller than */
/*     the sum of the number of unstable poles of G and the number of */
/*     unstable zeros of G1. */

/*     The reduction of the ALPHA-stable part G1 is done as follows. */

/*     If JOB = 'B', the square-root stochastic Balance & Truncate */
/*     method of [1] is used. */
/*     For an ALPHA-stable continuous-time system (DICO = 'C'), */
/*     the resulting reduced model is stochastically balanced. */

/*     If JOB = 'F', the balancing-free square-root version of the */
/*     stochastic Balance & Truncate method [1] is used to reduce */
/*     the ALPHA-stable part G1. */

/*     If JOB = 'S', the stochastic balancing method is used to reduce */
/*     the ALPHA-stable part G1, in conjunction with the square-root */
/*     version of the Singular Perturbation Approximation method [3,4]. */

/*     If JOB = 'P', the stochastic balancing method is used to reduce */
/*     the ALPHA-stable part G1, in conjunction with the balancing-free */
/*     square-root version of the Singular Perturbation Approximation */
/*     method [3,4]. */

/*     REFERENCES */

/*     [1] Varga A. and Fasol K.H. */
/*         A new square-root balancing-free stochastic truncation model */
/*         reduction algorithm. */
/*         Proc. 12th IFAC World Congress, Sydney, 1993. */

/*     [2] Safonov M. G. and Chiang R. Y. */
/*         Model reduction for robust control: a Schur relative error */
/*         method. */
/*         Int. J. Adapt. Contr. Sign. Proc., vol. 2, pp. 259-272, 1988. */

/*     [3] Green M. and Anderson B. D. O. */
/*         Generalized balanced stochastic truncation. */
/*         Proc. 29-th CDC, Honolulu, Hawaii, pp. 476-481, 1990. */

/*     [4] Varga A. */
/*         Balancing-free square-root algorithm for computing */
/*         singular perturbation approximations. */
/*         Proc. 30-th IEEE CDC,  Brighton, Dec. 11-13, 1991, */
/*         Vol. 2, pp. 1062-1065. */

/*     NUMERICAL ASPECTS */

/*     The implemented methods rely on accuracy enhancing square-root or */
/*     balancing-free square-root techniques. The effectiveness of the */
/*     accuracy enhancing technique depends on the accuracy of the */
/*     solution of a Riccati equation. An ill-conditioned Riccati */
/*     solution typically results when [D BETA*I] is nearly */
/*     rank deficient. */
/*                                      3 */
/*     The algorithm requires about 100N  floating point operations. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, May 2000. */
/*     D. Sima, University of Bucharest, May 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, May 2000. */
/*     Partly based on the RASP routine SRBFS, by A. Varga, 1992. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2000. */
/*              Oct. 2001. */

/*     KEYWORDS */

/*     Minimal realization, model reduction, multivariable system, */
/*     state-space model, state-space representation, */
/*     stochastic balancing. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 436 "AB09HD.f"
    /* Parameter adjustments */
#line 436 "AB09HD.f"
    a_dim1 = *lda;
#line 436 "AB09HD.f"
    a_offset = 1 + a_dim1;
#line 436 "AB09HD.f"
    a -= a_offset;
#line 436 "AB09HD.f"
    b_dim1 = *ldb;
#line 436 "AB09HD.f"
    b_offset = 1 + b_dim1;
#line 436 "AB09HD.f"
    b -= b_offset;
#line 436 "AB09HD.f"
    c_dim1 = *ldc;
#line 436 "AB09HD.f"
    c_offset = 1 + c_dim1;
#line 436 "AB09HD.f"
    c__ -= c_offset;
#line 436 "AB09HD.f"
    d_dim1 = *ldd;
#line 436 "AB09HD.f"
    d_offset = 1 + d_dim1;
#line 436 "AB09HD.f"
    d__ -= d_offset;
#line 436 "AB09HD.f"
    --hsv;
#line 436 "AB09HD.f"
    --iwork;
#line 436 "AB09HD.f"
    --dwork;
#line 436 "AB09HD.f"
    --bwork;
#line 436 "AB09HD.f"

#line 436 "AB09HD.f"
    /* Function Body */
#line 436 "AB09HD.f"
    *info = 0;
#line 437 "AB09HD.f"
    *iwarn = 0;
#line 438 "AB09HD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 439 "AB09HD.f"
    fixord = lsame_(ordsel, "F", (ftnlen)1, (ftnlen)1);
#line 440 "AB09HD.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 441 "AB09HD.f"
    bta = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || lsame_(job, "F", (ftnlen)
	    1, (ftnlen)1);
#line 442 "AB09HD.f"
    spa = lsame_(job, "S", (ftnlen)1, (ftnlen)1) || lsame_(job, "P", (ftnlen)
	    1, (ftnlen)1);
#line 443 "AB09HD.f"
    mb = *m;
#line 444 "AB09HD.f"
    if (*beta > 0.) {
#line 444 "AB09HD.f"
	mb = *m + *p;
#line 444 "AB09HD.f"
    }
/* Computing MAX */
/* Computing MAX */
#line 445 "AB09HD.f"
    i__3 = max(*n,mb);
/* Computing MAX */
#line 445 "AB09HD.f"
    i__4 = *p * (mb + 2), i__5 = *n * 10 * (*n + 1);
#line 445 "AB09HD.f"
    i__1 = 2, i__2 = *n * (max(i__3,*p) + 5), i__1 = max(i__1,i__2), i__2 = (*
	    n << 1) * *p + max(i__4,i__5);
#line 445 "AB09HD.f"
    lw = (*n << 1) * *n + mb * (*n + *p) + max(i__1,i__2);

/*     Test the input scalar arguments. */

#line 450 "AB09HD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 451 "AB09HD.f"
	*info = -1;
#line 452 "AB09HD.f"
    } else if (! (bta || spa)) {
#line 453 "AB09HD.f"
	*info = -2;
#line 454 "AB09HD.f"
    } else if (! (lequil || lsame_(equil, "N", (ftnlen)1, (ftnlen)1))) {
#line 455 "AB09HD.f"
	*info = -3;
#line 456 "AB09HD.f"
    } else if (! (fixord || lsame_(ordsel, "A", (ftnlen)1, (ftnlen)1))) {
#line 457 "AB09HD.f"
	*info = -4;
#line 458 "AB09HD.f"
    } else if (*n < 0) {
#line 459 "AB09HD.f"
	*info = -5;
#line 460 "AB09HD.f"
    } else if (*m < 0) {
#line 461 "AB09HD.f"
	*info = -6;
#line 462 "AB09HD.f"
    } else if (*p < 0 || *beta == 0. && *p > *m) {
#line 463 "AB09HD.f"
	*info = -7;
#line 464 "AB09HD.f"
    } else if (fixord && (*nr < 0 || *nr > *n)) {
#line 465 "AB09HD.f"
	*info = -8;
#line 466 "AB09HD.f"
    } else if (discr && (*alpha < 0. || *alpha > 1.) || ! discr && *alpha > 
	    0.) {
#line 468 "AB09HD.f"
	*info = -9;
#line 469 "AB09HD.f"
    } else if (*beta < 0.) {
#line 470 "AB09HD.f"
	*info = -10;
#line 471 "AB09HD.f"
    } else if (*lda < max(1,*n)) {
#line 472 "AB09HD.f"
	*info = -12;
#line 473 "AB09HD.f"
    } else if (*ldb < max(1,*n)) {
#line 474 "AB09HD.f"
	*info = -14;
#line 475 "AB09HD.f"
    } else if (*ldc < max(1,*p)) {
#line 476 "AB09HD.f"
	*info = -16;
#line 477 "AB09HD.f"
    } else if (*ldd < max(1,*p)) {
#line 478 "AB09HD.f"
	*info = -18;
#line 479 "AB09HD.f"
    } else if (*tol1 >= 1.) {
#line 480 "AB09HD.f"
	*info = -21;
#line 481 "AB09HD.f"
    } else if (*tol2 > 0. && ! fixord && *tol2 > *tol1 || *tol2 >= 1.) {
#line 483 "AB09HD.f"
	*info = -22;
#line 484 "AB09HD.f"
    } else if (*ldwork < lw) {
#line 485 "AB09HD.f"
	*info = -25;
#line 486 "AB09HD.f"
    }

#line 488 "AB09HD.f"
    if (*info != 0) {

/*        Error return. */

#line 492 "AB09HD.f"
	i__1 = -(*info);
#line 492 "AB09HD.f"
	xerbla_("AB09HD", &i__1, (ftnlen)6);
#line 493 "AB09HD.f"
	return 0;
#line 494 "AB09HD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 498 "AB09HD.f"
    i__1 = min(*n,*m);
#line 498 "AB09HD.f"
    if (min(i__1,*p) == 0 || bta && fixord && *nr == 0) {
#line 500 "AB09HD.f"
	*nr = 0;
#line 501 "AB09HD.f"
	*ns = 0;
#line 502 "AB09HD.f"
	iwork[1] = 0;
#line 503 "AB09HD.f"
	dwork[1] = 2.;
#line 504 "AB09HD.f"
	dwork[2] = 1.;
#line 505 "AB09HD.f"
	return 0;
#line 506 "AB09HD.f"
    }

#line 508 "AB09HD.f"
    if (lequil) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D, B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */
/*        Workspace: N. */

#line 515 "AB09HD.f"
	maxred = 100.;
#line 516 "AB09HD.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 518 "AB09HD.f"
    }

/*     Allocate working storage. */

#line 522 "AB09HD.f"
    nn = *n * *n;
#line 523 "AB09HD.f"
    ku = 1;
#line 524 "AB09HD.f"
    kwr = ku + nn;
#line 525 "AB09HD.f"
    kwi = kwr + *n;
#line 526 "AB09HD.f"
    kw = kwi + *n;
#line 527 "AB09HD.f"
    lwr = *ldwork - kw + 1;

/*     Reduce A to a block-diagonal real Schur form, with the */
/*     ALPHA-unstable part in the leading diagonal position, using a */
/*     non-orthogonal similarity transformation A <- inv(T)*A*T and */
/*     apply the transformation to B and C: B <- inv(T)*B and C <- C*T. */

/*     Workspace needed:      N*(N+2); */
/*     Additional workspace:  need   3*N; */
/*                            prefer larger. */

#line 538 "AB09HD.f"
    tb01kd_(dico, "Unstable", "General", n, m, p, alpha, &a[a_offset], lda, &
	    b[b_offset], ldb, &c__[c_offset], ldc, &nu, &dwork[ku], n, &dwork[
	    kwr], &dwork[kwi], &dwork[kw], &lwr, &ierr, (ftnlen)1, (ftnlen)8, 
	    (ftnlen)7);

#line 542 "AB09HD.f"
    if (ierr != 0) {
#line 543 "AB09HD.f"
	if (ierr != 3) {
#line 544 "AB09HD.f"
	    *info = 1;
#line 545 "AB09HD.f"
	} else {
#line 546 "AB09HD.f"
	    *info = 8;
#line 547 "AB09HD.f"
	}
#line 548 "AB09HD.f"
	return 0;
#line 549 "AB09HD.f"
    }

#line 551 "AB09HD.f"
    wrkopt = (integer) dwork[kw] + kw - 1;

#line 553 "AB09HD.f"
    iwarnl = 0;
#line 554 "AB09HD.f"
    *ns = *n - nu;
#line 555 "AB09HD.f"
    if (fixord) {
/* Computing MAX */
#line 556 "AB09HD.f"
	i__1 = 0, i__2 = *nr - nu;
#line 556 "AB09HD.f"
	nra = max(i__1,i__2);
#line 557 "AB09HD.f"
	if (*nr < nu) {
#line 557 "AB09HD.f"
	    iwarnl = 3;
#line 557 "AB09HD.f"
	}
#line 559 "AB09HD.f"
    } else {
#line 560 "AB09HD.f"
	nra = 0;
#line 561 "AB09HD.f"
    }

/*     Finish if the system is completely unstable. */

#line 565 "AB09HD.f"
    if (*ns == 0) {
#line 566 "AB09HD.f"
	*nr = nu;
#line 567 "AB09HD.f"
	iwork[1] = *ns;
#line 568 "AB09HD.f"
	dwork[1] = (doublereal) wrkopt;
#line 569 "AB09HD.f"
	dwork[2] = 1.;
#line 570 "AB09HD.f"
	return 0;
#line 571 "AB09HD.f"
    }

#line 573 "AB09HD.f"
    nu1 = nu + 1;

/*     Allocate working storage. */

#line 577 "AB09HD.f"
    n2 = *n + *n;
#line 578 "AB09HD.f"
    kb = 1;
#line 579 "AB09HD.f"
    kd = kb + *n * mb;
#line 580 "AB09HD.f"
    kt = kd + *p * mb;
#line 581 "AB09HD.f"
    kti = kt + *n * *n;
#line 582 "AB09HD.f"
    kw = kti + *n * *n;

/*     Form [B 0] and [D BETA*I]. */

#line 586 "AB09HD.f"
    dlacpy_("F", ns, m, &b[nu1 + b_dim1], ldb, &dwork[kb], n, (ftnlen)1);
#line 587 "AB09HD.f"
    dlacpy_("F", p, m, &d__[d_offset], ldd, &dwork[kd], p, (ftnlen)1);
#line 588 "AB09HD.f"
    if (*beta > 0.) {
#line 589 "AB09HD.f"
	dlaset_("F", ns, p, &c_b19, &c_b19, &dwork[kb + *n * *m], n, (ftnlen)
		1);
#line 590 "AB09HD.f"
	dlaset_("F", p, p, &c_b19, beta, &dwork[kd + *p * *m], p, (ftnlen)1);
#line 591 "AB09HD.f"
    }

/*     For discrete-time case, apply the discrete-to-continuous bilinear */
/*     transformation to the stable part. */

#line 596 "AB09HD.f"
    if (discr) {

/*        Real workspace:    need  N, prefer larger; */
/*        Integer workspace: need  N. */

#line 601 "AB09HD.f"
	i__1 = *ldwork - kt + 1;
#line 601 "AB09HD.f"
	ab04md_("Discrete", ns, &mb, p, &c_b24, &c_b24, &a[nu1 + nu1 * a_dim1]
		, lda, &dwork[kb], n, &c__[nu1 * c_dim1 + 1], ldc, &dwork[kd],
		 p, &iwork[1], &dwork[kt], &i__1, &ierr, (ftnlen)8);
/* Computing MAX */
#line 604 "AB09HD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[kt] + kt - 1;
#line 604 "AB09HD.f"
	wrkopt = max(i__1,i__2);
#line 605 "AB09HD.f"
    }

/*     Compute in DWORK(KTI) and DWORK(KT) the Cholesky factors S and R */
/*     of the controllability and observability Grammians, respectively. */
/*     Real workspace:    need  2*N*N + MB*(N+P)+ */
/*                              MAX( 2, N*(MAX(N,MB,P)+5), */
/*                                   2*N*P+MAX(P*(MB+2), 10*N*(N+1) ) ); */
/*                        prefer larger. */
/*     Integer workspace: need  2*N. */

#line 615 "AB09HD.f"
    i__1 = *ldwork - kw + 1;
#line 615 "AB09HD.f"
    ab09hy_(ns, &mb, p, &a[nu1 + nu1 * a_dim1], lda, &dwork[kb], n, &c__[nu1 *
	     c_dim1 + 1], ldc, &dwork[kd], p, &scalec, &scaleo, &dwork[kti], 
	    n, &dwork[kt], n, &iwork[1], &dwork[kw], &i__1, &bwork[1], info);
#line 619 "AB09HD.f"
    if (*info != 0) {
#line 619 "AB09HD.f"
	return 0;
#line 619 "AB09HD.f"
    }
/* Computing MAX */
#line 621 "AB09HD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 621 "AB09HD.f"
    wrkopt = max(i__1,i__2);
#line 622 "AB09HD.f"
    ricond = dwork[kw + 1];

/*     Compute a BTA or SPA of the stable part. */
/*     Real workspace:  need  2*N*N + MB*(N+P)+ */
/*                            MAX( 1, 2*N*N+5*N, N*MAX(MB,P) ). */

#line 628 "AB09HD.f"
    epsm = dlamch_("Epsilon", (ftnlen)7);
/* Computing MAX */
#line 629 "AB09HD.f"
    d__2 = *tol1, d__3 = *n * epsm;
#line 629 "AB09HD.f"
    d__1 = max(d__2,d__3);
#line 629 "AB09HD.f"
    i__1 = *ldwork - kw + 1;
#line 629 "AB09HD.f"
    ab09ix_("C", job, "Schur", ordsel, ns, &mb, p, &nra, &scalec, &scaleo, &a[
	    nu1 + nu1 * a_dim1], lda, &dwork[kb], n, &c__[nu1 * c_dim1 + 1], 
	    ldc, &dwork[kd], p, &dwork[kti], n, &dwork[kt], n, &nmr, &hsv[1], 
	    &d__1, tol2, &iwork[1], &dwork[kw], &i__1, iwarn, &ierr, (ftnlen)
	    1, (ftnlen)1, (ftnlen)5, (ftnlen)1);
#line 634 "AB09HD.f"
    *iwarn = max(*iwarn,iwarnl);
#line 635 "AB09HD.f"
    if (ierr != 0) {
#line 636 "AB09HD.f"
	*info = 7;
#line 637 "AB09HD.f"
	return 0;
#line 638 "AB09HD.f"
    }
/* Computing MAX */
#line 639 "AB09HD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 639 "AB09HD.f"
    wrkopt = max(i__1,i__2);

/*     Check if the resulting order is greater than the number of */
/*     unstable zeros (this check is implicit by looking at Hankel */
/*     singular values equal to 1). */

#line 645 "AB09HD.f"
    if (nra < *ns && hsv[nra + 1] >= 1. - pow_dd(&epsm, &c_b29)) {
#line 646 "AB09HD.f"
	*info = 9;
#line 647 "AB09HD.f"
	return 0;
#line 648 "AB09HD.f"
    }

/*     For discrete-time case, apply the continuous-to-discrete */
/*     bilinear transformation. */

#line 653 "AB09HD.f"
    if (discr) {
#line 654 "AB09HD.f"
	ab04md_("Continuous", &nra, &mb, p, &c_b24, &c_b24, &a[nu1 + nu1 * 
		a_dim1], lda, &dwork[kb], n, &c__[nu1 * c_dim1 + 1], ldc, &
		dwork[kd], p, &iwork[1], &dwork[1], ldwork, &ierr, (ftnlen)10)
		;
/* Computing MAX */
#line 657 "AB09HD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[1];
#line 657 "AB09HD.f"
	wrkopt = max(i__1,i__2);
#line 658 "AB09HD.f"
    }

#line 660 "AB09HD.f"
    dlacpy_("F", &nra, m, &dwork[kb], n, &b[nu1 + b_dim1], ldb, (ftnlen)1);
#line 661 "AB09HD.f"
    dlacpy_("F", p, m, &dwork[kd], p, &d__[d_offset], ldd, (ftnlen)1);

#line 663 "AB09HD.f"
    *nr = nra + nu;

#line 665 "AB09HD.f"
    iwork[1] = nmr;
#line 666 "AB09HD.f"
    dwork[1] = (doublereal) wrkopt;
#line 667 "AB09HD.f"
    dwork[2] = ricond;

#line 669 "AB09HD.f"
    return 0;
/* *** Last line of AB09HD *** */
} /* ab09hd_ */

