#line 1 "AB13DD.f"
/* AB13DD.f -- translated by f2c (version 20100827).
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

#line 1 "AB13DD.f"
/* Table of constant values */

static doublereal c_b20 = 1.;
static doublereal c_b21 = 0.;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b163 = -1.;

/* Subroutine */ int ab13dd_(char *dico, char *jobe, char *equil, char *jobd, 
	integer *n, integer *m, integer *p, doublereal *fpeak, doublereal *a, 
	integer *lda, doublereal *e, integer *lde, doublereal *b, integer *
	ldb, doublereal *c__, integer *ldc, doublereal *d__, integer *ldd, 
	doublereal *gpeak, doublereal *tol, integer *iwork, doublereal *dwork,
	 integer *ldwork, doublecomplex *cwork, integer *lcwork, integer *
	info, ftnlen dico_len, ftnlen jobe_len, ftnlen equil_len, ftnlen 
	jobd_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, e_dim1, e_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal), log(doublereal), 
	    sin(doublereal), cos(doublereal), atan(doublereal);

    /* Local variables */
    static integer i__, j, k, n2, ia, ib, ic, id, ie, ih, ii, im;
    static doublereal pi;
    static integer ir, is, nn, iu, iv, pm;
    static doublereal tm;
    static integer lw, ny, ih12, ihi, ipa, iar, ias, ibs, ibt, ipe, ibv, icu, 
	    ies, ilo, isb, isc, isl, nei;
    static doublereal eps, rat;
    static integer nws, n2pm, imin;
    static doublereal anrm;
    static char vect[1];
    static integer ierr, itau, iter;
    static doublereal enrm, temp[1];
    static integer iwrk;
    static doublereal wmax, rtol;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    tg01ad_(char *, integer *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, ftnlen), tg01bd_(char *, 
	    char *, char *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static doublereal gamma;
    extern doublereal ab13dx_(char *, char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublecomplex *, integer *, integer *, ftnlen, ftnlen, ftnlen);
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal omega;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     mb01sd_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, ftnlen), dggev_(char *, char *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), mb03xd_(char *, char *, char *, char *, integer *
	    , doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    static doublereal rcond;
    static logical fulle;
    static doublereal bound, bnorm, cnorm;
    static logical withd;
    static integer minpm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *);
    static logical nodyn;
    static doublereal toler, wrmin;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen);
    extern doublereal dlapy2_(doublereal *, doublereal *);
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *), dgebal_(
	    char *, integer *, doublereal *, integer *, integer *, integer *, 
	    doublereal *, integer *, ftnlen), dggbal_(char *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgehrd_(integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *);
    static doublereal gammal, fpeaki;
    extern /* Subroutine */ int dlascl_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen);
    static doublereal gammas;
    static logical ilascl;
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static doublereal fpeaks;
    static logical ilescl;
    extern /* Subroutine */ int dgesvd_(char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen);
    static doublereal safmin;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen);
    static doublereal safmax, maxred, bignum;
    extern /* Subroutine */ int dhgeqz_(char *, char *, char *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), xerbla_(char *, integer *, 
	    ftnlen), dhseqr_(char *, char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *, integer *, integer *, ftnlen, ftnlen), 
	    dlasrt_(char *, integer *, doublereal *, integer *, ftnlen);
    static integer maxcwk;
    static logical lequil;
    extern /* Subroutine */ int dormhr_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen), dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dorgqr_(integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, integer *);
    static logical usepen;
    static integer mincwr;
    static doublereal anrmto, enrmto;
    static integer minwrk, maxwrk;
    static doublereal smlnum;


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

/*     To compute the L-infinity norm of a continuous-time or */
/*     discrete-time system, either standard or in the descriptor form, */

/*                                     -1 */
/*        G(lambda) = C*( lambda*E - A ) *B + D . */

/*     The norm is finite if and only if the matrix pair (A,E) has no */
/*     eigenvalue on the boundary of the stability domain, i.e., the */
/*     imaginary axis, or the unit circle, respectively. It is assumed */
/*     that the matrix E is nonsingular. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system, as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOBE    CHARACTER*1 */
/*             Specifies whether E is a general square or an identity */
/*             matrix, as follows: */
/*             = 'G':  E is a general square matrix; */
/*             = 'I':  E is the identity matrix. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the system (A,E,B,C) or (A,B,C), as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     JOBD    CHARACTER*1 */
/*             Specifies whether or not a non-zero matrix D appears in */
/*             the given state space model: */
/*             = 'D':  D is present; */
/*             = 'Z':  D is assumed a zero matrix. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The column size of the matrix B.  M >= 0. */

/*     P       (input) INTEGER */
/*             The row size of the matrix C.  P >= 0. */

/*     FPEAK   (input/output) DOUBLE PRECISION array, dimension (2) */
/*             On entry, this parameter must contain an estimate of the */
/*             frequency where the gain of the frequency response would */
/*             achieve its peak value. Setting FPEAK(2) = 0 indicates an */
/*             infinite frequency. An accurate estimate could reduce the */
/*             number of iterations of the iterative algorithm. If no */
/*             estimate is available, set FPEAK(1) = 0, and FPEAK(2) = 1. */
/*             FPEAK(1) >= 0, FPEAK(2) >= 0. */
/*             On exit, if INFO = 0, this array contains the frequency */
/*             OMEGA, where the gain of the frequency response achieves */
/*             its peak value GPEAK, i.e., */

/*                 || G ( j*OMEGA ) || = GPEAK ,  if DICO = 'C', or */

/*                         j*OMEGA */
/*                 || G ( e       ) || = GPEAK ,  if DICO = 'D', */

/*             where OMEGA = FPEAK(1), if FPEAK(2) > 0, and OMEGA is */
/*             infinite, if FPEAK(2) = 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state dynamics matrix A. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= max(1,N). */

/*     E       (input) DOUBLE PRECISION array, dimension (LDE,N) */
/*             If JOBE = 'G', the leading N-by-N part of this array must */
/*             contain the descriptor matrix E of the system. */
/*             If JOBE = 'I', then E is assumed to be the identity */
/*             matrix and is not referenced. */

/*     LDE     INTEGER */
/*             The leading dimension of the array E. */
/*             LDE >= MAX(1,N), if JOBE = 'G'; */
/*             LDE >= 1,        if JOBE = 'I'. */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             system input matrix B. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B.  LDB >= max(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             system output matrix C. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C.  LDC >= max(1,P). */

/*     D       (input) DOUBLE PRECISION array, dimension (LDD,M) */
/*             If JOBD = 'D', the leading P-by-M part of this array must */
/*             contain the direct transmission matrix D. */
/*             The array D is not referenced if JOBD = 'Z'. */

/*     LDD     INTEGER */
/*             The leading dimension of array D. */
/*             LDD >= MAX(1,P), if JOBD = 'D'; */
/*             LDD >= 1,        if JOBD = 'Z'. */

/*     GPEAK   (output) DOUBLE PRECISION array, dimension (2) */
/*             The L-infinity norm of the system, i.e., the peak gain */
/*             of the frequency response (as measured by the largest */
/*             singular value in the MIMO case), coded in the same way */
/*             as FPEAK. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             Tolerance used to set the accuracy in determining the */
/*             norm.  0 <= TOL < 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N) */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) contains the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The dimension of the array DWORK. */
/*             LDWORK >= K, where K can be computed using the following */
/*             pseudo-code (or the Fortran code included in the routine) */

/*                d = 6*MIN(P,M); */
/*                c = MAX( 4*MIN(P,M) + MAX(P,M), d ); */
/*                if ( MIN(P,M) = 0 ) then */
/*                   K = 1; */
/*                else if( N = 0 or B = 0 or C = 0 ) then */
/*                   if( JOBD = 'D' ) then */
/*                      K = P*M + c; */
/*                   else */
/*                      K = 1; */
/*                   end */
/*                else */
/*                   if ( DICO = 'D' ) then */
/*                      b = 0;  e = d; */
/*                   else */
/*                      b = N*(N+M);  e = c; */
/*                      if ( JOBD = Z' ) then  b = b + P*M;  end */
/*                   end */
/*                   if ( JOBD = 'D' ) then */
/*                      r = P*M; */
/*                      if ( JOBE = 'I', DICO = 'C', */
/*                           N > 0, B <> 0, C <> 0 ) then */
/*                         K = P*P + M*M; */
/*                         r = r + N*(P+M); */
/*                      else */
/*                         K = 0; */
/*                      end */
/*                      K = K + r + c;  r = r + MIN(P,M); */
/*                   else */
/*                      r = 0;  K = 0; */
/*                   end */
/*                   r = r + N*(N+P+M); */
/*                   if ( JOBE = 'G' ) then */
/*                      r = r + N*N; */
/*                      if ( EQUIL = 'S' ) then */
/*                         K = MAX( K, r + 9*N ); */
/*                      end */
/*                      K = MAX( K, r + 4*N + MAX( M, 2*N*N, N+b+e ) ); */
/*                   else */
/*                      K = MAX( K, r + N + */
/*                                  MAX( M, P, N*N+2*N, 3*N+b+e ) ); */
/*                   end */
/*                   w = 0; */
/*                   if ( JOBE = 'I', DICO = 'C' ) then */
/*                      w = r + 4*N*N + 11*N; */
/*                      if ( JOBD = 'D' ) then */
/*                         w = w + MAX(M,P) + N*(P+M); */
/*                      end */
/*                   end */
/*                   if ( JOBE = 'E' or DICO = 'D' or JOBD = 'D' ) then */
/*                      w = MAX( w, r + 6*N + (2*N+P+M)*(2*N+P+M) + */
/*                               MAX( 2*(N+P+M), 8*N*N + 16*N ) ); */
/*                   end */
/*                   K = MAX( 1, K, w, r + 2*N + e ); */
/*                end */

/*             For good performance, LDWORK must generally be larger. */

/*             An easily computable upper bound is */

/*             K = MAX( 1, 15*N*N + P*P + M*M + (6*N+3)*(P+M) + 4*P*M + */
/*                         N*M + 22*N + 7*MIN(P,M) ). */

/*             The smallest workspace is obtained for DICO = 'C', */
/*             JOBE = 'I', and JOBD = 'Z', namely */

/*             K = MAX( 1, N*N + N*P + N*M + N + */
/*                         MAX( N*N + N*M + P*M + 3*N + c, */
/*                              4*N*N + 10*N ) ). */

/*             for which an upper bound is */

/*             K = MAX( 1, 6*N*N + N*P + 2*N*M + P*M + 11*N + MAX(P,M) + */
/*                         6*MIN(P,M) ). */

/*     CWORK   COMPLEX*16 array, dimension (LCWORK) */
/*             On exit, if INFO = 0, CWORK(1) contains the optimal */
/*             LCWORK. */

/*     LCWORK  INTEGER */
/*             The dimension of the array CWORK. */
/*             LCWORK >= 1,  if N = 0, or B = 0, or C = 0; */
/*             LCWORK >= MAX(1, (N+M)*(N+P) + 2*MIN(P,M) + MAX(P,M)), */
/*                           otherwise. */
/*             For good performance, LCWORK must generally be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the matrix E is (numerically) singular; */
/*             = 2:  the (periodic) QR (or QZ) algorithm for computing */
/*                   eigenvalues did not converge; */
/*             = 3:  the SVD algorithm for computing singular values did */
/*                   not converge; */
/*             = 4:  the tolerance is too small and the algorithm did */
/*                   not converge. */

/*     METHOD */

/*     The routine implements the method presented in [1], with */
/*     extensions and refinements for improving numerical robustness and */
/*     efficiency. Structure-exploiting eigenvalue computations for */
/*     Hamiltonian matrices are used if JOBE = 'I', DICO = 'C', and the */
/*     symmetric matrices to be implicitly inverted are not too ill- */
/*     conditioned. Otherwise, generalized eigenvalue computations are */
/*     used in the iterative algorithm of [1]. */

/*     REFERENCES */

/*     [1] Bruinsma, N.A. and Steinbuch, M. */
/*         A fast algorithm to compute the Hinfinity-norm of a transfer */
/*         function matrix. */
/*         Systems & Control Letters, vol. 14, pp. 287-293, 1990. */

/*     NUMERICAL ASPECTS */

/*     If the algorithm does not converge in MAXIT = 30 iterations */
/*     (INFO = 4), the tolerance must be increased. */

/*     FURTHER COMMENTS */

/*     If the matrix E is singular, other SLICOT Library routines */
/*     could be used before calling AB13DD, for removing the singular */
/*     part of the system. */

/*     CONTRIBUTORS */

/*     D. Sima, University of Bucharest, May 2001. */
/*     V. Sima, Research Institute for Informatics, Bucharest, May 2001. */
/*     Partly based on SLICOT Library routine AB13CD by P.Hr. Petkov, */
/*     D.W. Gu and M.M. Konstantinov. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, June 2001, */
/*     May 2003, Aug. 2005, March 2008, May 2009, Sep. 2009. */

/*     KEYWORDS */

/*     H-infinity optimal control, robust control, system norm. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar parameters. */

#line 371 "AB13DD.f"
    /* Parameter adjustments */
#line 371 "AB13DD.f"
    --fpeak;
#line 371 "AB13DD.f"
    a_dim1 = *lda;
#line 371 "AB13DD.f"
    a_offset = 1 + a_dim1;
#line 371 "AB13DD.f"
    a -= a_offset;
#line 371 "AB13DD.f"
    e_dim1 = *lde;
#line 371 "AB13DD.f"
    e_offset = 1 + e_dim1;
#line 371 "AB13DD.f"
    e -= e_offset;
#line 371 "AB13DD.f"
    b_dim1 = *ldb;
#line 371 "AB13DD.f"
    b_offset = 1 + b_dim1;
#line 371 "AB13DD.f"
    b -= b_offset;
#line 371 "AB13DD.f"
    c_dim1 = *ldc;
#line 371 "AB13DD.f"
    c_offset = 1 + c_dim1;
#line 371 "AB13DD.f"
    c__ -= c_offset;
#line 371 "AB13DD.f"
    d_dim1 = *ldd;
#line 371 "AB13DD.f"
    d_offset = 1 + d_dim1;
#line 371 "AB13DD.f"
    d__ -= d_offset;
#line 371 "AB13DD.f"
    --gpeak;
#line 371 "AB13DD.f"
    --iwork;
#line 371 "AB13DD.f"
    --dwork;
#line 371 "AB13DD.f"
    --cwork;
#line 371 "AB13DD.f"

#line 371 "AB13DD.f"
    /* Function Body */
#line 371 "AB13DD.f"
    n2 = *n << 1;
#line 372 "AB13DD.f"
    nn = *n * *n;
#line 373 "AB13DD.f"
    pm = *p + *m;
#line 374 "AB13DD.f"
    n2pm = n2 + pm;
#line 375 "AB13DD.f"
    minpm = min(*p,*m);
#line 376 "AB13DD.f"
    *info = 0;
#line 377 "AB13DD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 378 "AB13DD.f"
    fulle = lsame_(jobe, "G", (ftnlen)1, (ftnlen)1);
#line 379 "AB13DD.f"
    lequil = lsame_(equil, "S", (ftnlen)1, (ftnlen)1);
#line 380 "AB13DD.f"
    withd = lsame_(jobd, "D", (ftnlen)1, (ftnlen)1);

#line 382 "AB13DD.f"
    if (! (discr || lsame_(dico, "C", (ftnlen)1, (ftnlen)1))) {
#line 383 "AB13DD.f"
	*info = -1;
#line 384 "AB13DD.f"
    } else if (! (fulle || lsame_(jobe, "I", (ftnlen)1, (ftnlen)1))) {
#line 385 "AB13DD.f"
	*info = -2;
#line 386 "AB13DD.f"
    } else if (! (lequil || lsame_(equil, "N", (ftnlen)1, (ftnlen)1))) {
#line 387 "AB13DD.f"
	*info = -3;
#line 388 "AB13DD.f"
    } else if (! (withd || lsame_(jobd, "Z", (ftnlen)1, (ftnlen)1))) {
#line 389 "AB13DD.f"
	*info = -4;
#line 390 "AB13DD.f"
    } else if (*n < 0) {
#line 391 "AB13DD.f"
	*info = -5;
#line 392 "AB13DD.f"
    } else if (*m < 0) {
#line 393 "AB13DD.f"
	*info = -6;
#line 394 "AB13DD.f"
    } else if (*p < 0) {
#line 395 "AB13DD.f"
	*info = -7;
#line 396 "AB13DD.f"
    } else if (min(fpeak[1],fpeak[2]) < 0.) {
#line 397 "AB13DD.f"
	*info = -8;
#line 398 "AB13DD.f"
    } else if (*lda < max(1,*n)) {
#line 399 "AB13DD.f"
	*info = -10;
#line 400 "AB13DD.f"
    } else if (*lde < 1 || fulle && *lde < *n) {
#line 401 "AB13DD.f"
	*info = -12;
#line 402 "AB13DD.f"
    } else if (*ldb < max(1,*n)) {
#line 403 "AB13DD.f"
	*info = -14;
#line 404 "AB13DD.f"
    } else if (*ldc < max(1,*p)) {
#line 405 "AB13DD.f"
	*info = -16;
#line 406 "AB13DD.f"
    } else if (*ldd < 1 || withd && *ldd < *p) {
#line 407 "AB13DD.f"
	*info = -18;
#line 408 "AB13DD.f"
    } else if (*tol < 0. || *tol >= 1.) {
#line 409 "AB13DD.f"
	*info = -20;
#line 410 "AB13DD.f"
    } else {
#line 411 "AB13DD.f"
	bnorm = dlange_("1-norm", n, m, &b[b_offset], ldb, &dwork[1], (ftnlen)
		6);
#line 412 "AB13DD.f"
	cnorm = dlange_("1-norm", p, n, &c__[c_offset], ldc, &dwork[1], (
		ftnlen)6);
#line 413 "AB13DD.f"
	nodyn = *n == 0 || min(bnorm,cnorm) == 0.;
#line 414 "AB13DD.f"
	usepen = fulle || discr;

/*        Compute workspace. */

#line 418 "AB13DD.f"
	id = minpm * 6;
/* Computing MAX */
#line 419 "AB13DD.f"
	i__1 = (minpm << 2) + max(*p,*m);
#line 419 "AB13DD.f"
	ic = max(i__1,id);
#line 420 "AB13DD.f"
	if (minpm == 0) {
#line 421 "AB13DD.f"
	    minwrk = 1;
#line 422 "AB13DD.f"
	} else if (nodyn) {
#line 423 "AB13DD.f"
	    if (withd) {
#line 424 "AB13DD.f"
		minwrk = *p * *m + ic;
#line 425 "AB13DD.f"
	    } else {
#line 426 "AB13DD.f"
		minwrk = 1;
#line 427 "AB13DD.f"
	    }
#line 428 "AB13DD.f"
	} else {
#line 429 "AB13DD.f"
	    if (discr) {
#line 430 "AB13DD.f"
		ib = 0;
#line 431 "AB13DD.f"
		ie = id;
#line 432 "AB13DD.f"
	    } else {
#line 433 "AB13DD.f"
		ib = *n * (*n + *m);
#line 434 "AB13DD.f"
		if (! withd) {
#line 434 "AB13DD.f"
		    ib += *p * *m;
#line 434 "AB13DD.f"
		}
#line 436 "AB13DD.f"
		ie = ic;
#line 437 "AB13DD.f"
	    }
#line 438 "AB13DD.f"
	    if (withd) {
#line 439 "AB13DD.f"
		ir = *p * *m;
#line 440 "AB13DD.f"
		if (! usepen) {
#line 441 "AB13DD.f"
		    minwrk = *p * *p + *m * *m;
#line 442 "AB13DD.f"
		    ir += *n * pm;
#line 443 "AB13DD.f"
		} else {
#line 444 "AB13DD.f"
		    minwrk = 0;
#line 445 "AB13DD.f"
		}
#line 446 "AB13DD.f"
		minwrk = minwrk + ir + ic;
#line 447 "AB13DD.f"
		ir += minpm;
#line 448 "AB13DD.f"
	    } else {
#line 449 "AB13DD.f"
		ir = 0;
#line 450 "AB13DD.f"
		minwrk = 0;
#line 451 "AB13DD.f"
	    }
#line 452 "AB13DD.f"
	    ir += *n * (*n + pm);
#line 453 "AB13DD.f"
	    if (fulle) {
#line 454 "AB13DD.f"
		ir += nn;
#line 455 "AB13DD.f"
		if (lequil) {
/* Computing MAX */
#line 455 "AB13DD.f"
		    i__1 = minwrk, i__2 = ir + *n * 9;
#line 455 "AB13DD.f"
		    minwrk = max(i__1,i__2);
#line 455 "AB13DD.f"
		}
/* Computing MAX */
/* Computing MAX */
#line 457 "AB13DD.f"
		i__3 = *m, i__4 = nn << 1, i__3 = max(i__3,i__4), i__4 = *n + 
			ib + ie;
#line 457 "AB13DD.f"
		i__1 = minwrk, i__2 = ir + (*n << 2) + max(i__3,i__4);
#line 457 "AB13DD.f"
		minwrk = max(i__1,i__2);
#line 459 "AB13DD.f"
	    } else {
/* Computing MAX */
/* Computing MAX */
#line 460 "AB13DD.f"
		i__3 = max(*m,*p), i__4 = nn + n2, i__3 = max(i__3,i__4), 
			i__4 = *n * 3 + ib + ie;
#line 460 "AB13DD.f"
		i__1 = minwrk, i__2 = ir + *n + max(i__3,i__4);
#line 460 "AB13DD.f"
		minwrk = max(i__1,i__2);
#line 462 "AB13DD.f"
	    }
#line 463 "AB13DD.f"
	    lw = 0;
#line 464 "AB13DD.f"
	    if (! usepen) {
#line 465 "AB13DD.f"
		lw = ir + (nn << 2) + *n * 11;
#line 466 "AB13DD.f"
		if (withd) {
#line 466 "AB13DD.f"
		    lw = lw + max(*m,*p) + *n * pm;
#line 466 "AB13DD.f"
		}
#line 468 "AB13DD.f"
	    }
#line 469 "AB13DD.f"
	    if (usepen || withd) {
/* Computing MAX */
/* Computing MAX */
#line 469 "AB13DD.f"
		i__3 = n2pm + pm, i__4 = nn + n2 << 3;
#line 469 "AB13DD.f"
		i__1 = lw, i__2 = ir + *n * 6 + n2pm * n2pm + max(i__3,i__4);
#line 469 "AB13DD.f"
		lw = max(i__1,i__2);
#line 469 "AB13DD.f"
	    }
/* Computing MAX */
#line 472 "AB13DD.f"
	    i__1 = max(1,minwrk), i__1 = max(i__1,lw), i__2 = ir + n2 + ie;
#line 472 "AB13DD.f"
	    minwrk = max(i__1,i__2);
#line 473 "AB13DD.f"
	}

#line 475 "AB13DD.f"
	if (*ldwork < minwrk) {
#line 476 "AB13DD.f"
	    *info = -23;
#line 477 "AB13DD.f"
	} else {
#line 478 "AB13DD.f"
	    if (nodyn) {
#line 479 "AB13DD.f"
		mincwr = 1;
#line 480 "AB13DD.f"
	    } else {
/* Computing MAX */
#line 481 "AB13DD.f"
		i__1 = 1, i__2 = (*n + *m) * (*n + *p) + (minpm << 1) + max(*
			p,*m);
#line 481 "AB13DD.f"
		mincwr = max(i__1,i__2);
#line 483 "AB13DD.f"
	    }
#line 484 "AB13DD.f"
	    if (*lcwork < mincwr) {
#line 484 "AB13DD.f"
		*info = -25;
#line 484 "AB13DD.f"
	    }
#line 486 "AB13DD.f"
	}
#line 487 "AB13DD.f"
    }

#line 489 "AB13DD.f"
    if (*info != 0) {
#line 490 "AB13DD.f"
	i__1 = -(*info);
#line 490 "AB13DD.f"
	xerbla_("AB13DD", &i__1, (ftnlen)6);
#line 491 "AB13DD.f"
	return 0;
#line 492 "AB13DD.f"
    }

/*     Quick return if possible. */

#line 496 "AB13DD.f"
    if (*m == 0 || *p == 0) {
#line 497 "AB13DD.f"
	gpeak[1] = 0.;
#line 498 "AB13DD.f"
	fpeak[1] = 0.;
#line 499 "AB13DD.f"
	gpeak[2] = 1.;
#line 500 "AB13DD.f"
	fpeak[2] = 1.;
#line 501 "AB13DD.f"
	dwork[1] = 1.;
#line 502 "AB13DD.f"
	cwork[1].r = 1., cwork[1].i = 0.;
#line 503 "AB13DD.f"
	return 0;
#line 504 "AB13DD.f"
    }

/*     Determine the maximum singular value of G(infinity) = D . */
/*     If JOBE = 'I' and DICO = 'C', the full SVD of D, D = U*S*V', is */
/*     computed and saved for later use. */

/*     (Note: Comments in the code beginning "Workspace:" describe the */
/*     minimal amount of real workspace needed at that point in the */
/*     code, as well as the preferred amount for good performance. */
/*     NB refers to the optimal block size for the immediately */
/*     following subroutine, as returned by ILAENV.) */

#line 516 "AB13DD.f"
    id = 1;
#line 517 "AB13DD.f"
    if (withd) {
#line 518 "AB13DD.f"
	is = id + *p * *m;
#line 519 "AB13DD.f"
	if (usepen || nodyn) {
#line 520 "AB13DD.f"
	    iu = is + minpm;
#line 521 "AB13DD.f"
	    iv = iu;
#line 522 "AB13DD.f"
	    iwrk = iv;
#line 523 "AB13DD.f"
	    *(unsigned char *)vect = 'N';
#line 524 "AB13DD.f"
	} else {
#line 525 "AB13DD.f"
	    ibv = is + minpm;
#line 526 "AB13DD.f"
	    icu = ibv + *n * *m;
#line 527 "AB13DD.f"
	    iu = icu + *p * *n;
#line 528 "AB13DD.f"
	    iv = iu + *p * *p;
#line 529 "AB13DD.f"
	    iwrk = iv + *m * *m;
#line 530 "AB13DD.f"
	    *(unsigned char *)vect = 'A';
#line 531 "AB13DD.f"
	}

/*        Workspace: need   P*M + MIN(P,M) + V + */
/*                          MAX( 3*MIN(P,M) + MAX(P,M), 5*MIN(P,M) ), */
/*                          where V = N*(M+P) + P*P + M*M, */
/*                                        if JOBE = 'I' and DICO = 'C', */
/*                                        and N > 0, B <> 0, C <> 0, */
/*                                V = 0,  otherwise; */
/*                   prefer larger. */

#line 541 "AB13DD.f"
	dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[id], p, (ftnlen)4);
#line 542 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 542 "AB13DD.f"
	dgesvd_(vect, vect, p, m, &dwork[id], p, &dwork[is], &dwork[iu], p, &
		dwork[iv], m, &dwork[iwrk], &i__1, &ierr, (ftnlen)1, (ftnlen)
		1);
#line 545 "AB13DD.f"
	if (ierr > 0) {
#line 546 "AB13DD.f"
	    *info = 3;
#line 547 "AB13DD.f"
	    return 0;
#line 548 "AB13DD.f"
	}
#line 549 "AB13DD.f"
	gammal = dwork[is];
#line 550 "AB13DD.f"
	maxwrk = (integer) dwork[iwrk] + iwrk - 1;

/*        Restore D for later calculations. */

#line 554 "AB13DD.f"
	dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[id], p, (ftnlen)4);
#line 555 "AB13DD.f"
    } else {
#line 556 "AB13DD.f"
	iwrk = 1;
#line 557 "AB13DD.f"
	gammal = 0.;
#line 558 "AB13DD.f"
	maxwrk = 1;
#line 559 "AB13DD.f"
    }

/*     Quick return if possible. */

#line 563 "AB13DD.f"
    if (nodyn) {
#line 564 "AB13DD.f"
	gpeak[1] = gammal;
#line 565 "AB13DD.f"
	fpeak[1] = 0.;
#line 566 "AB13DD.f"
	gpeak[2] = 1.;
#line 567 "AB13DD.f"
	fpeak[2] = 1.;
#line 568 "AB13DD.f"
	dwork[1] = (doublereal) maxwrk;
#line 569 "AB13DD.f"
	cwork[1].r = 1., cwork[1].i = 0.;
#line 570 "AB13DD.f"
	return 0;
#line 571 "AB13DD.f"
    }

#line 573 "AB13DD.f"
    if (! usepen && withd) {

/*        Standard continuous-time case, D <> 0: Compute B*V and C'*U . */

#line 577 "AB13DD.f"
	dgemm_("No Transpose", "Transpose", n, m, m, &c_b20, &b[b_offset], 
		ldb, &dwork[iv], m, &c_b21, &dwork[ibv], n, (ftnlen)12, (
		ftnlen)9);
#line 579 "AB13DD.f"
	dgemm_("Transpose", "No Transpose", n, p, p, &c_b20, &c__[c_offset], 
		ldc, &dwork[iu], p, &c_b21, &dwork[icu], n, (ftnlen)9, (
		ftnlen)12);

/*        U and V are no longer needed: free their memory space. */
/*        Total workspace here: need   P*M + MIN(P,M) + N*(M+P) */
/*        (JOBE = 'I', DICO = 'C', JOBD = 'D'). */

#line 586 "AB13DD.f"
	iwrk = iu;
#line 587 "AB13DD.f"
    }

/*     Get machine constants. */

#line 591 "AB13DD.f"
    eps = dlamch_("Epsilon", (ftnlen)7);
#line 592 "AB13DD.f"
    safmin = dlamch_("Safe minimum", (ftnlen)12);
#line 593 "AB13DD.f"
    safmax = 1. / safmin;
#line 594 "AB13DD.f"
    dlabad_(&safmin, &safmax);
#line 595 "AB13DD.f"
    smlnum = sqrt(safmin) / dlamch_("Precision", (ftnlen)9);
#line 596 "AB13DD.f"
    bignum = 1. / smlnum;
#line 597 "AB13DD.f"
    toler = sqrt(eps);

/*     Initiate the transformation of the system to an equivalent one, */
/*     to be used for eigenvalue computations. */

/*     Additional workspace: need   N*N + N*M + P*N + 2*N, if JOBE = 'I'; */
/*                                2*N*N + N*M + P*N + 2*N, if JOBE = 'G'. */

#line 605 "AB13DD.f"
    ia = iwrk;
#line 606 "AB13DD.f"
    ie = ia + nn;
#line 607 "AB13DD.f"
    if (fulle) {
#line 608 "AB13DD.f"
	ib = ie + nn;
#line 609 "AB13DD.f"
    } else {
#line 610 "AB13DD.f"
	ib = ie;
#line 611 "AB13DD.f"
    }
#line 612 "AB13DD.f"
    ic = ib + *n * *m;
#line 613 "AB13DD.f"
    ir = ic + *p * *n;
#line 614 "AB13DD.f"
    ii = ir + *n;
#line 615 "AB13DD.f"
    ibt = ii + *n;

#line 617 "AB13DD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)4);
#line 618 "AB13DD.f"
    dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ib], n, (ftnlen)4);
#line 619 "AB13DD.f"
    dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ic], p, (ftnlen)4);

/*     Scale A if maximum element is outside the range [SMLNUM,BIGNUM]. */

#line 623 "AB13DD.f"
    anrm = dlange_("Max", n, n, &dwork[ia], n, &dwork[1], (ftnlen)3);
#line 624 "AB13DD.f"
    ilascl = FALSE_;
#line 625 "AB13DD.f"
    if (anrm > 0. && anrm < smlnum) {
#line 626 "AB13DD.f"
	anrmto = smlnum;
#line 627 "AB13DD.f"
	ilascl = TRUE_;
#line 628 "AB13DD.f"
    } else if (anrm > bignum) {
#line 629 "AB13DD.f"
	anrmto = bignum;
#line 630 "AB13DD.f"
	ilascl = TRUE_;
#line 631 "AB13DD.f"
    }
#line 632 "AB13DD.f"
    if (ilascl) {
#line 632 "AB13DD.f"
	dlascl_("General", &c__0, &c__0, &anrm, &anrmto, n, n, &dwork[ia], n, 
		&ierr, (ftnlen)7);
#line 632 "AB13DD.f"
    }

#line 636 "AB13DD.f"
    if (fulle) {

/*        Descriptor system. */

/*        Additional workspace: need   N. */

#line 642 "AB13DD.f"
	iwrk = ibt + *n;
#line 643 "AB13DD.f"
	dlacpy_("Full", n, n, &e[e_offset], lde, &dwork[ie], n, (ftnlen)4);

/*        Scale E if maximum element is outside the range */
/*        [SMLNUM,BIGNUM]. */

#line 648 "AB13DD.f"
	enrm = dlange_("Max", n, n, &dwork[ie], n, &dwork[1], (ftnlen)3);
#line 649 "AB13DD.f"
	ilescl = FALSE_;
#line 650 "AB13DD.f"
	if (enrm > 0. && enrm < smlnum) {
#line 651 "AB13DD.f"
	    enrmto = smlnum;
#line 652 "AB13DD.f"
	    ilescl = TRUE_;
#line 653 "AB13DD.f"
	} else if (enrm > bignum) {
#line 654 "AB13DD.f"
	    enrmto = bignum;
#line 655 "AB13DD.f"
	    ilescl = TRUE_;
#line 656 "AB13DD.f"
	} else if (enrm == 0.) {

/*           Error return: Matrix E is 0. */

#line 660 "AB13DD.f"
	    *info = 1;
#line 661 "AB13DD.f"
	    return 0;
#line 662 "AB13DD.f"
	}
#line 663 "AB13DD.f"
	if (ilescl) {
#line 663 "AB13DD.f"
	    dlascl_("General", &c__0, &c__0, &enrm, &enrmto, n, n, &dwork[ie],
		     n, &ierr, (ftnlen)7);
#line 663 "AB13DD.f"
	}

/*        Equilibrate the system, if required. */

/*        Additional workspace: need   6*N. */

#line 671 "AB13DD.f"
	if (lequil) {
#line 671 "AB13DD.f"
	    tg01ad_("All", n, n, m, p, &c_b21, &dwork[ia], n, &dwork[ie], n, &
		    dwork[ib], n, &dwork[ic], p, &dwork[ii], &dwork[ir], &
		    dwork[iwrk], &ierr, (ftnlen)3);
#line 671 "AB13DD.f"
	}

/*        For efficiency of later calculations, the system (A,E,B,C) is */
/*        reduced to an equivalent one with the state matrix A in */
/*        Hessenberg form, and E upper triangular. */
/*        First, permute (A,E) to make it more nearly triangular. */

#line 682 "AB13DD.f"
	dggbal_("Permute", n, &dwork[ia], n, &dwork[ie], n, &ilo, &ihi, &
		dwork[ii], &dwork[ir], &dwork[iwrk], &ierr, (ftnlen)7);

/*        Apply the permutations to (the copies of) B and C. */

#line 688 "AB13DD.f"
	i__1 = ihi + 1;
#line 688 "AB13DD.f"
	for (i__ = *n; i__ >= i__1; --i__) {
#line 689 "AB13DD.f"
	    k = (integer) dwork[ii + i__ - 1];
#line 690 "AB13DD.f"
	    if (k != i__) {
#line 690 "AB13DD.f"
		dswap_(m, &dwork[ib + i__ - 1], n, &dwork[ib + k - 1], n);
#line 690 "AB13DD.f"
	    }
#line 693 "AB13DD.f"
	    k = (integer) dwork[ir + i__ - 1];
#line 694 "AB13DD.f"
	    if (k != i__) {
#line 694 "AB13DD.f"
		dswap_(p, &dwork[ic + (i__ - 1) * *p], &c__1, &dwork[ic + (k 
			- 1) * *p], &c__1);
#line 694 "AB13DD.f"
	    }
#line 697 "AB13DD.f"
/* L10: */
#line 697 "AB13DD.f"
	}

#line 699 "AB13DD.f"
	i__1 = ilo - 1;
#line 699 "AB13DD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 700 "AB13DD.f"
	    k = (integer) dwork[ii + i__ - 1];
#line 701 "AB13DD.f"
	    if (k != i__) {
#line 701 "AB13DD.f"
		dswap_(m, &dwork[ib + i__ - 1], n, &dwork[ib + k - 1], n);
#line 701 "AB13DD.f"
	    }
#line 704 "AB13DD.f"
	    k = (integer) dwork[ir + i__ - 1];
#line 705 "AB13DD.f"
	    if (k != i__) {
#line 705 "AB13DD.f"
		dswap_(p, &dwork[ic + (i__ - 1) * *p], &c__1, &dwork[ic + (k 
			- 1) * *p], &c__1);
#line 705 "AB13DD.f"
	    }
#line 708 "AB13DD.f"
/* L20: */
#line 708 "AB13DD.f"
	}

/*        Reduce (A,E) to generalized Hessenberg form and apply the */
/*        transformations to B and C. */
/*        Additional workspace: need   N + MAX(N,M); */
/*                              prefer N + MAX(N,M)*NB. */

#line 715 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 715 "AB13DD.f"
	tg01bd_("General", "No Q", "No Z", n, m, p, &ilo, &ihi, &dwork[ia], n,
		 &dwork[ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[1], &
		c__1, &dwork[1], &c__1, &dwork[iwrk], &i__1, &ierr, (ftnlen)7,
		 (ftnlen)4, (ftnlen)4);
/* Computing MAX */
#line 719 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 719 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);

/*        Check whether matrix E is nonsingular. */
/*        Additional workspace: need   3*N. */

#line 724 "AB13DD.f"
	dtrcon_("1-norm", "Upper", "Non Unit", n, &dwork[ie], n, &rcond, &
		dwork[iwrk], &iwork[1], &ierr, (ftnlen)6, (ftnlen)5, (ftnlen)
		8);
#line 726 "AB13DD.f"
	if (rcond <= (doublereal) (*n) * 10. * eps) {

/*           Error return: Matrix E is numerically singular. */

#line 730 "AB13DD.f"
	    *info = 1;
#line 731 "AB13DD.f"
	    return 0;
#line 732 "AB13DD.f"
	}

/*        Perform QZ algorithm, computing eigenvalues. The generalized */
/*        Hessenberg form is saved for later use. */
/*        Additional workspace: need   2*N*N + N; */
/*                              prefer larger. */

#line 739 "AB13DD.f"
	ias = iwrk;
#line 740 "AB13DD.f"
	ies = ias + nn;
#line 741 "AB13DD.f"
	iwrk = ies + nn;
#line 742 "AB13DD.f"
	dlacpy_("Full", n, n, &dwork[ia], n, &dwork[ias], n, (ftnlen)4);
#line 743 "AB13DD.f"
	dlacpy_("Full", n, n, &dwork[ie], n, &dwork[ies], n, (ftnlen)4);
#line 744 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 744 "AB13DD.f"
	dhgeqz_("Eigenvalues", "No Vectors", "No Vectors", n, &ilo, &ihi, &
		dwork[ias], n, &dwork[ies], n, &dwork[ir], &dwork[ii], &dwork[
		ibt], &dwork[1], n, &dwork[1], n, &dwork[iwrk], &i__1, &ierr, 
		(ftnlen)11, (ftnlen)10, (ftnlen)10);
#line 748 "AB13DD.f"
	if (ierr != 0) {
#line 749 "AB13DD.f"
	    *info = 2;
#line 750 "AB13DD.f"
	    return 0;
#line 751 "AB13DD.f"
	}
/* Computing MAX */
#line 752 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 752 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);

/*        Check if unscaling would cause over/underflow; if so, rescale */
/*        eigenvalues (DWORK( IR+I-1 ),DWORK( II+I-1 ),DWORK( IBT+I-1 )) */
/*        so DWORK( IBT+I-1 ) is on the order of E(I,I) and */
/*        DWORK( IR+I-1 ) and DWORK( II+I-1 ) are on the order of A(I,I). */

#line 759 "AB13DD.f"
	if (ilascl) {

#line 761 "AB13DD.f"
	    i__1 = *n;
#line 761 "AB13DD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 762 "AB13DD.f"
		if (dwork[ii + i__ - 1] != 0.) {
#line 763 "AB13DD.f"
		    if (dwork[ir + i__ - 1] / safmax > anrmto / anrm || 
			    safmin / dwork[ir + i__ - 1] > anrm / anrmto) {
#line 767 "AB13DD.f"
			tm = (d__1 = dwork[ia + (i__ - 1) * *n + i__] / dwork[
				ir + i__ - 1], abs(d__1));
#line 768 "AB13DD.f"
			dwork[ibt + i__ - 1] *= tm;
#line 769 "AB13DD.f"
			dwork[ir + i__ - 1] *= tm;
#line 770 "AB13DD.f"
			dwork[ii + i__ - 1] *= tm;
#line 771 "AB13DD.f"
		    } else if (dwork[ii + i__ - 1] / safmax > anrmto / anrm ||
			     safmin / dwork[ii + i__ - 1] > anrm / anrmto) {
#line 775 "AB13DD.f"
			tm = (d__1 = dwork[ia + i__ * *n + i__] / dwork[ii + 
				i__ - 1], abs(d__1));
#line 776 "AB13DD.f"
			dwork[ibt + i__ - 1] *= tm;
#line 777 "AB13DD.f"
			dwork[ir + i__ - 1] *= tm;
#line 778 "AB13DD.f"
			dwork[ii + i__ - 1] *= tm;
#line 779 "AB13DD.f"
		    }
#line 780 "AB13DD.f"
		}
#line 781 "AB13DD.f"
/* L30: */
#line 781 "AB13DD.f"
	    }

#line 783 "AB13DD.f"
	}

#line 785 "AB13DD.f"
	if (ilescl) {

#line 787 "AB13DD.f"
	    i__1 = *n;
#line 787 "AB13DD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 788 "AB13DD.f"
		if (dwork[ii + i__ - 1] != 0.) {
#line 789 "AB13DD.f"
		    if (dwork[ibt + i__ - 1] / safmax > enrmto / enrm || 
			    safmin / dwork[ibt + i__ - 1] > enrm / enrmto) {
#line 793 "AB13DD.f"
			tm = (d__1 = dwork[ie + (i__ - 1) * *n + i__] / dwork[
				ibt + i__ - 1], abs(d__1));
#line 794 "AB13DD.f"
			dwork[ibt + i__ - 1] *= tm;
#line 795 "AB13DD.f"
			dwork[ir + i__ - 1] *= tm;
#line 796 "AB13DD.f"
			dwork[ii + i__ - 1] *= tm;
#line 797 "AB13DD.f"
		    }
#line 798 "AB13DD.f"
		}
#line 799 "AB13DD.f"
/* L40: */
#line 799 "AB13DD.f"
	    }

#line 801 "AB13DD.f"
	}

/*        Undo scaling. */

#line 805 "AB13DD.f"
	if (ilascl) {
#line 806 "AB13DD.f"
	    dlascl_("Hessenberg", &c__0, &c__0, &anrmto, &anrm, n, n, &dwork[
		    ia], n, &ierr, (ftnlen)10);
#line 808 "AB13DD.f"
	    dlascl_("General", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &dwork[
		    ir], n, &ierr, (ftnlen)7);
#line 810 "AB13DD.f"
	    dlascl_("General", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &dwork[
		    ii], n, &ierr, (ftnlen)7);
#line 812 "AB13DD.f"
	}

#line 814 "AB13DD.f"
	if (ilescl) {
#line 815 "AB13DD.f"
	    dlascl_("Upper", &c__0, &c__0, &enrmto, &enrm, n, n, &dwork[ie], 
		    n, &ierr, (ftnlen)5);
#line 817 "AB13DD.f"
	    dlascl_("General", &c__0, &c__0, &enrmto, &enrm, n, &c__1, &dwork[
		    ibt], n, &ierr, (ftnlen)7);
#line 819 "AB13DD.f"
	}

#line 821 "AB13DD.f"
    } else {

/*        Standard state-space system. */

#line 825 "AB13DD.f"
	if (lequil) {

/*           Equilibrate the system. */

#line 829 "AB13DD.f"
	    maxred = 100.;
#line 830 "AB13DD.f"
	    tb01id_("All", n, m, p, &maxred, &dwork[ia], n, &dwork[ib], n, &
		    dwork[ic], p, &dwork[ii], &ierr, (ftnlen)3);
#line 833 "AB13DD.f"
	}

/*        For efficiency of later calculations, the system (A,B,C) is */
/*        reduced to a similar one with the state matrix in Hessenberg */
/*        form. */

/*        First, permute the matrix A to make it more nearly triangular */
/*        and apply the permutations to B and C. */

#line 842 "AB13DD.f"
	dgebal_("Permute", n, &dwork[ia], n, &ilo, &ihi, &dwork[ir], &ierr, (
		ftnlen)7);

#line 845 "AB13DD.f"
	i__1 = ihi + 1;
#line 845 "AB13DD.f"
	for (i__ = *n; i__ >= i__1; --i__) {
#line 846 "AB13DD.f"
	    k = (integer) dwork[ir + i__ - 1];
#line 847 "AB13DD.f"
	    if (k != i__) {
#line 848 "AB13DD.f"
		dswap_(m, &dwork[ib + i__ - 1], n, &dwork[ib + k - 1], n);
#line 850 "AB13DD.f"
		dswap_(p, &dwork[ic + (i__ - 1) * *p], &c__1, &dwork[ic + (k 
			- 1) * *p], &c__1);
#line 852 "AB13DD.f"
	    }
#line 853 "AB13DD.f"
/* L50: */
#line 853 "AB13DD.f"
	}

#line 855 "AB13DD.f"
	i__1 = ilo - 1;
#line 855 "AB13DD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 856 "AB13DD.f"
	    k = (integer) dwork[ir + i__ - 1];
#line 857 "AB13DD.f"
	    if (k != i__) {
#line 858 "AB13DD.f"
		dswap_(m, &dwork[ib + i__ - 1], n, &dwork[ib + k - 1], n);
#line 860 "AB13DD.f"
		dswap_(p, &dwork[ic + (i__ - 1) * *p], &c__1, &dwork[ic + (k 
			- 1) * *p], &c__1);
#line 862 "AB13DD.f"
	    }
#line 863 "AB13DD.f"
/* L60: */
#line 863 "AB13DD.f"
	}

/*        Reduce A to upper Hessenberg form and apply the transformations */
/*        to B and C. */
/*        Additional workspace: need   N;   (from II) */
/*                              prefer N*NB. */

#line 870 "AB13DD.f"
	itau = ir;
#line 871 "AB13DD.f"
	iwrk = itau + *n;
#line 872 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 872 "AB13DD.f"
	dgehrd_(n, &ilo, &ihi, &dwork[ia], n, &dwork[itau], &dwork[iwrk], &
		i__1, &ierr);
/* Computing MAX */
#line 874 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 874 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);

/*        Additional workspace: need   M; */
/*                              prefer M*NB. */

#line 879 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 879 "AB13DD.f"
	dormhr_("Left", "Transpose", n, m, &ilo, &ihi, &dwork[ia], n, &dwork[
		itau], &dwork[ib], n, &dwork[iwrk], &i__1, &ierr, (ftnlen)4, (
		ftnlen)9);
/* Computing MAX */
#line 882 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 882 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);

/*        Additional workspace: need   P; */
/*                              prefer P*NB. */

#line 887 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 887 "AB13DD.f"
	dormhr_("Right", "NoTranspose", p, n, &ilo, &ihi, &dwork[ia], n, &
		dwork[itau], &dwork[ic], p, &dwork[iwrk], &i__1, &ierr, (
		ftnlen)5, (ftnlen)11);
/* Computing MAX */
#line 890 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 890 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);

/*        Compute the eigenvalues. The Hessenberg form is saved for */
/*        later use. */
/*        Additional workspace:  need   N*N + N;   (from IBT) */
/*                               prefer larger. */

#line 897 "AB13DD.f"
	ias = ibt;
#line 898 "AB13DD.f"
	iwrk = ias + nn;
#line 899 "AB13DD.f"
	dlacpy_("Full", n, n, &dwork[ia], n, &dwork[ias], n, (ftnlen)4);
#line 900 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 900 "AB13DD.f"
	dhseqr_("Eigenvalues", "No Vectors", n, &ilo, &ihi, &dwork[ias], n, &
		dwork[ir], &dwork[ii], &dwork[1], n, &dwork[iwrk], &i__1, &
		ierr, (ftnlen)11, (ftnlen)10);
#line 903 "AB13DD.f"
	if (ierr > 0) {
#line 904 "AB13DD.f"
	    *info = 2;
#line 905 "AB13DD.f"
	    return 0;
#line 906 "AB13DD.f"
	}
/* Computing MAX */
#line 907 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 907 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);

#line 909 "AB13DD.f"
	if (ilascl) {

/*           Undo scaling for the Hessenberg form of A and eigenvalues. */

#line 913 "AB13DD.f"
	    dlascl_("Hessenberg", &c__0, &c__0, &anrmto, &anrm, n, n, &dwork[
		    ia], n, &ierr, (ftnlen)10);
#line 915 "AB13DD.f"
	    dlascl_("General", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &dwork[
		    ir], n, &ierr, (ftnlen)7);
#line 917 "AB13DD.f"
	    dlascl_("General", &c__0, &c__0, &anrmto, &anrm, n, &c__1, &dwork[
		    ii], n, &ierr, (ftnlen)7);
#line 919 "AB13DD.f"
	}

#line 921 "AB13DD.f"
    }

/*     Look for (generalized) eigenvalues on the boundary of the */
/*     stability domain. (Their existence implies an infinite norm.) */
/*     Additional workspace:  need   2*N.   (from IAS) */

#line 927 "AB13DD.f"
    im = ias;
#line 928 "AB13DD.f"
    iar = im + *n;
#line 929 "AB13DD.f"
    imin = ii;
#line 930 "AB13DD.f"
    wrmin = safmax;
#line 931 "AB13DD.f"
    bound = eps * 1e3;

#line 933 "AB13DD.f"
    if (discr) {
#line 934 "AB13DD.f"
	gammal = 0.;

/*        For discrete-time case, compute the logarithm of the non-zero */
/*        eigenvalues and save their moduli and absolute real parts. */
/*        (The logarithms are overwritten on the eigenvalues.) */
/*        Also, find the minimum distance to the unit circle. */

#line 941 "AB13DD.f"
	if (fulle) {

#line 943 "AB13DD.f"
	    i__1 = *n - 1;
#line 943 "AB13DD.f"
	    for (i__ = 0; i__ <= i__1; ++i__) {
#line 944 "AB13DD.f"
		tm = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
#line 945 "AB13DD.f"
		if (dwork[ibt + i__] >= 1. || dwork[ibt + i__] < 1. && tm < 
			safmax * dwork[ibt + i__]) {
#line 948 "AB13DD.f"
		    tm /= dwork[ibt + i__];
#line 949 "AB13DD.f"
		} else {

/*                 The pencil has too large eigenvalues. SAFMAX is used. */

#line 953 "AB13DD.f"
		    tm = safmax;
#line 954 "AB13DD.f"
		}
#line 955 "AB13DD.f"
		if (tm != 0.) {
#line 956 "AB13DD.f"
		    dwork[ii + i__] = atan2(dwork[ii + i__], dwork[ir + i__]);
#line 957 "AB13DD.f"
		    dwork[ir + i__] = log(tm);
#line 958 "AB13DD.f"
		}
#line 959 "AB13DD.f"
		dwork[im] = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
#line 960 "AB13DD.f"
		tm = (d__1 = 1. - tm, abs(d__1));
#line 961 "AB13DD.f"
		if (tm < wrmin) {
#line 962 "AB13DD.f"
		    imin = ii + i__;
#line 963 "AB13DD.f"
		    wrmin = tm;
#line 964 "AB13DD.f"
		}
#line 965 "AB13DD.f"
		++im;
#line 966 "AB13DD.f"
		dwork[iar + i__] = (d__1 = dwork[ir + i__], abs(d__1));
#line 967 "AB13DD.f"
/* L70: */
#line 967 "AB13DD.f"
	    }

#line 969 "AB13DD.f"
	} else {

#line 971 "AB13DD.f"
	    i__1 = *n - 1;
#line 971 "AB13DD.f"
	    for (i__ = 0; i__ <= i__1; ++i__) {
#line 972 "AB13DD.f"
		tm = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
#line 973 "AB13DD.f"
		if (tm != 0.) {
#line 974 "AB13DD.f"
		    dwork[ii + i__] = atan2(dwork[ii + i__], dwork[ir + i__]);
#line 975 "AB13DD.f"
		    dwork[ir + i__] = log(tm);
#line 976 "AB13DD.f"
		}
#line 977 "AB13DD.f"
		dwork[im] = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
#line 978 "AB13DD.f"
		tm = (d__1 = 1. - tm, abs(d__1));
#line 979 "AB13DD.f"
		if (tm < wrmin) {
#line 980 "AB13DD.f"
		    imin = ii + i__;
#line 981 "AB13DD.f"
		    wrmin = tm;
#line 982 "AB13DD.f"
		}
#line 983 "AB13DD.f"
		++im;
#line 984 "AB13DD.f"
		dwork[iar + i__] = (d__1 = dwork[ir + i__], abs(d__1));
#line 985 "AB13DD.f"
/* L80: */
#line 985 "AB13DD.f"
	    }

#line 987 "AB13DD.f"
	}

#line 989 "AB13DD.f"
    } else {

/*        For continuous-time case, save moduli of eigenvalues and */
/*        absolute real parts and find the maximum modulus and minimum */
/*        absolute real part. */

#line 995 "AB13DD.f"
	wmax = 0.;

#line 997 "AB13DD.f"
	if (fulle) {

#line 999 "AB13DD.f"
	    i__1 = *n - 1;
#line 999 "AB13DD.f"
	    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1000 "AB13DD.f"
		tm = (d__1 = dwork[ir + i__], abs(d__1));
#line 1001 "AB13DD.f"
		dwork[im] = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
#line 1002 "AB13DD.f"
		if (dwork[ibt + i__] >= 1. || dwork[ibt + i__] < 1. && dwork[
			im] < safmax * dwork[ibt + i__]) {
#line 1006 "AB13DD.f"
		    tm /= dwork[ibt + i__];
#line 1007 "AB13DD.f"
		    dwork[im] /= dwork[ibt + i__];
#line 1008 "AB13DD.f"
		} else {
#line 1009 "AB13DD.f"
		    if (tm < safmax * dwork[ibt + i__]) {
#line 1010 "AB13DD.f"
			tm /= dwork[ibt + i__];
#line 1011 "AB13DD.f"
		    } else {

/*                    The pencil has too large eigenvalues. */
/*                    SAFMAX is used. */

#line 1016 "AB13DD.f"
			tm = safmax;
#line 1017 "AB13DD.f"
		    }
#line 1018 "AB13DD.f"
		    dwork[im] = safmax;
#line 1019 "AB13DD.f"
		}
#line 1020 "AB13DD.f"
		if (tm < wrmin) {
#line 1021 "AB13DD.f"
		    imin = ii + i__;
#line 1022 "AB13DD.f"
		    wrmin = tm;
#line 1023 "AB13DD.f"
		}
#line 1024 "AB13DD.f"
		dwork[iar + i__] = tm;
#line 1025 "AB13DD.f"
		if (dwork[im] > wmax) {
#line 1025 "AB13DD.f"
		    wmax = dwork[im];
#line 1025 "AB13DD.f"
		}
#line 1027 "AB13DD.f"
		++im;
#line 1028 "AB13DD.f"
/* L90: */
#line 1028 "AB13DD.f"
	    }

#line 1030 "AB13DD.f"
	} else {

#line 1032 "AB13DD.f"
	    i__1 = *n - 1;
#line 1032 "AB13DD.f"
	    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1033 "AB13DD.f"
		tm = (d__1 = dwork[ir + i__], abs(d__1));
#line 1034 "AB13DD.f"
		if (tm < wrmin) {
#line 1035 "AB13DD.f"
		    imin = ii + i__;
#line 1036 "AB13DD.f"
		    wrmin = tm;
#line 1037 "AB13DD.f"
		}
#line 1038 "AB13DD.f"
		dwork[im] = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
#line 1039 "AB13DD.f"
		if (dwork[im] > wmax) {
#line 1039 "AB13DD.f"
		    wmax = dwork[im];
#line 1039 "AB13DD.f"
		}
#line 1041 "AB13DD.f"
		++im;
#line 1042 "AB13DD.f"
		dwork[iar + i__] = tm;
#line 1043 "AB13DD.f"
/* L100: */
#line 1043 "AB13DD.f"
	    }

#line 1045 "AB13DD.f"
	}

#line 1047 "AB13DD.f"
	bound += eps * wmax;

#line 1049 "AB13DD.f"
    }

#line 1051 "AB13DD.f"
    im -= *n;

#line 1053 "AB13DD.f"
    if (wrmin < bound) {

/*        The L-infinity norm was found as infinite. */

#line 1057 "AB13DD.f"
	gpeak[1] = 1.;
#line 1058 "AB13DD.f"
	gpeak[2] = 0.;
#line 1059 "AB13DD.f"
	tm = (d__1 = dwork[imin], abs(d__1));
#line 1060 "AB13DD.f"
	if (discr) {
#line 1060 "AB13DD.f"
	    tm = (d__1 = atan2(sin(tm), cos(tm)), abs(d__1));
#line 1060 "AB13DD.f"
	}
#line 1062 "AB13DD.f"
	fpeak[1] = tm;
#line 1063 "AB13DD.f"
	if (tm < safmax) {
#line 1064 "AB13DD.f"
	    fpeak[2] = 1.;
#line 1065 "AB13DD.f"
	} else {
#line 1066 "AB13DD.f"
	    fpeak[2] = 0.;
#line 1067 "AB13DD.f"
	}

#line 1069 "AB13DD.f"
	dwork[1] = (doublereal) maxwrk;
#line 1070 "AB13DD.f"
	cwork[1].r = 1., cwork[1].i = 0.;
#line 1071 "AB13DD.f"
	return 0;
#line 1072 "AB13DD.f"
    }

/*     Determine the maximum singular value of */
/*        G(lambda) = C*inv(lambda*E - A)*B + D, */
/*     over a selected set of frequencies. Besides the frequencies w = 0, */
/*     w = pi (if DICO = 'D'), and the given value FPEAK, this test set */
/*     contains the peak frequency for each mode (or an approximation */
/*     of it). The (generalized) Hessenberg form of the system is used. */

/*     First, determine the maximum singular value of G(0) and set FPEAK */
/*     accordingly. */
/*     Additional workspace: */
/*           complex: need   1, if DICO = 'C'; */
/*                           (N+M)*(N+P)+2*MIN(P,M)+MAX(P,M)), otherwise; */
/*                    prefer larger; */
/*           real:    need   LDW0+LDW1+LDW2, where */
/*                           LDW0 = N*N+N*M, if DICO = 'C'; */
/*                           LDW0 = 0,       if DICO = 'D'; */
/*                           LDW1 = P*M,     if DICO = 'C', JOBD = 'Z'; */
/*                           LDW1 = 0,       otherwise; */
/*                           LDW2 = MIN(P,M)+MAX(3*MIN(P,M)+MAX(P,M), */
/*                                               5*MIN(P,M)), */
/*                                              if DICO = 'C'; */
/*                           LDW2 = 6*MIN(P,M), otherwise. */
/*                    prefer larger. */

#line 1098 "AB13DD.f"
    if (discr) {
#line 1099 "AB13DD.f"
	ias = ia;
#line 1100 "AB13DD.f"
	ibs = ib;
#line 1101 "AB13DD.f"
	iwrk = iar + *n;
#line 1102 "AB13DD.f"
    } else {
#line 1103 "AB13DD.f"
	ias = iar + *n;
#line 1104 "AB13DD.f"
	ibs = ias + nn;
#line 1105 "AB13DD.f"
	iwrk = ibs + *n * *m;
#line 1106 "AB13DD.f"
	dlacpy_("Upper", n, n, &dwork[ia], n, &dwork[ias], n, (ftnlen)5);
#line 1107 "AB13DD.f"
	i__1 = *n - 1;
#line 1107 "AB13DD.f"
	i__2 = *n + 1;
#line 1107 "AB13DD.f"
	i__3 = *n + 1;
#line 1107 "AB13DD.f"
	dcopy_(&i__1, &dwork[ia + 1], &i__2, &dwork[ias + 1], &i__3);
#line 1108 "AB13DD.f"
	dlacpy_("Full", n, m, &dwork[ib], n, &dwork[ibs], n, (ftnlen)4);
#line 1109 "AB13DD.f"
    }
#line 1110 "AB13DD.f"
    i__1 = *ldwork - iwrk + 1;
#line 1110 "AB13DD.f"
    gamma = ab13dx_(dico, jobe, jobd, n, m, p, &c_b21, &dwork[ias], n, &dwork[
	    ie], n, &dwork[ibs], n, &dwork[ic], p, &dwork[id], p, &iwork[1], &
	    dwork[iwrk], &i__1, &cwork[1], lcwork, &ierr, (ftnlen)1, (ftnlen)
	    1, (ftnlen)1);
/* Computing MAX */
#line 1114 "AB13DD.f"
    i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 1114 "AB13DD.f"
    maxwrk = max(i__1,maxwrk);
#line 1115 "AB13DD.f"
    if (ierr >= 1 && ierr <= *n) {
#line 1116 "AB13DD.f"
	gpeak[1] = 1.;
#line 1117 "AB13DD.f"
	fpeak[1] = 0.;
#line 1118 "AB13DD.f"
	gpeak[2] = 0.;
#line 1119 "AB13DD.f"
	fpeak[2] = 1.;
#line 1120 "AB13DD.f"
	goto L340;
#line 1121 "AB13DD.f"
    } else if (ierr == *n + 1) {
#line 1122 "AB13DD.f"
	*info = 3;
#line 1123 "AB13DD.f"
	return 0;
#line 1124 "AB13DD.f"
    }

#line 1126 "AB13DD.f"
    fpeaks = fpeak[1];
#line 1127 "AB13DD.f"
    fpeaki = fpeak[2];
#line 1128 "AB13DD.f"
    if (gammal < gamma) {
#line 1129 "AB13DD.f"
	gammal = gamma;
#line 1130 "AB13DD.f"
	fpeak[1] = 0.;
#line 1131 "AB13DD.f"
	fpeak[2] = 1.;
#line 1132 "AB13DD.f"
    } else if (! discr) {
#line 1133 "AB13DD.f"
	fpeak[1] = 1.;
#line 1134 "AB13DD.f"
	fpeak[2] = 0.;
#line 1135 "AB13DD.f"
    }

#line 1137 "AB13DD.f"
    maxcwk = (integer) cwork[1].r;

#line 1139 "AB13DD.f"
    if (discr) {

/*        Try the frequency w = pi. */

#line 1143 "AB13DD.f"
	pi = atan(1.) * 4.;
#line 1144 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 1144 "AB13DD.f"
	gamma = ab13dx_(dico, jobe, jobd, n, m, p, &pi, &dwork[ia], n, &dwork[
		ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[id], p, &iwork[1]
		, &dwork[iwrk], &i__1, &cwork[1], lcwork, &ierr, (ftnlen)1, (
		ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 1148 "AB13DD.f"
	i__1 = (integer) cwork[1].r;
#line 1148 "AB13DD.f"
	maxcwk = max(i__1,maxcwk);
/* Computing MAX */
#line 1149 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 1149 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);
#line 1150 "AB13DD.f"
	if (ierr >= 1 && ierr <= *n) {
#line 1151 "AB13DD.f"
	    gpeak[1] = 1.;
#line 1152 "AB13DD.f"
	    fpeak[1] = pi;
#line 1153 "AB13DD.f"
	    gpeak[2] = 0.;
#line 1154 "AB13DD.f"
	    fpeak[2] = 1.;
#line 1155 "AB13DD.f"
	    goto L340;
#line 1156 "AB13DD.f"
	} else if (ierr == *n + 1) {
#line 1157 "AB13DD.f"
	    *info = 3;
#line 1158 "AB13DD.f"
	    return 0;
#line 1159 "AB13DD.f"
	}

#line 1161 "AB13DD.f"
	if (gammal < gamma) {
#line 1162 "AB13DD.f"
	    gammal = gamma;
#line 1163 "AB13DD.f"
	    fpeak[1] = pi;
#line 1164 "AB13DD.f"
	    fpeak[2] = 1.;
#line 1165 "AB13DD.f"
	}

#line 1167 "AB13DD.f"
    } else {
#line 1168 "AB13DD.f"
	iwrk = ias;

/*        Restore D, if needed. */

#line 1172 "AB13DD.f"
	if (withd) {
#line 1172 "AB13DD.f"
	    dlacpy_("Full", p, m, &d__[d_offset], ldd, &dwork[id], p, (ftnlen)
		    4);
#line 1172 "AB13DD.f"
	}
#line 1174 "AB13DD.f"
    }

/*     Build the remaining set of frequencies. */
/*     Complex workspace:  need   (N+M)*(N+P)+2*MIN(P,M)+MAX(P,M)); */
/*                         prefer larger. */
/*     Real workspace:     need   LDW2, see above; */
/*                         prefer larger. */

#line 1182 "AB13DD.f"
    if (min(fpeaks,fpeaki) != 0.) {

/*        Compute also the norm at the given (finite) frequency. */

#line 1186 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 1186 "AB13DD.f"
	gamma = ab13dx_(dico, jobe, jobd, n, m, p, &fpeaks, &dwork[ia], n, &
		dwork[ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[id], p, &
		iwork[1], &dwork[iwrk], &i__1, &cwork[1], lcwork, &ierr, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 1190 "AB13DD.f"
	i__1 = (integer) cwork[1].r;
#line 1190 "AB13DD.f"
	maxcwk = max(i__1,maxcwk);
/* Computing MAX */
#line 1191 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 1191 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);
#line 1192 "AB13DD.f"
	if (discr) {
#line 1193 "AB13DD.f"
	    tm = (d__1 = atan2(sin(fpeaks), cos(fpeaks)), abs(d__1));
#line 1194 "AB13DD.f"
	} else {
#line 1195 "AB13DD.f"
	    tm = fpeaks;
#line 1196 "AB13DD.f"
	}
#line 1197 "AB13DD.f"
	if (ierr >= 1 && ierr <= *n) {
#line 1198 "AB13DD.f"
	    gpeak[1] = 1.;
#line 1199 "AB13DD.f"
	    fpeak[1] = tm;
#line 1200 "AB13DD.f"
	    gpeak[2] = 0.;
#line 1201 "AB13DD.f"
	    fpeak[2] = 1.;
#line 1202 "AB13DD.f"
	    goto L340;
#line 1203 "AB13DD.f"
	} else if (ierr == *n + 1) {
#line 1204 "AB13DD.f"
	    *info = 3;
#line 1205 "AB13DD.f"
	    return 0;
#line 1206 "AB13DD.f"
	}

#line 1208 "AB13DD.f"
	if (gammal < gamma) {
#line 1209 "AB13DD.f"
	    gammal = gamma;
#line 1210 "AB13DD.f"
	    fpeak[1] = tm;
#line 1211 "AB13DD.f"
	    fpeak[2] = 1.;
#line 1212 "AB13DD.f"
	}

#line 1214 "AB13DD.f"
    }

#line 1216 "AB13DD.f"
    i__1 = *n - 1;
#line 1216 "AB13DD.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1217 "AB13DD.f"
	if (dwork[ii + i__] >= 0. && dwork[im + i__] > 0.) {
#line 1218 "AB13DD.f"
	    if (dwork[im + i__] >= 1. || dwork[im + i__] < 1. && dwork[iar + 
		    i__] < safmax * dwork[im + i__]) {
#line 1220 "AB13DD.f"
		rat = dwork[iar + i__] / dwork[im + i__];
#line 1221 "AB13DD.f"
	    } else {
#line 1222 "AB13DD.f"
		rat = 1.;
#line 1223 "AB13DD.f"
	    }
/* Computing MAX */
/* Computing 2nd power */
#line 1224 "AB13DD.f"
	    d__3 = rat;
#line 1224 "AB13DD.f"
	    d__1 = .25, d__2 = 1. - d__3 * d__3 * 2.;
#line 1224 "AB13DD.f"
	    omega = dwork[im + i__] * sqrt((max(d__1,d__2)));

#line 1226 "AB13DD.f"
	    i__2 = *ldwork - iwrk + 1;
#line 1226 "AB13DD.f"
	    gamma = ab13dx_(dico, jobe, jobd, n, m, p, &omega, &dwork[ia], n, 
		    &dwork[ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[id], 
		    p, &iwork[1], &dwork[iwrk], &i__2, &cwork[1], lcwork, &
		    ierr, (ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 1231 "AB13DD.f"
	    i__2 = (integer) cwork[1].r;
#line 1231 "AB13DD.f"
	    maxcwk = max(i__2,maxcwk);
/* Computing MAX */
#line 1232 "AB13DD.f"
	    i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 1232 "AB13DD.f"
	    maxwrk = max(i__2,maxwrk);
#line 1233 "AB13DD.f"
	    if (discr) {
#line 1234 "AB13DD.f"
		tm = (d__1 = atan2(sin(omega), cos(omega)), abs(d__1));
#line 1235 "AB13DD.f"
	    } else {
#line 1236 "AB13DD.f"
		tm = omega;
#line 1237 "AB13DD.f"
	    }
#line 1238 "AB13DD.f"
	    if (ierr >= 1 && ierr <= *n) {
#line 1239 "AB13DD.f"
		gpeak[1] = 1.;
#line 1240 "AB13DD.f"
		fpeak[1] = tm;
#line 1241 "AB13DD.f"
		gpeak[2] = 0.;
#line 1242 "AB13DD.f"
		fpeak[2] = 1.;
#line 1243 "AB13DD.f"
		goto L340;
#line 1244 "AB13DD.f"
	    } else if (ierr == *n + 1) {
#line 1245 "AB13DD.f"
		*info = 3;
#line 1246 "AB13DD.f"
		return 0;
#line 1247 "AB13DD.f"
	    }

#line 1249 "AB13DD.f"
	    if (gammal < gamma) {
#line 1250 "AB13DD.f"
		gammal = gamma;
#line 1251 "AB13DD.f"
		fpeak[1] = tm;
#line 1252 "AB13DD.f"
		fpeak[2] = 1.;
#line 1253 "AB13DD.f"
	    }

#line 1255 "AB13DD.f"
	}
#line 1256 "AB13DD.f"
/* L110: */
#line 1256 "AB13DD.f"
    }

/*     Return if the lower bound is zero. */

#line 1260 "AB13DD.f"
    if (gammal == 0.) {
#line 1261 "AB13DD.f"
	gpeak[1] = 0.;
#line 1262 "AB13DD.f"
	fpeak[1] = 0.;
#line 1263 "AB13DD.f"
	gpeak[2] = 1.;
#line 1264 "AB13DD.f"
	fpeak[2] = 1.;
#line 1265 "AB13DD.f"
	goto L340;
#line 1266 "AB13DD.f"
    }

/*     Start the modified gamma iteration for the Bruinsma-Steinbuch */
/*     algorithm. */

#line 1271 "AB13DD.f"
    if (! discr) {
#line 1271 "AB13DD.f"
	rtol = toler * 100.;
#line 1271 "AB13DD.f"
    }
#line 1273 "AB13DD.f"
    iter = 0;

/*     WHILE ( Iteration may continue ) DO */

#line 1277 "AB13DD.f"
L120:

#line 1279 "AB13DD.f"
    ++iter;
#line 1280 "AB13DD.f"
    gamma = (*tol + 1.) * gammal;
#line 1281 "AB13DD.f"
    usepen = fulle || discr;
#line 1282 "AB13DD.f"
    if (! usepen && withd) {

/*           Check whether one can use an explicit Hamiltonian matrix: */
/*           compute */
/*           min(rcond(GAMMA**2*Im - S'*S), rcond(GAMMA**2*Ip - S*S')). */
/*           If P = M = 1, then GAMMA**2 - S(1)**2 is used instead. */

#line 1289 "AB13DD.f"
	if (*m != *p) {
/* Computing 2nd power */
#line 1290 "AB13DD.f"
	    d__1 = dwork[is] / gamma;
#line 1290 "AB13DD.f"
	    rcond = 1. - d__1 * d__1;
#line 1291 "AB13DD.f"
	} else if (minpm > 1) {
/* Computing 2nd power */
#line 1292 "AB13DD.f"
	    d__1 = gamma;
/* Computing 2nd power */
#line 1292 "AB13DD.f"
	    d__2 = dwork[is];
/* Computing 2nd power */
#line 1292 "AB13DD.f"
	    d__3 = gamma;
/* Computing 2nd power */
#line 1292 "AB13DD.f"
	    d__4 = dwork[is + *p - 1];
#line 1292 "AB13DD.f"
	    rcond = (d__1 * d__1 - d__2 * d__2) / (d__3 * d__3 - d__4 * d__4);
#line 1294 "AB13DD.f"
	} else {
/* Computing 2nd power */
#line 1295 "AB13DD.f"
	    d__1 = gamma;
/* Computing 2nd power */
#line 1295 "AB13DD.f"
	    d__2 = dwork[is];
#line 1295 "AB13DD.f"
	    rcond = d__1 * d__1 - d__2 * d__2;
#line 1296 "AB13DD.f"
	}

#line 1298 "AB13DD.f"
	usepen = rcond < rtol;
#line 1299 "AB13DD.f"
    }

#line 1301 "AB13DD.f"
    if (usepen) {

/*           Use the QZ algorithm on a pencil. */
/*           Additional workspace here:  need   6*N.   (from IR) */

#line 1306 "AB13DD.f"
	ii = ir + n2;
#line 1307 "AB13DD.f"
	ibt = ii + n2;
#line 1308 "AB13DD.f"
	ih12 = ibt + n2;
#line 1309 "AB13DD.f"
	im = ih12;

/*           Set up the needed parts of the Hamiltonian pencil (H,J), */

/*                  ( H11  H12 ) */
/*              H = (          ) , */
/*                  ( H21  H22 ) */

/*           with */

/*                 ( A  0  )            ( 0  B )            ( E  0  ) */
/*           H11 = (       ),     H12 = (      )/nB,  J11 = (       ), */
/*                 ( 0 -A' )            ( C' 0 )            ( 0  E' ) */

/*                 ( C  0  )            ( Ip  D/g ) */
/*           H21 = (       )*nB,  H22 = (         ), */
/*                 ( 0 -B' )            ( D'/g Im ) */

/*           if DICO = 'C', and */

/*                 ( A  0  )            ( B  0  )            ( E  0 ) */
/*           H11 = (       ),     H12 = (       )/nB,  J11 = (      ), */
/*                 ( 0  E' )            ( 0  C' )            ( 0  A') */

/*                 ( 0  0  )            ( Im  D'/g )         ( 0  B') */
/*           H21 = (       )*nB,  H22 = (          ),  J21 = (      )*nB, */
/*                 ( C  0  )            ( D/g  Ip  )         ( 0  0 ) */

/*           if DICO = 'D', where g = GAMMA, and nB = norm(B,1). */
/*           First build [H12; H22]. */

#line 1340 "AB13DD.f"
	temp[0] = 0.;
#line 1341 "AB13DD.f"
	ih = ih12;

#line 1343 "AB13DD.f"
	if (discr) {

#line 1345 "AB13DD.f"
	    i__1 = *m;
#line 1345 "AB13DD.f"
	    for (j = 1; j <= i__1; ++j) {

#line 1347 "AB13DD.f"
		i__2 = *n;
#line 1347 "AB13DD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1348 "AB13DD.f"
		    dwork[ih] = b[i__ + j * b_dim1] / bnorm;
#line 1349 "AB13DD.f"
		    ++ih;
#line 1350 "AB13DD.f"
/* L130: */
#line 1350 "AB13DD.f"
		}

#line 1352 "AB13DD.f"
		i__2 = *n + *m;
#line 1352 "AB13DD.f"
		dcopy_(&i__2, temp, &c__0, &dwork[ih], &c__1);
#line 1353 "AB13DD.f"
		dwork[ih + *n + j - 1] = 1.;
#line 1354 "AB13DD.f"
		ih = ih + *n + *m;

#line 1356 "AB13DD.f"
		i__2 = *p;
#line 1356 "AB13DD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1357 "AB13DD.f"
		    dwork[ih] = d__[i__ + j * d_dim1] / gamma;
#line 1358 "AB13DD.f"
		    ++ih;
#line 1359 "AB13DD.f"
/* L140: */
#line 1359 "AB13DD.f"
		}

#line 1361 "AB13DD.f"
/* L150: */
#line 1361 "AB13DD.f"
	    }

#line 1363 "AB13DD.f"
	    i__1 = *p;
#line 1363 "AB13DD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 1364 "AB13DD.f"
		dcopy_(n, temp, &c__0, &dwork[ih], &c__1);
#line 1365 "AB13DD.f"
		ih += *n;

#line 1367 "AB13DD.f"
		i__2 = *n;
#line 1367 "AB13DD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1368 "AB13DD.f"
		    dwork[ih] = c__[j + i__ * c_dim1] / bnorm;
#line 1369 "AB13DD.f"
		    ++ih;
#line 1370 "AB13DD.f"
/* L160: */
#line 1370 "AB13DD.f"
		}

#line 1372 "AB13DD.f"
		i__2 = *m;
#line 1372 "AB13DD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1373 "AB13DD.f"
		    dwork[ih] = d__[j + i__ * d_dim1] / gamma;
#line 1374 "AB13DD.f"
		    ++ih;
#line 1375 "AB13DD.f"
/* L170: */
#line 1375 "AB13DD.f"
		}

#line 1377 "AB13DD.f"
		dcopy_(p, temp, &c__0, &dwork[ih], &c__1);
#line 1378 "AB13DD.f"
		dwork[ih + j - 1] = 1.;
#line 1379 "AB13DD.f"
		ih += *p;
#line 1380 "AB13DD.f"
/* L180: */
#line 1380 "AB13DD.f"
	    }

#line 1382 "AB13DD.f"
	} else {

#line 1384 "AB13DD.f"
	    i__1 = *p;
#line 1384 "AB13DD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 1385 "AB13DD.f"
		dcopy_(n, temp, &c__0, &dwork[ih], &c__1);
#line 1386 "AB13DD.f"
		ih += *n;

#line 1388 "AB13DD.f"
		i__2 = *n;
#line 1388 "AB13DD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1389 "AB13DD.f"
		    dwork[ih] = c__[j + i__ * c_dim1] / bnorm;
#line 1390 "AB13DD.f"
		    ++ih;
#line 1391 "AB13DD.f"
/* L190: */
#line 1391 "AB13DD.f"
		}

#line 1393 "AB13DD.f"
		dcopy_(p, temp, &c__0, &dwork[ih], &c__1);
#line 1394 "AB13DD.f"
		dwork[ih + j - 1] = 1.;
#line 1395 "AB13DD.f"
		ih += *p;

#line 1397 "AB13DD.f"
		i__2 = *m;
#line 1397 "AB13DD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1398 "AB13DD.f"
		    dwork[ih] = d__[j + i__ * d_dim1] / gamma;
#line 1399 "AB13DD.f"
		    ++ih;
#line 1400 "AB13DD.f"
/* L200: */
#line 1400 "AB13DD.f"
		}

#line 1402 "AB13DD.f"
/* L210: */
#line 1402 "AB13DD.f"
	    }

#line 1404 "AB13DD.f"
	    i__1 = *m;
#line 1404 "AB13DD.f"
	    for (j = 1; j <= i__1; ++j) {

#line 1406 "AB13DD.f"
		i__2 = *n;
#line 1406 "AB13DD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1407 "AB13DD.f"
		    dwork[ih] = b[i__ + j * b_dim1] / bnorm;
#line 1408 "AB13DD.f"
		    ++ih;
#line 1409 "AB13DD.f"
/* L220: */
#line 1409 "AB13DD.f"
		}

#line 1411 "AB13DD.f"
		dcopy_(n, temp, &c__0, &dwork[ih], &c__1);
#line 1412 "AB13DD.f"
		ih += *n;

#line 1414 "AB13DD.f"
		i__2 = *p;
#line 1414 "AB13DD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 1415 "AB13DD.f"
		    dwork[ih] = d__[i__ + j * d_dim1] / gamma;
#line 1416 "AB13DD.f"
		    ++ih;
#line 1417 "AB13DD.f"
/* L230: */
#line 1417 "AB13DD.f"
		}

#line 1419 "AB13DD.f"
		dcopy_(m, temp, &c__0, &dwork[ih], &c__1);
#line 1420 "AB13DD.f"
		dwork[ih + j - 1] = 1.;
#line 1421 "AB13DD.f"
		ih += *m;
#line 1422 "AB13DD.f"
/* L240: */
#line 1422 "AB13DD.f"
	    }

#line 1424 "AB13DD.f"
	}

/*           Compute the QR factorization of [H12; H22]. */
/*           For large P and M, it could be more efficient to exploit the */
/*           structure of [H12; H22] and use the factored form of Q. */
/*           Additional workspace: need   (2*N+P+M)*(2*N+P+M)+2*(P+M); */
/*                                 prefer (2*N+P+M)*(2*N+P+M)+P+M+ */
/*                                                           (P+M)*NB. */

#line 1433 "AB13DD.f"
	itau = ih12 + n2pm * n2pm;
#line 1434 "AB13DD.f"
	iwrk = itau + pm;
#line 1435 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 1435 "AB13DD.f"
	dgeqrf_(&n2pm, &pm, &dwork[ih12], &n2pm, &dwork[itau], &dwork[iwrk], &
		i__1, &ierr);
/* Computing MAX */
#line 1437 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 1437 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);

/*           Apply part of the orthogonal transformation: */
/*           Q1 = Q(:,P+M+(1:2*N))' to the matrix [H11; H21/GAMMA]. */
/*           If DICO = 'C', apply Q(1:2*N,P+M+(1:2*N))' to the */
/*           matrix J11. */
/*           If DICO = 'D', apply Q1 to the matrix [J11; J21/GAMMA]. */
/*           H11, H21, J11, and J21 are not fully built. */
/*           First, build the (2*N+P+M)-by-(2*N+P+M) matrix Q. */
/*           Using Q will often provide better efficiency than the direct */
/*           use of the factored form of Q, especially when P+M < N. */
/*           Additional workspace: need   P+M+2*N+P+M; */
/*                                 prefer P+M+(2*N+P+M)*NB. */

#line 1451 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 1451 "AB13DD.f"
	dorgqr_(&n2pm, &n2pm, &pm, &dwork[ih12], &n2pm, &dwork[itau], &dwork[
		iwrk], &i__1, &ierr);
/* Computing MAX */
#line 1454 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 1454 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);

/*           Additional workspace: need   8*N*N. */

#line 1458 "AB13DD.f"
	ipa = itau;
#line 1459 "AB13DD.f"
	ipe = ipa + (nn << 2);
#line 1460 "AB13DD.f"
	iwrk = ipe + (nn << 2);
#line 1461 "AB13DD.f"
	dgemm_("Transpose", "No Transpose", &n2, n, n, &c_b20, &dwork[ih12 + 
		pm * n2pm], &n2pm, &a[a_offset], lda, &c_b21, &dwork[ipa], &
		n2, (ftnlen)9, (ftnlen)12);
#line 1464 "AB13DD.f"
	if (discr) {
#line 1465 "AB13DD.f"
	    d__1 = bnorm / gamma;
#line 1465 "AB13DD.f"
	    dgemm_("Transpose", "No Transpose", &n2, n, p, &d__1, &dwork[ih12 
		    + pm * n2pm + n2 + *m], &n2pm, &c__[c_offset], ldc, &
		    c_b20, &dwork[ipa], &n2, (ftnlen)9, (ftnlen)12);
#line 1468 "AB13DD.f"
	    if (fulle) {
#line 1469 "AB13DD.f"
		dgemm_("Transpose", "Transpose", &n2, n, n, &c_b20, &dwork[
			ih12 + pm * n2pm + *n], &n2pm, &e[e_offset], lde, &
			c_b21, &dwork[ipa + (nn << 1)], &n2, (ftnlen)9, (
			ftnlen)9);
#line 1472 "AB13DD.f"
	    } else {
#line 1473 "AB13DD.f"
		ma02ad_("Full", n, &n2, &dwork[ih12 + pm * n2pm + *n], &n2pm, 
			&dwork[ipa + (nn << 1)], &n2, (ftnlen)4);
#line 1475 "AB13DD.f"
		ny = *n;
#line 1476 "AB13DD.f"
	    }
#line 1477 "AB13DD.f"
	} else {
#line 1478 "AB13DD.f"
	    d__1 = bnorm / gamma;
#line 1478 "AB13DD.f"
	    dgemm_("Transpose", "No Transpose", &n2, n, p, &d__1, &dwork[ih12 
		    + pm * n2pm + n2], &n2pm, &c__[c_offset], ldc, &c_b20, &
		    dwork[ipa], &n2, (ftnlen)9, (ftnlen)12);
#line 1481 "AB13DD.f"
	    dgemm_("Transpose", "Transpose", &n2, n, n, &c_b163, &dwork[ih12 
		    + pm * n2pm + *n], &n2pm, &a[a_offset], lda, &c_b21, &
		    dwork[ipa + (nn << 1)], &n2, (ftnlen)9, (ftnlen)9);
#line 1484 "AB13DD.f"
	    d__1 = -bnorm / gamma;
#line 1484 "AB13DD.f"
	    dgemm_("Transpose", "Transpose", &n2, n, m, &d__1, &dwork[ih12 + 
		    pm * n2pm + n2 + *p], &n2pm, &b[b_offset], ldb, &c_b20, &
		    dwork[ipa + (nn << 1)], &n2, (ftnlen)9, (ftnlen)9);
#line 1487 "AB13DD.f"
	    ny = n2;
#line 1488 "AB13DD.f"
	}

#line 1490 "AB13DD.f"
	if (fulle) {
#line 1491 "AB13DD.f"
	    dgemm_("Transpose", "No Transpose", &n2, n, n, &c_b20, &dwork[
		    ih12 + pm * n2pm], &n2pm, &e[e_offset], lde, &c_b21, &
		    dwork[ipe], &n2, (ftnlen)9, (ftnlen)12);
#line 1494 "AB13DD.f"
	} else {
#line 1495 "AB13DD.f"
	    ma02ad_("Full", &ny, &n2, &dwork[ih12 + pm * n2pm], &n2pm, &dwork[
		    ipe], &n2, (ftnlen)4);
#line 1497 "AB13DD.f"
	}
#line 1498 "AB13DD.f"
	if (discr) {
#line 1499 "AB13DD.f"
	    dgemm_("Transpose", "Transpose", &n2, n, n, &c_b20, &dwork[ih12 + 
		    pm * n2pm + *n], &n2pm, &a[a_offset], lda, &c_b21, &dwork[
		    ipe + (nn << 1)], &n2, (ftnlen)9, (ftnlen)9);
#line 1502 "AB13DD.f"
	    d__1 = bnorm / gamma;
#line 1502 "AB13DD.f"
	    dgemm_("Transpose", "Transpose", &n2, n, m, &d__1, &dwork[ih12 + 
		    pm * n2pm + n2], &n2pm, &b[b_offset], ldb, &c_b20, &dwork[
		    ipe + (nn << 1)], &n2, (ftnlen)9, (ftnlen)9);
#line 1505 "AB13DD.f"
	} else {
#line 1506 "AB13DD.f"
	    if (fulle) {
#line 1506 "AB13DD.f"
		dgemm_("Transpose", "Transpose", &n2, n, n, &c_b20, &dwork[
			ih12 + pm * n2pm + *n], &n2pm, &e[e_offset], lde, &
			c_b21, &dwork[ipe + (nn << 1)], &n2, (ftnlen)9, (
			ftnlen)9);
#line 1506 "AB13DD.f"
	    }
#line 1510 "AB13DD.f"
	}

/*           Compute the eigenvalues of the Hamiltonian pencil. */
/*           Additional workspace: need   16*N; */
/*                                 prefer larger. */

#line 1516 "AB13DD.f"
	i__1 = *ldwork - iwrk + 1;
#line 1516 "AB13DD.f"
	dggev_("No Vectors", "No Vectors", &n2, &dwork[ipa], &n2, &dwork[ipe],
		 &n2, &dwork[ir], &dwork[ii], &dwork[ibt], &dwork[1], &n2, &
		dwork[1], &n2, &dwork[iwrk], &i__1, &ierr, (ftnlen)10, (
		ftnlen)10);
#line 1520 "AB13DD.f"
	if (ierr > 0) {
#line 1521 "AB13DD.f"
	    *info = 2;
#line 1522 "AB13DD.f"
	    return 0;
#line 1523 "AB13DD.f"
	}
/* Computing MAX */
#line 1524 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk - 1;
#line 1524 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);

#line 1526 "AB13DD.f"
    } else if (! withd) {

/*           Standard continuous-time case with D = 0. */
/*           Form the needed part of the Hamiltonian matrix explicitly: */
/*              H = H11 - H12*inv(H22)*H21/g. */
/*           Additional workspace: need   2*N*N+N.   (from IBT) */

#line 1533 "AB13DD.f"
	ih = ibt;
#line 1534 "AB13DD.f"
	ih12 = ih + nn;
#line 1535 "AB13DD.f"
	isl = ih12 + nn + *n;
#line 1536 "AB13DD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ih], n, (ftnlen)4);

/*           Compute triangles of -C'*C/GAMMA and B*B'/GAMMA. */

#line 1540 "AB13DD.f"
	d__1 = -1. / gamma;
#line 1540 "AB13DD.f"
	dsyrk_("Lower", "Transpose", n, p, &d__1, &c__[c_offset], ldc, &c_b21,
		 &dwork[ih12], n, (ftnlen)5, (ftnlen)9);
#line 1542 "AB13DD.f"
	d__1 = 1. / gamma;
#line 1542 "AB13DD.f"
	dsyrk_("Upper", "No Transpose", n, m, &d__1, &b[b_offset], ldb, &
		c_b21, &dwork[ih12 + *n], n, (ftnlen)5, (ftnlen)12);

#line 1545 "AB13DD.f"
    } else {

/*           Standard continuous-time case with D <> 0 and the SVD of D */
/*           can be used. Compute explicitly the needed part of the */
/*           Hamiltonian matrix: */

/*               (A+B1*S'*inv(g^2*Ip-S*S')*C1' g*B1*inv(g^2*Im-S'*S)*B1') */
/*           H = (                                                      ) */
/*               (  -g*C1*inv(g^2*Ip-S*S')*C1'            -H11'         ) */

/*           where g = GAMMA, B1 = B*V, C1 = C'*U, and H11 is the first */
/*           block of H. */
/*           Primary additional workspace: need   2*N*N+N   (from IBT) */
/*           (for building the relevant part of the Hamiltonian matrix). */

/*           Compute C1*sqrt(inv(g^2*Ip-S*S')) . */
/*           Additional workspace: need   MAX(M,P)+N*P. */

#line 1563 "AB13DD.f"
	ih = ibt;
#line 1564 "AB13DD.f"
	ih12 = ih + nn;
#line 1565 "AB13DD.f"
	isl = ih12 + nn + *n;

#line 1567 "AB13DD.f"
	i__1 = minpm - 1;
#line 1567 "AB13DD.f"
	for (i__ = 0; i__ <= i__1; ++i__) {
/* Computing 2nd power */
#line 1568 "AB13DD.f"
	    d__1 = gamma;
/* Computing 2nd power */
#line 1568 "AB13DD.f"
	    d__2 = dwork[is + i__];
#line 1568 "AB13DD.f"
	    dwork[isl + i__] = 1. / sqrt(d__1 * d__1 - d__2 * d__2);
#line 1569 "AB13DD.f"
/* L250: */
#line 1569 "AB13DD.f"
	}

#line 1571 "AB13DD.f"
	if (*m < *p) {
#line 1572 "AB13DD.f"
	    dwork[isl + *m] = 1. / gamma;
#line 1573 "AB13DD.f"
	    i__1 = *p - *m - 1;
#line 1573 "AB13DD.f"
	    dcopy_(&i__1, &dwork[isl + *m], &c__0, &dwork[isl + *m + 1], &
		    c__1);
#line 1575 "AB13DD.f"
	}
#line 1576 "AB13DD.f"
	isc = isl + max(*m,*p);
#line 1577 "AB13DD.f"
	dlacpy_("Full", n, p, &dwork[icu], n, &dwork[isc], n, (ftnlen)4);
#line 1579 "AB13DD.f"
	mb01sd_("Column", n, p, &dwork[isc], n, &dwork[1], &dwork[isl], (
		ftnlen)6);

/*           Compute B1*S' . */
/*           Additional workspace: need   N*M. */

#line 1585 "AB13DD.f"
	isb = isc + *p * *n;
#line 1586 "AB13DD.f"
	dlacpy_("Full", n, m, &dwork[ibv], n, &dwork[isb], n, (ftnlen)4);
#line 1588 "AB13DD.f"
	mb01sd_("Column", n, &minpm, &dwork[isb], n, &dwork[1], &dwork[is], (
		ftnlen)6);

/*           Compute B1*S'*sqrt(inv(g^2*Ip-S*S')) . */

#line 1593 "AB13DD.f"
	mb01sd_("Column", n, &minpm, &dwork[isb], n, &dwork[1], &dwork[isl], (
		ftnlen)6);

/*           Compute H11 . */

#line 1598 "AB13DD.f"
	dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ih], n, (ftnlen)4);
#line 1599 "AB13DD.f"
	dgemm_("No Transpose", "Transpose", n, n, &minpm, &c_b20, &dwork[isb],
		 n, &dwork[isc], n, &c_b20, &dwork[ih], n, (ftnlen)12, (
		ftnlen)9);

/*           Compute B1*sqrt(inv(g^2*Im-S'*S)) . */

#line 1605 "AB13DD.f"
	if (*p < *m) {
#line 1606 "AB13DD.f"
	    dwork[isl + *p] = 1. / gamma;
#line 1607 "AB13DD.f"
	    i__1 = *m - *p - 1;
#line 1607 "AB13DD.f"
	    dcopy_(&i__1, &dwork[isl + *p], &c__0, &dwork[isl + *p + 1], &
		    c__1);
#line 1609 "AB13DD.f"
	}
#line 1610 "AB13DD.f"
	dlacpy_("Full", n, m, &dwork[ibv], n, &dwork[isb], n, (ftnlen)4);
#line 1612 "AB13DD.f"
	mb01sd_("Column", n, m, &dwork[isb], n, &dwork[1], &dwork[isl], (
		ftnlen)6);

/*           Compute the lower triangle of H21 and the upper triangle */
/*           of H12. */

#line 1618 "AB13DD.f"
	d__1 = -gamma;
#line 1618 "AB13DD.f"
	dsyrk_("Lower", "No Transpose", n, p, &d__1, &dwork[isc], n, &c_b21, &
		dwork[ih12], n, (ftnlen)5, (ftnlen)12);
#line 1620 "AB13DD.f"
	dsyrk_("Upper", "No Transpose", n, m, &gamma, &dwork[isb], n, &c_b21, 
		&dwork[ih12 + *n], n, (ftnlen)5, (ftnlen)12);
#line 1622 "AB13DD.f"
    }

#line 1624 "AB13DD.f"
    if (! usepen) {

/*           Compute the eigenvalues of the Hamiltonian matrix by the */
/*           symplectic URV and the periodic Schur decompositions. */
/*           Additional workspace: need   (2*N+8)*N; */
/*                                 prefer larger. */

#line 1631 "AB13DD.f"
	iwrk = isl + nn;
#line 1632 "AB13DD.f"
	i__1 = *ldwork - iwrk - *n + 1;
#line 1632 "AB13DD.f"
	mb03xd_("Both", "Eigenvalues", "No vectors", "No vectors", n, &dwork[
		ih], n, &dwork[ih12], n, &dwork[isl], n, temp, &c__1, temp, &
		c__1, temp, &c__1, temp, &c__1, &dwork[ir], &dwork[ii], &ilo, 
		&dwork[iwrk], &dwork[iwrk + *n], &i__1, &ierr, (ftnlen)4, (
		ftnlen)11, (ftnlen)10, (ftnlen)10);
#line 1638 "AB13DD.f"
	if (ierr > 0) {
#line 1639 "AB13DD.f"
	    *info = 2;
#line 1640 "AB13DD.f"
	    return 0;
#line 1641 "AB13DD.f"
	}
/* Computing MAX */
#line 1642 "AB13DD.f"
	i__1 = (integer) dwork[iwrk] + iwrk + *n - 1;
#line 1642 "AB13DD.f"
	maxwrk = max(i__1,maxwrk);
#line 1643 "AB13DD.f"
    }

/*        Detect eigenvalues on the boundary of the stability domain, */
/*        if any. The test is based on a round-off level of eps*rho(H) */
/*        (after balancing) resulting in worst-case perturbations of */
/*        order sqrt(eps*rho(H)), for continuous-time systems, on the */
/*        real part of poles of multiplicity two (typical as GAMMA */
/*        approaches the infinity norm). Similarly, in the discrete-time */
/*        case. Above, rho(H) is the maximum modulus of eigenvalues */
/*        (continuous-time case). */

/*        Compute maximum eigenvalue modulus and check the absolute real */
/*        parts (if DICO = 'C'), or moduli (if DICO = 'D'). */

#line 1657 "AB13DD.f"
    wmax = 0.;

#line 1659 "AB13DD.f"
    if (usepen) {

/*           Additional workspace: need   2*N, if DICO = 'D';   (from IM) */
/*                                        0,   if DICO = 'C'. */

#line 1664 "AB13DD.f"
	i__1 = n2 - 1;
#line 1664 "AB13DD.f"
	for (i__ = 0; i__ <= i__1; ++i__) {
#line 1665 "AB13DD.f"
	    tm = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
#line 1666 "AB13DD.f"
	    if (dwork[ibt + i__] >= 1. || dwork[ibt + i__] < 1. && tm < 
		    safmax * dwork[ibt + i__]) {
#line 1669 "AB13DD.f"
		tm /= dwork[ibt + i__];
#line 1670 "AB13DD.f"
	    } else {

/*                 The pencil has too large eigenvalues. SAFMAX is used. */

#line 1674 "AB13DD.f"
		tm = safmax;
#line 1675 "AB13DD.f"
	    }
#line 1676 "AB13DD.f"
	    wmax = max(wmax,tm);
#line 1677 "AB13DD.f"
	    if (discr) {
#line 1677 "AB13DD.f"
		dwork[im + i__] = tm;
#line 1677 "AB13DD.f"
	    }
#line 1679 "AB13DD.f"
/* L260: */
#line 1679 "AB13DD.f"
	}

#line 1681 "AB13DD.f"
    } else {

#line 1683 "AB13DD.f"
	i__1 = *n - 1;
#line 1683 "AB13DD.f"
	for (i__ = 0; i__ <= i__1; ++i__) {
#line 1684 "AB13DD.f"
	    tm = dlapy2_(&dwork[ir + i__], &dwork[ii + i__]);
#line 1685 "AB13DD.f"
	    wmax = max(wmax,tm);
#line 1686 "AB13DD.f"
/* L270: */
#line 1686 "AB13DD.f"
	}

#line 1688 "AB13DD.f"
    }

#line 1690 "AB13DD.f"
    nei = 0;

#line 1692 "AB13DD.f"
    if (usepen) {

#line 1694 "AB13DD.f"
	i__1 = n2 - 1;
#line 1694 "AB13DD.f"
	for (i__ = 0; i__ <= i__1; ++i__) {
#line 1695 "AB13DD.f"
	    if (discr) {
#line 1696 "AB13DD.f"
		tm = (d__1 = 1. - dwork[im + i__], abs(d__1));
#line 1697 "AB13DD.f"
	    } else {
#line 1698 "AB13DD.f"
		tm = (d__1 = dwork[ir + i__], abs(d__1));
#line 1699 "AB13DD.f"
		if (dwork[ibt + i__] >= 1. || dwork[ibt + i__] < 1. && tm < 
			safmax * dwork[ibt + i__]) {
#line 1702 "AB13DD.f"
		    tm /= dwork[ibt + i__];
#line 1703 "AB13DD.f"
		} else {

/*                    The pencil has too large eigenvalues. */
/*                    SAFMAX is used. */

#line 1708 "AB13DD.f"
		    tm = safmax;
#line 1709 "AB13DD.f"
		}
#line 1710 "AB13DD.f"
	    }
#line 1711 "AB13DD.f"
	    if (tm <= toler * sqrt(wmax + 100.)) {
#line 1712 "AB13DD.f"
		dwork[ir + nei] = dwork[ir + i__] / dwork[ibt + i__];
#line 1713 "AB13DD.f"
		dwork[ii + nei] = dwork[ii + i__] / dwork[ibt + i__];
#line 1714 "AB13DD.f"
		++nei;
#line 1715 "AB13DD.f"
	    }
#line 1716 "AB13DD.f"
/* L280: */
#line 1716 "AB13DD.f"
	}

#line 1718 "AB13DD.f"
    } else {

#line 1720 "AB13DD.f"
	i__1 = *n - 1;
#line 1720 "AB13DD.f"
	for (i__ = 0; i__ <= i__1; ++i__) {
#line 1721 "AB13DD.f"
	    tm = (d__1 = dwork[ir + i__], abs(d__1));
#line 1722 "AB13DD.f"
	    if (tm <= toler * sqrt(wmax + 100.)) {
#line 1723 "AB13DD.f"
		dwork[ir + nei] = dwork[ir + i__];
#line 1724 "AB13DD.f"
		dwork[ii + nei] = dwork[ii + i__];
#line 1725 "AB13DD.f"
		++nei;
#line 1726 "AB13DD.f"
	    }
#line 1727 "AB13DD.f"
/* L290: */
#line 1727 "AB13DD.f"
	}

#line 1729 "AB13DD.f"
    }

#line 1731 "AB13DD.f"
    if (nei == 0) {

/*           There is no eigenvalue on the boundary of the stability */
/*           domain for G = ( ONE + TOL )*GAMMAL. The norm was found. */

#line 1736 "AB13DD.f"
	gpeak[1] = gammal;
#line 1737 "AB13DD.f"
	gpeak[2] = 1.;
#line 1738 "AB13DD.f"
	goto L340;
#line 1739 "AB13DD.f"
    }

/*        Compute the frequencies where the gain G is attained and */
/*        generate new test frequencies. */

#line 1744 "AB13DD.f"
    nws = 0;

#line 1746 "AB13DD.f"
    if (discr) {

#line 1748 "AB13DD.f"
	i__1 = nei - 1;
#line 1748 "AB13DD.f"
	for (i__ = 0; i__ <= i__1; ++i__) {
#line 1749 "AB13DD.f"
	    tm = atan2(dwork[ii + i__], dwork[ir + i__]);
#line 1750 "AB13DD.f"
	    dwork[ir + i__] = max(eps,tm);
#line 1751 "AB13DD.f"
	    ++nws;
#line 1752 "AB13DD.f"
/* L300: */
#line 1752 "AB13DD.f"
	}

#line 1754 "AB13DD.f"
    } else {

#line 1756 "AB13DD.f"
	j = 0;

#line 1758 "AB13DD.f"
	i__1 = nei - 1;
#line 1758 "AB13DD.f"
	for (i__ = 0; i__ <= i__1; ++i__) {
#line 1759 "AB13DD.f"
	    if (dwork[ii + i__] > eps) {
#line 1760 "AB13DD.f"
		dwork[ir + nws] = dwork[ii + i__];
#line 1761 "AB13DD.f"
		++nws;
#line 1762 "AB13DD.f"
	    } else if (dwork[ii + i__] == eps) {
#line 1763 "AB13DD.f"
		++j;
#line 1764 "AB13DD.f"
		if (j == 1) {
#line 1765 "AB13DD.f"
		    dwork[ir + nws] = eps;
#line 1766 "AB13DD.f"
		    ++nws;
#line 1767 "AB13DD.f"
		}
#line 1768 "AB13DD.f"
	    }
#line 1769 "AB13DD.f"
/* L310: */
#line 1769 "AB13DD.f"
	}

#line 1771 "AB13DD.f"
    }

#line 1773 "AB13DD.f"
    dlasrt_("Increasing", &nws, &dwork[ir], &ierr, (ftnlen)10);
#line 1774 "AB13DD.f"
    lw = 1;

#line 1776 "AB13DD.f"
    i__1 = nws - 1;
#line 1776 "AB13DD.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1777 "AB13DD.f"
	if (dwork[ir + lw - 1] != dwork[ir + i__]) {
#line 1778 "AB13DD.f"
	    dwork[ir + lw] = dwork[ir + i__];
#line 1779 "AB13DD.f"
	    ++lw;
#line 1780 "AB13DD.f"
	}
#line 1781 "AB13DD.f"
/* L320: */
#line 1781 "AB13DD.f"
    }

#line 1783 "AB13DD.f"
    if (lw == 1) {
#line 1784 "AB13DD.f"
	if (iter == 1 && nws >= 1) {

/*              Duplicate the frequency trying to force iteration. */

#line 1788 "AB13DD.f"
	    dwork[ir + 1] = dwork[ir];
#line 1789 "AB13DD.f"
	    ++lw;
#line 1790 "AB13DD.f"
	} else {

/*              The norm was found. */

#line 1794 "AB13DD.f"
	    gpeak[1] = gammal;
#line 1795 "AB13DD.f"
	    gpeak[2] = 1.;
#line 1796 "AB13DD.f"
	    goto L340;
#line 1797 "AB13DD.f"
	}
#line 1798 "AB13DD.f"
    }

/*        Form the vector of mid-points and compute the gain at new test */
/*        frequencies. Save the current lower bound. */

#line 1803 "AB13DD.f"
    iwrk = ir + lw;
#line 1804 "AB13DD.f"
    gammas = gammal;

#line 1806 "AB13DD.f"
    i__1 = lw - 2;
#line 1806 "AB13DD.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 1807 "AB13DD.f"
	if (discr) {
#line 1808 "AB13DD.f"
	    omega = (dwork[ir + i__] + dwork[ir + i__ + 1]) / 2.;
#line 1809 "AB13DD.f"
	} else {
#line 1810 "AB13DD.f"
	    omega = sqrt(dwork[ir + i__] * dwork[ir + i__ + 1]);
#line 1811 "AB13DD.f"
	}

/*           Additional workspace:  need   LDW2, see above; */
/*                                  prefer larger. */

#line 1816 "AB13DD.f"
	i__2 = *ldwork - iwrk + 1;
#line 1816 "AB13DD.f"
	gamma = ab13dx_(dico, jobe, jobd, n, m, p, &omega, &dwork[ia], n, &
		dwork[ie], n, &dwork[ib], n, &dwork[ic], p, &dwork[id], p, &
		iwork[1], &dwork[iwrk], &i__2, &cwork[1], lcwork, &ierr, (
		ftnlen)1, (ftnlen)1, (ftnlen)1);
/* Computing MAX */
#line 1821 "AB13DD.f"
	i__2 = (integer) cwork[1].r;
#line 1821 "AB13DD.f"
	maxcwk = max(i__2,maxcwk);
/* Computing MAX */
#line 1822 "AB13DD.f"
	i__2 = (integer) dwork[iwrk] + iwrk - 1;
#line 1822 "AB13DD.f"
	maxwrk = max(i__2,maxwrk);
#line 1823 "AB13DD.f"
	if (discr) {
#line 1824 "AB13DD.f"
	    tm = (d__1 = atan2(sin(omega), cos(omega)), abs(d__1));
#line 1825 "AB13DD.f"
	} else {
#line 1826 "AB13DD.f"
	    tm = omega;
#line 1827 "AB13DD.f"
	}
#line 1828 "AB13DD.f"
	if (ierr >= 1 && ierr <= *n) {
#line 1829 "AB13DD.f"
	    gpeak[1] = 1.;
#line 1830 "AB13DD.f"
	    fpeak[1] = tm;
#line 1831 "AB13DD.f"
	    gpeak[2] = 0.;
#line 1832 "AB13DD.f"
	    fpeak[2] = 1.;
#line 1833 "AB13DD.f"
	    goto L340;
#line 1834 "AB13DD.f"
	} else if (ierr == *n + 1) {
#line 1835 "AB13DD.f"
	    *info = 3;
#line 1836 "AB13DD.f"
	    return 0;
#line 1837 "AB13DD.f"
	}

#line 1839 "AB13DD.f"
	if (gammal < gamma) {
#line 1840 "AB13DD.f"
	    gammal = gamma;
#line 1841 "AB13DD.f"
	    fpeak[1] = tm;
#line 1842 "AB13DD.f"
	    fpeak[2] = 1.;
#line 1843 "AB13DD.f"
	}
#line 1844 "AB13DD.f"
/* L330: */
#line 1844 "AB13DD.f"
    }

/*        If the lower bound has not been improved, return. (This is a */
/*        safeguard against undetected modes of Hamiltonian matrix on the */
/*        boundary of the stability domain.) */

#line 1850 "AB13DD.f"
    if (gammal < gammas * (*tol / 10. + 1.)) {
#line 1851 "AB13DD.f"
	gpeak[1] = gammal;
#line 1852 "AB13DD.f"
	gpeak[2] = 1.;
#line 1853 "AB13DD.f"
	goto L340;
#line 1854 "AB13DD.f"
    }

/*     END WHILE */

#line 1858 "AB13DD.f"
    if (iter <= 30) {
#line 1859 "AB13DD.f"
	goto L120;
#line 1860 "AB13DD.f"
    } else {
#line 1861 "AB13DD.f"
	*info = 4;
#line 1862 "AB13DD.f"
	return 0;
#line 1863 "AB13DD.f"
    }

#line 1865 "AB13DD.f"
L340:
#line 1866 "AB13DD.f"
    dwork[1] = (doublereal) maxwrk;
#line 1867 "AB13DD.f"
    cwork[1].r = (doublereal) maxcwk, cwork[1].i = 0.;
#line 1868 "AB13DD.f"
    return 0;
/* *** Last line of AB13DD *** */
} /* ab13dd_ */

