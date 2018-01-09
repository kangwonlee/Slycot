#line 1 "IB01BD.f"
/* IB01BD.f -- translated by f2c (version 20100827).
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

#line 1 "IB01BD.f"
/* Table of constant values */

static doublereal c_b21 = 0.;
static integer c__0 = 0;

/* Subroutine */ int ib01bd_(char *meth, char *job, char *jobck, integer *
	nobr, integer *n, integer *m, integer *l, integer *nsmpl, doublereal *
	r__, integer *ldr, doublereal *a, integer *lda, doublereal *c__, 
	integer *ldc, doublereal *b, integer *ldb, doublereal *d__, integer *
	ldd, doublereal *q, integer *ldq, doublereal *ry, integer *ldry, 
	doublereal *s, integer *lds, doublereal *k, integer *ldk, doublereal *
	tol, integer *iwork, doublereal *dwork, integer *ldwork, logical *
	bwork, integer *iwarn, integer *info, ftnlen meth_len, ftnlen job_len,
	 ftnlen jobck_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, k_dim1, k_offset, q_dim1, q_offset, r_dim1, r_offset, 
	    ry_dim1, ry_offset, s_dim1, s_offset, i__1, i__2, i__3, i__4, 
	    i__5;

    /* Local variables */
    static integer i__, n2, ia, ic, id, ig, ik, io, ll, iq, ir, is, it, nl, 
	    iv, nn, ix, nr, iaw;
    static doublereal sep;
    static integer iwi, npl, iwr;
    static doublereal rcnd[8], ferr;
    static integer ierr;
    static logical n4sid;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ib01pd_(char *, char *, char *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);
    static char jobbd[1];
    static integer ifact;
    extern /* Subroutine */ int sb02nd_(char *, char *, char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, ftnlen, ftnlen, ftnlen, ftnlen), sb02rd_(
	    char *, char *, char *, char *, char *, char *, char *, char *, 
	    char *, integer *, doublereal *, integer *, doublereal *, integer 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, integer *, logical *, integer *, ftnlen,
	     ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static char jobcv[1];
    static doublereal rcond;
    extern /* Subroutine */ int sb02mt_(char *, char *, char *, char *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    integer *, integer *, doublereal *, integer *, integer *, ftnlen, 
	    ftnlen, ftnlen, ftnlen);
    static integer lnobr, mnobr;
    static logical withb, withc;
    static integer ldunn;
    static logical withd, moesp, withk;
    static integer jwork;
    static doublereal rnorm;
    static logical combin;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlaset_(char *, integer *, integer *, doublereal *, doublereal *, 
	    doublereal *, integer *, ftnlen), xerbla_(char *, integer *, 
	    ftnlen);
    static integer oufact[2];
    static char jobcov[1];
    static doublereal rcondr;
    static logical withal;
    static integer lmnobr, mnobrn, iwarnl;
    static logical withco;
    static integer lmmnol, minwrk, maxwrk;


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

/*     To estimate the system matrices A, C, B, and D, the noise */
/*     covariance matrices Q, Ry, and S, and the Kalman gain matrix K */
/*     of a linear time-invariant state space model, using the */
/*     processed triangular factor R of the concatenated block Hankel */
/*     matrices, provided by SLICOT Library routine IB01AD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     METH    CHARACTER*1 */
/*             Specifies the subspace identification method to be used, */
/*             as follows: */
/*             = 'M':  MOESP  algorithm with past inputs and outputs; */
/*             = 'N':  N4SID  algorithm; */
/*             = 'C':  combined method:  MOESP  algorithm for finding the */
/*                     matrices A and C, and  N4SID  algorithm for */
/*                     finding the matrices B and D. */

/*     JOB     CHARACTER*1 */
/*             Specifies which matrices should be computed, as follows: */
/*             = 'A':  compute all system matrices, A, B, C, and D; */
/*             = 'C':  compute the matrices A and C only; */
/*             = 'B':  compute the matrix B only; */
/*             = 'D':  compute the matrices B and D only. */

/*     JOBCK   CHARACTER*1 */
/*             Specifies whether or not the covariance matrices and the */
/*             Kalman gain matrix are to be computed, as follows: */
/*             = 'C':  the covariance matrices only should be computed; */
/*             = 'K':  the covariance matrices and the Kalman gain */
/*                     matrix should be computed; */
/*             = 'N':  the covariance matrices and the Kalman gain matrix */
/*                     should not be computed. */

/*     Input/Output Parameters */

/*     NOBR    (input) INTEGER */
/*             The number of block rows,  s,  in the input and output */
/*             Hankel matrices processed by other routines.  NOBR > 1. */

/*     N       (input) INTEGER */
/*             The order of the system.  NOBR > N > 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     NSMPL   (input) INTEGER */
/*             If  JOBCK = 'C' or 'K',  the total number of samples used */
/*             for calculating the covariance matrices. */
/*             NSMPL >= 2*(M+L)*NOBR. */
/*             This parameter is not meaningful if  JOBCK = 'N'. */

/*     R       (input/workspace) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             On entry, the leading  2*(M+L)*NOBR-by-2*(M+L)*NOBR  part */
/*             of this array must contain the relevant data for the MOESP */
/*             or N4SID algorithms, as constructed by SLICOT Library */
/*             routine IB01AD. Let  R_ij,  i,j = 1:4,  be the */
/*             ij submatrix of  R  (denoted  S  in IB01AD),  partitioned */
/*             by  M*NOBR,  L*NOBR,  M*NOBR,  and  L*NOBR  rows and */
/*             columns. The submatrix  R_22  contains the matrix of left */
/*             singular vectors used. Also needed, for  METH = 'N'  or */
/*             JOBCK <> 'N',  are the submatrices  R_11,  R_14 : R_44, */
/*             and, for  METH = 'M' or 'C'  and  JOB <> 'C', the */
/*             submatrices  R_31  and  R_12,  containing the processed */
/*             matrices  R_1c  and  R_2c,  respectively, as returned by */
/*             SLICOT Library routine IB01AD. */
/*             Moreover, if  METH = 'N'  and  JOB = 'A' or 'C',  the */
/*             block-row  R_41 : R_43  must contain the transpose of the */
/*             block-column  R_14 : R_34  as returned by SLICOT Library */
/*             routine IB01AD. */
/*             The remaining part of  R  is used as workspace. */
/*             On exit, part of this array is overwritten. Specifically, */
/*             if  METH = 'M',  R_22  and  R_31  are overwritten if */
/*                 JOB = 'B' or 'D',  and  R_12,  R_22,  R_14 : R_34, */
/*                 and possibly  R_11  are overwritten if  JOBCK <> 'N'; */
/*             if  METH = 'N',  all needed submatrices are overwritten. */
/*             The details of the contents of  R  need not be known if */
/*             this routine is called once just after calling the SLICOT */
/*             Library routine IB01AD. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= 2*(M+L)*NOBR. */

/*     A       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDA,N) */
/*             On entry, if  METH = 'N' or 'C'  and  JOB = 'B' or 'D', */
/*             the leading N-by-N part of this array must contain the */
/*             system state matrix. */
/*             If  METH = 'M'  or  (METH = 'N' or 'C'  and JOB = 'A' */
/*             or 'C'),  this array need not be set on input. */
/*             On exit, if  JOB = 'A' or 'C'  and  INFO = 0,  the */
/*             leading N-by-N part of this array contains the system */
/*             state matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= N,  if  JOB = 'A' or 'C',  or  METH = 'N' or 'C' */
/*                            and  JOB = 'B' or 'D'; */
/*             LDA >= 1,  otherwise. */

/*     C       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDC,N) */
/*             On entry, if  METH = 'N' or 'C'  and  JOB = 'B' or 'D', */
/*             the leading L-by-N part of this array must contain the */
/*             system output matrix. */
/*             If  METH = 'M'  or  (METH = 'N' or 'C'  and JOB = 'A' */
/*             or 'C'),  this array need not be set on input. */
/*             On exit, if  JOB = 'A' or 'C'  and  INFO = 0,  or */
/*             INFO = 3  (or  INFO >= 0,  for  METH = 'M'),  the leading */
/*             L-by-N part of this array contains the system output */
/*             matrix. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= L,  if  JOB = 'A' or 'C',  or  METH = 'N' or 'C' */
/*                            and  JOB = 'B' or 'D'; */
/*             LDC >= 1,  otherwise. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             If  M > 0,  JOB = 'A', 'B', or 'D'  and  INFO = 0,  the */
/*             leading N-by-M part of this array contains the system */
/*             input matrix. If  M = 0  or  JOB = 'C',  this array is */
/*             not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= N,  if M > 0 and JOB = 'A', 'B', or 'D'; */
/*             LDB >= 1,  if M = 0 or  JOB = 'C'. */

/*     D       (output) DOUBLE PRECISION array, dimension (LDD,M) */
/*             If  M > 0,  JOB = 'A' or 'D'  and  INFO = 0,  the leading */
/*             L-by-M part of this array contains the system input-output */
/*             matrix. If  M = 0  or  JOB = 'C' or 'B',  this array is */
/*             not referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D. */
/*             LDD >= L,  if M > 0 and JOB = 'A' or 'D'; */
/*             LDD >= 1,  if M = 0 or  JOB = 'C' or 'B'. */

/*     Q       (output) DOUBLE PRECISION array, dimension (LDQ,N) */
/*             If  JOBCK = 'C' or 'K',  the leading N-by-N part of this */
/*             array contains the positive semidefinite state covariance */
/*             matrix. If  JOBCK = 'K',  this matrix has been used as */
/*             state weighting matrix for computing the Kalman gain. */
/*             This parameter is not referenced if JOBCK = 'N'. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q. */
/*             LDQ >= N,  if JOBCK = 'C' or 'K'; */
/*             LDQ >= 1,  if JOBCK = 'N'. */

/*     RY      (output) DOUBLE PRECISION array, dimension (LDRY,L) */
/*             If  JOBCK = 'C' or 'K',  the leading L-by-L part of this */
/*             array contains the positive (semi)definite output */
/*             covariance matrix. If  JOBCK = 'K',  this matrix has been */
/*             used as output weighting matrix for computing the Kalman */
/*             gain. */
/*             This parameter is not referenced if JOBCK = 'N'. */

/*     LDRY    INTEGER */
/*             The leading dimension of the array RY. */
/*             LDRY >= L,  if JOBCK = 'C' or 'K'; */
/*             LDRY >= 1,  if JOBCK = 'N'. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,L) */
/*             If  JOBCK = 'C' or 'K',  the leading N-by-L part of this */
/*             array contains the state-output cross-covariance matrix. */
/*             If  JOBCK = 'K',  this matrix has been used as state- */
/*             output weighting matrix for computing the Kalman gain. */
/*             This parameter is not referenced if JOBCK = 'N'. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S. */
/*             LDS >= N,  if JOBCK = 'C' or 'K'; */
/*             LDS >= 1,  if JOBCK = 'N'. */

/*     K       (output) DOUBLE PRECISION array, dimension ( LDK,L ) */
/*             If  JOBCK = 'K',  the leading  N-by-L  part of this array */
/*             contains the estimated Kalman gain matrix. */
/*             If  JOBCK = 'C' or 'N',  this array is not referenced. */

/*     LDK     INTEGER */
/*             The leading dimension of the array  K. */
/*             LDK >= N,  if JOBCK = 'K'; */
/*             LDK >= 1,  if JOBCK = 'C' or 'N'. */

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

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK) */
/*             LIWORK >= max(LIW1,LIW2), where */
/*             LIW1 = N,                     if METH <> 'N' and M = 0 */
/*                                        or JOB = 'C' and JOBCK = 'N'; */
/*             LIW1 = M*NOBR+N,              if METH <> 'N', JOB = 'C', */
/*                                           and JOBCK <> 'N'; */
/*             LIW1 = max(L*NOBR,M*NOBR),    if METH = 'M', JOB <> 'C', */
/*                                           and JOBCK = 'N'; */
/*             LIW1 = max(L*NOBR,M*NOBR+N),  if METH = 'M', JOB <> 'C', */
/*                                           and JOBCK = 'C' or 'K'; */
/*             LIW1 = max(M*NOBR+N,M*(N+L)), if METH = 'N', or METH = 'C' */
/*                                           and JOB  <> 'C'; */
/*             LIW2 = 0,                     if JOBCK <> 'K'; */
/*             LIW2 = N*N,                   if JOBCK =  'K'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK,  and  DWORK(2),  DWORK(3),  DWORK(4),  and */
/*             DWORK(5)  contain the reciprocal condition numbers of the */
/*             triangular factors of the following matrices (defined in */
/*             SLICOT Library routine IB01PD and in the lower level */
/*             routines): */
/*                GaL  (GaL = Un(1:(s-1)*L,1:n)), */
/*                R_1c (if  METH = 'M' or 'C'), */
/*                M    (if  JOBCK = 'C' or 'K'  or  METH = 'N'),  and */
/*                Q or T  (see SLICOT Library routine IB01PY or IB01PX), */
/*             respectively. */
/*             If  METH = 'N',  DWORK(3)  is set to one without any */
/*             calculations. Similarly, if  METH = 'M'  and  JOBCK = 'N', */
/*             DWORK(4)  is set to one. If  M = 0  or  JOB = 'C', */
/*             DWORK(3)  and  DWORK(5)  are set to one. */
/*             If  JOBCK = 'K'  and  INFO = 0,  DWORK(6)  to  DWORK(13) */
/*             contain information about the accuracy of the results when */
/*             computing the Kalman gain matrix, as follows: */
/*                DWORK(6)  - reciprocal condition number of the matrix */
/*                            U11  of the Nth order system of algebraic */
/*                            equations from which the solution matrix  X */
/*                            of the Riccati equation is obtained; */
/*                DWORK(7)  - reciprocal pivot growth factor for the LU */
/*                            factorization of the matrix  U11; */
/*                DWORK(8)  - reciprocal condition number of the matrix */
/*                            As = A - S*inv(Ry)*C,  which is inverted by */
/*                            the standard Riccati solver; */
/*                DWORK(9)  - reciprocal pivot growth factor for the LU */
/*                            factorization of the matrix  As; */
/*                DWORK(10) - reciprocal condition number of the matrix */
/*                            Ry; */
/*                DWORK(11) - reciprocal condition number of the matrix */
/*                            Ry + C*X*C'; */
/*                DWORK(12) - reciprocal condition number for the Riccati */
/*                            equation solution; */
/*                DWORK(13) - forward error bound for the Riccati */
/*                            equation solution. */
/*             On exit, if  INFO = -30,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( LDW1,LDW2,LDW3 ), where, if METH = 'M', */
/*             LDW1 >= max( 2*(L*NOBR-L)*N+2*N, (L*NOBR-L)*N+N*N+7*N ), */
/*                     if JOB = 'C' or JOB = 'A' and M = 0; */
/*             LDW1 >= max( 2*(L*NOBR-L)*N+N*N+7*N, */
/*                          (L*NOBR-L)*N+N+6*M*NOBR, (L*NOBR-L)*N+N+ */
/*                          max( L+M*NOBR, L*NOBR + */
/*                                         max( 3*L*NOBR+1, M ) ) ), */
/*                     if M > 0 and JOB = 'A', 'B', or 'D'; */
/*             LDW2 >= 0,                          if JOBCK = 'N'; */
/*             LDW2 >= L*NOBR*N+ */
/*                     max( (L*NOBR-L)*N+Aw+2*N+max(5*N,(2*M+L)*NOBR+L), */
/*                          4*(M*NOBR+N)+1, M*NOBR+2*N+L ), */
/*                                                 if JOBCK = 'C' or 'K', */
/*             where Aw = N+N*N, if M = 0 or JOB = 'C'; */
/*                   Aw = 0,     otherwise; */
/*             if METH = 'N', */
/*             LDW1 >= L*NOBR*N+max( (L*NOBR-L)*N+2*N+(2*M+L)*NOBR+L, */
/*                                   2*(L*NOBR-L)*N+N*N+8*N, */
/*                                   N+4*(M*NOBR+N)+1, M*NOBR+3*N+L ); */
/*             LDW2 >= 0, if M = 0 or JOB = 'C'; */
/*             LDW2 >= L*NOBR*N+M*NOBR*(N+L)*(M*(N+L)+1)+ */
/*                                max( (N+L)**2, 4*M*(N+L)+1 ), */
/*                     if M > 0 and JOB = 'A', 'B', or 'D'; */
/*             and, if METH = 'C', LDW1 as */
/*             max( LDW1 for METH = 'M', JOB = 'C', LDW1 for METH = 'N'), */
/*             and LDW2 for METH = 'N' are used; */
/*             LDW3 >= 0,                     if JOBCK <> 'K'; */
/*             LDW3 >= max(  4*N*N+2*N*L+L*L+max( 3*L,N*L ), */
/*                          14*N*N+12*N+5 ),  if JOBCK =  'K'. */
/*             For good performance,  LDWORK  should be larger. */

/*     BWORK   LOGICAL array, dimension (LBWORK) */
/*             LBWORK = 2*N, if JOBCK =  'K'; */
/*             LBWORK = 0,   if JOBCK <> 'K'. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 4:  a least squares problem to be solved has a */
/*                   rank-deficient coefficient matrix; */
/*             = 5:  the computed covariance matrices are too small. */
/*                   The problem seems to be a deterministic one; the */
/*                   gain matrix is set to zero. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  the singular value decomposition (SVD) algorithm did */
/*                   not converge; */
/*             = 3:  a singular upper triangular matrix was found; */
/*             = 3+i:  if  JOBCK = 'K'  and the associated Riccati */
/*                   equation could not be solved, where i = 1,...,6; */
/*                   (see the description of the parameter INFO for the */
/*                   SLICOT Library routine SB02RD for the meaning of */
/*                   the i values); */
/*             = 10: the QR algorithm did not converge. */

/*     METHOD */

/*     In the MOESP approach, the matrices  A  and  C  are first */
/*     computed from an estimated extended observability matrix [1], */
/*     and then, the matrices  B  and  D  are obtained by solving an */
/*     extended linear system in a least squares sense. */
/*     In the N4SID approach, besides the estimated extended */
/*     observability matrix, the solutions of two least squares problems */
/*     are used to build another least squares problem, whose solution */
/*     is needed to compute the system matrices  A,  C,  B,  and  D.  The */
/*     solutions of the two least squares problems are also optionally */
/*     used by both approaches to find the covariance matrices. */
/*     The Kalman gain matrix is obtained by solving a discrete-time */
/*     algebraic Riccati equation. */

/*     REFERENCES */

/*     [1] Verhaegen M., and Dewilde, P. */
/*         Subspace Model Identification. Part 1: The output-error */
/*         state-space model identification class of algorithms. */
/*         Int. J. Control, 56, pp. 1187-1210, 1992. */

/*     [2] Van Overschee, P., and De Moor, B. */
/*         N4SID: Two Subspace Algorithms for the Identification */
/*         of Combined Deterministic-Stochastic Systems. */
/*         Automatica, Vol.30, No.1, pp. 75-93, 1994. */

/*     [3] Van Overschee, P. */
/*         Subspace Identification : Theory - Implementation - */
/*         Applications. */
/*         Ph. D. Thesis, Department of Electrical Engineering, */
/*         Katholieke Universiteit Leuven, Belgium, Feb. 1995. */

/*     [4] Sima, V. */
/*         Subspace-based Algorithms for Multivariable System */
/*         Identification. */
/*         Studies in Informatics and Control, 5, pp. 335-344, 1996. */

/*     NUMERICAL ASPECTS */

/*     The implemented method consists in numerically stable steps. */

/*     FURTHER COMMENTS */

/*     The covariance matrices are computed using the N4SID approach. */
/*     Therefore, for efficiency reasons, it is advisable to set */
/*     METH = 'N',  if the Kalman gain matrix or covariance matrices */
/*     are needed  (JOBCK = 'K', or 'C').  When  JOBCK = 'N',  it could */
/*     be more efficient to use the combined method,  METH = 'C'. */
/*     Often, this combination will also provide better accuracy than */
/*     MOESP algorithm. */
/*     In some applications, it is useful to compute the system matrices */
/*     using two calls to this routine, the first one with  JOB = 'C', */
/*     and the second one with  JOB = 'B' or 'D'.  This is slightly less */
/*     efficient than using a single call with  JOB = 'A',  because some */
/*     calculations are repeated. If  METH = 'N',  all the calculations */
/*     at the first call are performed again at the second call; */
/*     moreover, it is required to save the needed submatrices of  R */
/*     before the first call and restore them before the second call. */
/*     If the covariance matrices and/or the Kalman gain are desired, */
/*     JOBCK  should be set to  'C'  or  'K'  at the second call. */
/*     If  B  and  D  are both needed, they should be computed at once. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 1999. */

/*     REVISIONS */

/*     March 2000, August 2000, Sept. 2001, March 2005. */

/*     KEYWORDS */

/*     Identification methods; least squares solutions; multivariable */
/*     systems; QR decomposition; singular value decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Arrays .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 471 "IB01BD.f"
    /* Parameter adjustments */
#line 471 "IB01BD.f"
    r_dim1 = *ldr;
#line 471 "IB01BD.f"
    r_offset = 1 + r_dim1;
#line 471 "IB01BD.f"
    r__ -= r_offset;
#line 471 "IB01BD.f"
    a_dim1 = *lda;
#line 471 "IB01BD.f"
    a_offset = 1 + a_dim1;
#line 471 "IB01BD.f"
    a -= a_offset;
#line 471 "IB01BD.f"
    c_dim1 = *ldc;
#line 471 "IB01BD.f"
    c_offset = 1 + c_dim1;
#line 471 "IB01BD.f"
    c__ -= c_offset;
#line 471 "IB01BD.f"
    b_dim1 = *ldb;
#line 471 "IB01BD.f"
    b_offset = 1 + b_dim1;
#line 471 "IB01BD.f"
    b -= b_offset;
#line 471 "IB01BD.f"
    d_dim1 = *ldd;
#line 471 "IB01BD.f"
    d_offset = 1 + d_dim1;
#line 471 "IB01BD.f"
    d__ -= d_offset;
#line 471 "IB01BD.f"
    q_dim1 = *ldq;
#line 471 "IB01BD.f"
    q_offset = 1 + q_dim1;
#line 471 "IB01BD.f"
    q -= q_offset;
#line 471 "IB01BD.f"
    ry_dim1 = *ldry;
#line 471 "IB01BD.f"
    ry_offset = 1 + ry_dim1;
#line 471 "IB01BD.f"
    ry -= ry_offset;
#line 471 "IB01BD.f"
    s_dim1 = *lds;
#line 471 "IB01BD.f"
    s_offset = 1 + s_dim1;
#line 471 "IB01BD.f"
    s -= s_offset;
#line 471 "IB01BD.f"
    k_dim1 = *ldk;
#line 471 "IB01BD.f"
    k_offset = 1 + k_dim1;
#line 471 "IB01BD.f"
    k -= k_offset;
#line 471 "IB01BD.f"
    --iwork;
#line 471 "IB01BD.f"
    --dwork;
#line 471 "IB01BD.f"
    --bwork;
#line 471 "IB01BD.f"

#line 471 "IB01BD.f"
    /* Function Body */
#line 471 "IB01BD.f"
    moesp = lsame_(meth, "M", (ftnlen)1, (ftnlen)1);
#line 472 "IB01BD.f"
    n4sid = lsame_(meth, "N", (ftnlen)1, (ftnlen)1);
#line 473 "IB01BD.f"
    combin = lsame_(meth, "C", (ftnlen)1, (ftnlen)1);
#line 474 "IB01BD.f"
    withal = lsame_(job, "A", (ftnlen)1, (ftnlen)1);
#line 475 "IB01BD.f"
    withc = lsame_(job, "C", (ftnlen)1, (ftnlen)1) || withal;
#line 476 "IB01BD.f"
    withd = lsame_(job, "D", (ftnlen)1, (ftnlen)1) || withal;
#line 477 "IB01BD.f"
    withb = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || withd;
#line 478 "IB01BD.f"
    withk = lsame_(jobck, "K", (ftnlen)1, (ftnlen)1);
#line 479 "IB01BD.f"
    withco = lsame_(jobck, "C", (ftnlen)1, (ftnlen)1) || withk;
#line 480 "IB01BD.f"
    mnobr = *m * *nobr;
#line 481 "IB01BD.f"
    lnobr = *l * *nobr;
#line 482 "IB01BD.f"
    lmnobr = lnobr + mnobr;
#line 483 "IB01BD.f"
    mnobrn = mnobr + *n;
#line 484 "IB01BD.f"
    ldunn = (lnobr - *l) * *n;
#line 485 "IB01BD.f"
    lmmnol = lnobr + (mnobr << 1) + *l;
#line 486 "IB01BD.f"
    nr = lmnobr + lmnobr;
#line 487 "IB01BD.f"
    npl = *n + *l;
#line 488 "IB01BD.f"
    n2 = *n + *n;
#line 489 "IB01BD.f"
    nn = *n * *n;
#line 490 "IB01BD.f"
    nl = *n * *l;
#line 491 "IB01BD.f"
    ll = *l * *l;
#line 492 "IB01BD.f"
    minwrk = 1;
#line 493 "IB01BD.f"
    *iwarn = 0;
#line 494 "IB01BD.f"
    *info = 0;

/*     Check the scalar input parameters. */

#line 498 "IB01BD.f"
    if (! (moesp || n4sid || combin)) {
#line 499 "IB01BD.f"
	*info = -1;
#line 500 "IB01BD.f"
    } else if (! (withb || withc)) {
#line 501 "IB01BD.f"
	*info = -2;
#line 502 "IB01BD.f"
    } else if (! (withco || lsame_(jobck, "N", (ftnlen)1, (ftnlen)1))) {
#line 503 "IB01BD.f"
	*info = -3;
#line 504 "IB01BD.f"
    } else if (*nobr <= 1) {
#line 505 "IB01BD.f"
	*info = -4;
#line 506 "IB01BD.f"
    } else if (*n <= 0 || *n >= *nobr) {
#line 507 "IB01BD.f"
	*info = -5;
#line 508 "IB01BD.f"
    } else if (*m < 0) {
#line 509 "IB01BD.f"
	*info = -6;
#line 510 "IB01BD.f"
    } else if (*l <= 0) {
#line 511 "IB01BD.f"
	*info = -7;
#line 512 "IB01BD.f"
    } else if (withco && *nsmpl < nr) {
#line 513 "IB01BD.f"
	*info = -8;
#line 514 "IB01BD.f"
    } else if (*ldr < nr) {
#line 515 "IB01BD.f"
	*info = -10;
#line 516 "IB01BD.f"
    } else if (*lda < 1 || (withc || withb && ! moesp) && *lda < *n) {
#line 518 "IB01BD.f"
	*info = -12;
#line 519 "IB01BD.f"
    } else if (*ldc < 1 || (withc || withb && ! moesp) && *ldc < *l) {
#line 521 "IB01BD.f"
	*info = -14;
#line 522 "IB01BD.f"
    } else if (*ldb < 1 || withb && *ldb < *n && *m > 0) {
#line 524 "IB01BD.f"
	*info = -16;
#line 525 "IB01BD.f"
    } else if (*ldd < 1 || withd && *ldd < *l && *m > 0) {
#line 527 "IB01BD.f"
	*info = -18;
#line 528 "IB01BD.f"
    } else if (*ldq < 1 || withco && *ldq < *n) {
#line 529 "IB01BD.f"
	*info = -20;
#line 530 "IB01BD.f"
    } else if (*ldry < 1 || withco && *ldry < *l) {
#line 531 "IB01BD.f"
	*info = -22;
#line 532 "IB01BD.f"
    } else if (*lds < 1 || withco && *lds < *n) {
#line 533 "IB01BD.f"
	*info = -24;
#line 534 "IB01BD.f"
    } else if (*ldk < 1 || withk && *ldk < *n) {
#line 535 "IB01BD.f"
	*info = -26;
#line 536 "IB01BD.f"
    } else {

/*        Compute workspace. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*         minimal amount of workspace needed at that point in the code, */
/*         as well as the preferred amount for good performance.) */

#line 543 "IB01BD.f"
	iaw = 0;
#line 544 "IB01BD.f"
	minwrk = ldunn + (*n << 2);
#line 545 "IB01BD.f"
	if (! n4sid) {
#line 546 "IB01BD.f"
	    id = 0;
#line 547 "IB01BD.f"
	    if (withc) {
/* Computing MAX */
#line 548 "IB01BD.f"
		i__1 = minwrk, i__2 = (ldunn << 1) + n2, i__1 = max(i__1,i__2)
			, i__2 = ldunn + nn + *n * 7;
#line 548 "IB01BD.f"
		minwrk = max(i__1,i__2);
#line 549 "IB01BD.f"
	    }
#line 550 "IB01BD.f"
	} else {
#line 551 "IB01BD.f"
	    id = *n;
#line 552 "IB01BD.f"
	}

#line 554 "IB01BD.f"
	if (*m > 0 && withb || ! moesp) {
/* Computing MAX */
#line 555 "IB01BD.f"
	    i__1 = minwrk, i__2 = (ldunn << 1) + nn + id + *n * 7;
#line 555 "IB01BD.f"
	    minwrk = max(i__1,i__2);
#line 556 "IB01BD.f"
	    if (moesp) {
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 556 "IB01BD.f"
		i__5 = lnobr * 3 + 1;
#line 556 "IB01BD.f"
		i__3 = *l + mnobr, i__4 = lnobr + max(i__5,*m);
#line 556 "IB01BD.f"
		i__1 = minwrk, i__2 = ldunn + *n + mnobr * 6, i__1 = max(i__1,
			i__2), i__2 = ldunn + *n + max(i__3,i__4);
#line 556 "IB01BD.f"
		minwrk = max(i__1,i__2);
#line 556 "IB01BD.f"
	    }
#line 560 "IB01BD.f"
	} else {
#line 561 "IB01BD.f"
	    if (! n4sid) {
#line 561 "IB01BD.f"
		iaw = *n + nn;
#line 561 "IB01BD.f"
	    }
#line 563 "IB01BD.f"
	}

#line 565 "IB01BD.f"
	if (! moesp || withco) {
/* Computing MAX */
/* Computing MAX */
#line 566 "IB01BD.f"
	    i__3 = *n * 5;
#line 566 "IB01BD.f"
	    i__1 = minwrk, i__2 = ldunn + iaw + n2 + max(i__3,lmmnol), i__1 = 
		    max(i__1,i__2), i__2 = id + (mnobrn << 2) + 1, i__1 = max(
		    i__1,i__2), i__2 = id + mnobrn + npl;
#line 566 "IB01BD.f"
	    minwrk = max(i__1,i__2);
#line 568 "IB01BD.f"
	    if (! moesp && *m > 0 && withb) {
/* Computing MAX */
/* Computing MAX */
/* Computing 2nd power */
#line 568 "IB01BD.f"
		i__5 = npl;
#line 568 "IB01BD.f"
		i__3 = i__5 * i__5, i__4 = (*m << 2) * npl + 1;
#line 568 "IB01BD.f"
		i__1 = minwrk, i__2 = mnobr * npl * (*m * npl + 1) + max(i__3,
			i__4);
#line 568 "IB01BD.f"
		minwrk = max(i__1,i__2);
#line 568 "IB01BD.f"
	    }
#line 571 "IB01BD.f"
	    minwrk = lnobr * *n + minwrk;
#line 572 "IB01BD.f"
	}

#line 574 "IB01BD.f"
	if (withk) {
/* Computing MAX */
/* Computing MAX */
#line 575 "IB01BD.f"
	    i__3 = *l * 3;
#line 575 "IB01BD.f"
	    i__1 = minwrk, i__2 = (nn << 2) + (nl << 1) + ll + max(i__3,nl), 
		    i__1 = max(i__1,i__2), i__2 = nn * 14 + *n * 12 + 5;
#line 575 "IB01BD.f"
	    minwrk = max(i__1,i__2);
#line 577 "IB01BD.f"
	}

#line 579 "IB01BD.f"
	if (*ldwork < minwrk) {
#line 580 "IB01BD.f"
	    *info = -30;
#line 581 "IB01BD.f"
	    dwork[1] = (doublereal) minwrk;
#line 582 "IB01BD.f"
	}
#line 583 "IB01BD.f"
    }

/*     Return if there are illegal arguments. */

#line 587 "IB01BD.f"
    if (*info != 0) {
#line 588 "IB01BD.f"
	i__1 = -(*info);
#line 588 "IB01BD.f"
	xerbla_("IB01BD", &i__1, (ftnlen)6);
#line 589 "IB01BD.f"
	return 0;
#line 590 "IB01BD.f"
    }

#line 592 "IB01BD.f"
    if (! withk) {
#line 593 "IB01BD.f"
	*(unsigned char *)jobcv = *(unsigned char *)jobck;
#line 594 "IB01BD.f"
    } else {
#line 595 "IB01BD.f"
	*(unsigned char *)jobcv = 'C';
#line 596 "IB01BD.f"
    }

#line 598 "IB01BD.f"
    io = 1;
#line 599 "IB01BD.f"
    if (! moesp || withco) {
#line 600 "IB01BD.f"
	jwork = io + lnobr * *n;
#line 601 "IB01BD.f"
    } else {
#line 602 "IB01BD.f"
	jwork = io;
#line 603 "IB01BD.f"
    }
#line 604 "IB01BD.f"
    maxwrk = minwrk;

/*     Call the computational routine for estimating system matrices. */

#line 608 "IB01BD.f"
    if (! combin) {
#line 609 "IB01BD.f"
	i__1 = *ldwork - jwork + 1;
#line 609 "IB01BD.f"
	ib01pd_(meth, job, jobcv, nobr, n, m, l, nsmpl, &r__[r_offset], ldr, &
		a[a_offset], lda, &c__[c_offset], ldc, &b[b_offset], ldb, &
		d__[d_offset], ldd, &q[q_offset], ldq, &ry[ry_offset], ldry, &
		s[s_offset], lds, &dwork[io], &lnobr, tol, &iwork[1], &dwork[
		jwork], &i__1, iwarn, info, (ftnlen)1, (ftnlen)1, (ftnlen)1);

#line 614 "IB01BD.f"
    } else {

#line 616 "IB01BD.f"
	if (withc) {
#line 617 "IB01BD.f"
	    if (withal) {
#line 618 "IB01BD.f"
		*(unsigned char *)jobcov = 'N';
#line 619 "IB01BD.f"
	    } else {
#line 620 "IB01BD.f"
		*(unsigned char *)jobcov = *(unsigned char *)jobcv;
#line 621 "IB01BD.f"
	    }
#line 622 "IB01BD.f"
	    i__1 = *ldwork - jwork + 1;
#line 622 "IB01BD.f"
	    ib01pd_("MOESP", "C and A", jobcov, nobr, n, m, l, nsmpl, &r__[
		    r_offset], ldr, &a[a_offset], lda, &c__[c_offset], ldc, &
		    b[b_offset], ldb, &d__[d_offset], ldd, &q[q_offset], ldq, 
		    &ry[ry_offset], ldry, &s[s_offset], lds, &dwork[io], &
		    lnobr, tol, &iwork[1], &dwork[jwork], &i__1, &iwarnl, 
		    info, (ftnlen)5, (ftnlen)7, (ftnlen)1);
#line 627 "IB01BD.f"
	    if (*info != 0) {
#line 627 "IB01BD.f"
		return 0;
#line 627 "IB01BD.f"
	    }
#line 629 "IB01BD.f"
	    *iwarn = max(*iwarn,iwarnl);
/* Computing MAX */
#line 630 "IB01BD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 630 "IB01BD.f"
	    maxwrk = max(i__1,i__2);
#line 631 "IB01BD.f"
	}

#line 633 "IB01BD.f"
	if (withb) {
#line 634 "IB01BD.f"
	    if (! withal) {
#line 635 "IB01BD.f"
		*(unsigned char *)jobbd = *(unsigned char *)job;
#line 636 "IB01BD.f"
	    } else {
#line 637 "IB01BD.f"
		*(unsigned char *)jobbd = 'D';
#line 638 "IB01BD.f"
	    }
#line 639 "IB01BD.f"
	    i__1 = *ldwork - jwork + 1;
#line 639 "IB01BD.f"
	    ib01pd_("N4SID", jobbd, jobcv, nobr, n, m, l, nsmpl, &r__[
		    r_offset], ldr, &a[a_offset], lda, &c__[c_offset], ldc, &
		    b[b_offset], ldb, &d__[d_offset], ldd, &q[q_offset], ldq, 
		    &ry[ry_offset], ldry, &s[s_offset], lds, &dwork[io], &
		    lnobr, tol, &iwork[1], &dwork[jwork], &i__1, &iwarnl, 
		    info, (ftnlen)5, (ftnlen)1, (ftnlen)1);
#line 643 "IB01BD.f"
	    *iwarn = max(*iwarn,iwarnl);
#line 644 "IB01BD.f"
	}
#line 645 "IB01BD.f"
    }

#line 647 "IB01BD.f"
    if (*info != 0) {
#line 647 "IB01BD.f"
	return 0;
#line 647 "IB01BD.f"
    }
/* Computing MAX */
#line 649 "IB01BD.f"
    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 649 "IB01BD.f"
    maxwrk = max(i__1,i__2);

#line 651 "IB01BD.f"
    for (i__ = 1; i__ <= 4; ++i__) {
#line 652 "IB01BD.f"
	rcnd[i__ - 1] = dwork[jwork + i__];
#line 653 "IB01BD.f"
/* L10: */
#line 653 "IB01BD.f"
    }

#line 655 "IB01BD.f"
    if (withk) {
#line 656 "IB01BD.f"
	if (*iwarn == 5) {

/*           The problem seems to be a deterministic one. Set the Kalman */
/*           gain to zero, set accuracy parameters and return. */

#line 661 "IB01BD.f"
	    dlaset_("Full", n, l, &c_b21, &c_b21, &k[k_offset], ldk, (ftnlen)
		    4);

#line 663 "IB01BD.f"
	    for (i__ = 6; i__ <= 12; ++i__) {
#line 664 "IB01BD.f"
		dwork[i__] = 1.;
#line 665 "IB01BD.f"
/* L20: */
#line 665 "IB01BD.f"
	    }

#line 667 "IB01BD.f"
	    dwork[13] = 0.;
#line 668 "IB01BD.f"
	} else {

/*           Compute the Kalman gain matrix. */

/*           Convert the optimal problem with coupling weighting terms */
/*           to a standard problem. */
/*           Workspace:  need   4*N*N+2*N*L+L*L+max( 3*L,N*L ); */
/*                       prefer larger. */

#line 677 "IB01BD.f"
	    ix = 1;
#line 678 "IB01BD.f"
	    iq = ix + nn;
#line 679 "IB01BD.f"
	    ia = iq + nn;
#line 680 "IB01BD.f"
	    ig = ia + nn;
#line 681 "IB01BD.f"
	    ic = ig + nn;
#line 682 "IB01BD.f"
	    ir = ic + nl;
#line 683 "IB01BD.f"
	    is = ir + ll;
#line 684 "IB01BD.f"
	    jwork = is + nl;

#line 686 "IB01BD.f"
	    ma02ad_("Full", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)4)
		    ;
#line 687 "IB01BD.f"
	    ma02ad_("Full", l, n, &c__[c_offset], ldc, &dwork[ic], n, (ftnlen)
		    4);
#line 688 "IB01BD.f"
	    dlacpy_("Upper", n, n, &q[q_offset], ldq, &dwork[iq], n, (ftnlen)
		    5);
#line 689 "IB01BD.f"
	    dlacpy_("Upper", l, l, &ry[ry_offset], ldry, &dwork[ir], l, (
		    ftnlen)5);
#line 690 "IB01BD.f"
	    dlacpy_("Full", n, l, &s[s_offset], lds, &dwork[is], n, (ftnlen)4)
		    ;

#line 692 "IB01BD.f"
	    i__1 = *ldwork - jwork + 1;
#line 692 "IB01BD.f"
	    sb02mt_("G needed", "Nonzero S", "Not factored", "Upper", n, l, &
		    dwork[ia], n, &dwork[ic], n, &dwork[iq], n, &dwork[ir], l,
		     &dwork[is], n, &iwork[1], &ifact, &dwork[ig], n, &iwork[*
		    l + 1], &dwork[jwork], &i__1, &ierr, (ftnlen)8, (ftnlen)9,
		     (ftnlen)12, (ftnlen)5);
#line 697 "IB01BD.f"
	    if (ierr != 0) {
#line 698 "IB01BD.f"
		*info = 3;
#line 699 "IB01BD.f"
		return 0;
#line 700 "IB01BD.f"
	    }
/* Computing MAX */
#line 701 "IB01BD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 701 "IB01BD.f"
	    maxwrk = max(i__1,i__2);
#line 702 "IB01BD.f"
	    rcondr = dwork[jwork + 1];

/*           Solve the Riccati equation. */
/*           Workspace:  need   14*N*N+12*N+5; */
/*                       prefer larger. */

#line 708 "IB01BD.f"
	    it = ic;
#line 709 "IB01BD.f"
	    iv = it + nn;
#line 710 "IB01BD.f"
	    iwr = iv + nn;
#line 711 "IB01BD.f"
	    iwi = iwr + n2;
#line 712 "IB01BD.f"
	    is = iwi + n2;
#line 713 "IB01BD.f"
	    jwork = is + n2 * n2;

#line 715 "IB01BD.f"
	    i__1 = *ldwork - jwork + 1;
#line 715 "IB01BD.f"
	    sb02rd_("All", "Discrete", "Direct", "NoTranspose", "Upper", 
		    "General scaling", "Unstable first", "Not factored", 
		    "Reduced", n, &dwork[ia], n, &dwork[it], n, &dwork[iv], n,
		     &dwork[ig], n, &dwork[iq], n, &dwork[ix], n, &sep, &
		    rcond, &ferr, &dwork[iwr], &dwork[iwi], &dwork[is], &n2, &
		    iwork[1], &dwork[jwork], &i__1, &bwork[1], &ierr, (ftnlen)
		    3, (ftnlen)8, (ftnlen)6, (ftnlen)11, (ftnlen)5, (ftnlen)
		    15, (ftnlen)14, (ftnlen)12, (ftnlen)7);

#line 723 "IB01BD.f"
	    if (ierr != 0 && ierr < 7) {
#line 724 "IB01BD.f"
		*info = ierr + 3;
#line 725 "IB01BD.f"
		return 0;
#line 726 "IB01BD.f"
	    }
/* Computing MAX */
#line 727 "IB01BD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 727 "IB01BD.f"
	    maxwrk = max(i__1,i__2);

#line 729 "IB01BD.f"
	    for (i__ = 1; i__ <= 4; ++i__) {
#line 730 "IB01BD.f"
		rcnd[i__ + 3] = dwork[jwork + i__];
#line 731 "IB01BD.f"
/* L30: */
#line 731 "IB01BD.f"
	    }

/*           Compute the gain matrix. */
/*           Workspace:  need   2*N*N+2*N*L+L*L+3*L; */
/*                       prefer larger. */

#line 737 "IB01BD.f"
	    ia = ix + nn;
#line 738 "IB01BD.f"
	    ic = ia + nn;
#line 739 "IB01BD.f"
	    ir = ic + nl;
#line 740 "IB01BD.f"
	    ik = ir + ll;
#line 741 "IB01BD.f"
	    jwork = ik + nl;

#line 743 "IB01BD.f"
	    ma02ad_("Full", n, n, &a[a_offset], lda, &dwork[ia], n, (ftnlen)4)
		    ;
#line 744 "IB01BD.f"
	    ma02ad_("Full", l, n, &c__[c_offset], ldc, &dwork[ic], n, (ftnlen)
		    4);
#line 745 "IB01BD.f"
	    dlacpy_("Upper", l, l, &ry[ry_offset], ldry, &dwork[ir], l, (
		    ftnlen)5);

#line 747 "IB01BD.f"
	    i__1 = *ldwork - jwork + 1;
#line 747 "IB01BD.f"
	    sb02nd_("Discrete", "NotFactored", "Upper", "Nonzero S", n, l, &
		    c__0, &dwork[ia], n, &dwork[ic], n, &dwork[ir], l, &iwork[
		    1], &s[s_offset], lds, &dwork[ix], n, &rnorm, &dwork[ik], 
		    l, oufact, &iwork[*l + 1], &dwork[jwork], &i__1, &ierr, (
		    ftnlen)8, (ftnlen)11, (ftnlen)5, (ftnlen)9);

#line 753 "IB01BD.f"
	    if (ierr != 0) {
#line 754 "IB01BD.f"
		if (ierr <= *l + 1) {
#line 755 "IB01BD.f"
		    *info = 3;
#line 756 "IB01BD.f"
		} else if (ierr == *l + 2) {
#line 757 "IB01BD.f"
		    *info = 10;
#line 758 "IB01BD.f"
		}
#line 759 "IB01BD.f"
		return 0;
#line 760 "IB01BD.f"
	    }
/* Computing MAX */
#line 761 "IB01BD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 761 "IB01BD.f"
	    maxwrk = max(i__1,i__2);

#line 763 "IB01BD.f"
	    ma02ad_("Full", l, n, &dwork[ik], l, &k[k_offset], ldk, (ftnlen)4)
		    ;

/*           Set the accuracy parameters. */

#line 767 "IB01BD.f"
	    dwork[11] = dwork[jwork + 1];

#line 769 "IB01BD.f"
	    for (i__ = 6; i__ <= 9; ++i__) {
#line 770 "IB01BD.f"
		dwork[i__] = rcnd[i__ - 2];
#line 771 "IB01BD.f"
/* L40: */
#line 771 "IB01BD.f"
	    }

#line 773 "IB01BD.f"
	    dwork[10] = rcondr;
#line 774 "IB01BD.f"
	    dwork[12] = rcond;
#line 775 "IB01BD.f"
	    dwork[13] = ferr;
#line 776 "IB01BD.f"
	}
#line 777 "IB01BD.f"
    }

/*     Return optimal workspace in  DWORK(1)  and the remaining */
/*     reciprocal condition numbers in the next locations. */

#line 782 "IB01BD.f"
    dwork[1] = (doublereal) maxwrk;

#line 784 "IB01BD.f"
    for (i__ = 2; i__ <= 5; ++i__) {
#line 785 "IB01BD.f"
	dwork[i__] = rcnd[i__ - 2];
#line 786 "IB01BD.f"
/* L50: */
#line 786 "IB01BD.f"
    }

#line 788 "IB01BD.f"
    return 0;

/* *** Last line of IB01BD *** */
} /* ib01bd_ */

