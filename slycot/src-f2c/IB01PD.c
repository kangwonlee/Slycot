#line 1 "IB01PD.f"
/* IB01PD.f -- translated by f2c (version 20100827).
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

#line 1 "IB01PD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static doublereal c_b43 = .66666666666666663;
static doublereal c_b63 = 1.;
static doublereal c_b66 = 0.;
static doublereal c_b173 = -1.;

/* Subroutine */ int ib01pd_(char *meth, char *job, char *jobcv, integer *
	nobr, integer *n, integer *m, integer *l, integer *nsmpl, doublereal *
	r__, integer *ldr, doublereal *a, integer *lda, doublereal *c__, 
	integer *ldc, doublereal *b, integer *ldb, doublereal *d__, integer *
	ldd, doublereal *q, integer *ldq, doublereal *ry, integer *ldry, 
	doublereal *s, integer *lds, doublereal *o, integer *ldo, doublereal *
	tol, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	iwarn, integer *info, ftnlen meth_len, ftnlen job_len, ftnlen 
	jobcv_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, o_dim1, o_offset, q_dim1, q_offset, r_dim1, r_offset, 
	    ry_dim1, ry_offset, s_dim1, s_offset, i__1, i__2, i__3, i__4, 
	    i__5;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, n2, id, nn, iu, nr, ix, nr2, nr3, nr4, iaw, ldw;
    static doublereal eps;
    static integer npl, isv, iun2, igal;
    static char fact[1], jobp[1];
    static integer ncol, rank, ierr, itau;
    static doublereal sval[3], toll, rnrm;
    static integer nrow;
    static logical n4sid;
    static integer itau1, itau2, ldun2;
    static doublereal toll1;
    static integer nr4mn, nr4pl;
    extern /* Subroutine */ int ma02ad_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    ma02ed_(char *, integer *, doublereal *, integer *, ftnlen), 
	    mb03od_(char *, integer *, integer *, doublereal *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, ftnlen), dgemm_(
	    char *, char *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, ftnlen, ftnlen), mb02ud_(char *, char *, 
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, ftnlen, ftnlen, ftnlen, 
	    ftnlen);
    static integer rank11;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int ib01px_(char *, integer *, integer *, integer 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen), ib01py_(
	    char *, char *, integer *, integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *, ftnlen, ftnlen);
    static integer rankm;
    extern /* Subroutine */ int mb02qy_(integer *, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *);
    static integer lnobr, mnobr;
    static logical shift, withb;
    static integer ldunn;
    static logical withc;
    static char jobpy[1];
    static logical fullr, moesp;
    static integer ihous;
    static logical withd;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static integer jwork;
    extern /* Subroutine */ int dsyrk_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen, ftnlen), dtrsm_(char *, char *, char *, char *
	    , integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen, ftnlen, ftnlen, ftnlen);
    static doublereal rcond1, rcond2, rcond3, rcond4;
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dgeqrf_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *), 
	    dlacpy_(char *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, ftnlen), dlaset_(char *, integer *, 
	    integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lmmnob;
    static logical withal;
    static integer lmnobr, lnobrn, mnobrn, iwarnl;
    static doublereal thresh;
    static integer lmmnol;
    static logical withco;
    extern /* Subroutine */ int dtrcon_(char *, char *, char *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen), dormqr_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen);
    static integer minwrk, maxwrk;
    static doublereal svlmax;
    extern /* Subroutine */ int dtrtrs_(char *, char *, char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *, ftnlen, ftnlen, ftnlen);


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

/*     To estimate the matrices A, C, B, and D of a linear time-invariant */
/*     (LTI) state space model, using the singular value decomposition */
/*     information provided by other routines. Optionally, the system and */
/*     noise covariance matrices, needed for the Kalman gain, are also */
/*     determined. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     METH    CHARACTER*1 */
/*             Specifies the subspace identification method to be used, */
/*             as follows: */
/*             = 'M':  MOESP  algorithm with past inputs and outputs; */
/*             = 'N':  N4SID  algorithm. */

/*     JOB     CHARACTER*1 */
/*             Specifies which matrices should be computed, as follows: */
/*             = 'A':  compute all system matrices, A, B, C, and D; */
/*             = 'C':  compute the matrices A and C only; */
/*             = 'B':  compute the matrix B only; */
/*             = 'D':  compute the matrices B and D only. */

/*     JOBCV   CHARACTER*1 */
/*             Specifies whether or not the covariance matrices are to */
/*             be computed, as follows: */
/*             = 'C':  the covariance matrices should be computed; */
/*             = 'N':  the covariance matrices should not be computed. */

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
/*             If JOBCV = 'C', the total number of samples used for */
/*             calculating the covariance matrices. */
/*             NSMPL >= 2*(M+L)*NOBR. */
/*             This parameter is not meaningful if  JOBCV = 'N'. */

/*     R       (input/workspace) DOUBLE PRECISION array, dimension */
/*             ( LDR,2*(M+L)*NOBR ) */
/*             On entry, the leading  2*(M+L)*NOBR-by-2*(M+L)*NOBR  part */
/*             of this array must contain the relevant data for the MOESP */
/*             or N4SID algorithms, as constructed by SLICOT Library */
/*             routines IB01AD or IB01ND. Let  R_ij,  i,j = 1:4,  be the */
/*             ij submatrix of  R  (denoted  S  in IB01AD and IB01ND), */
/*             partitioned by  M*NOBR,  L*NOBR,  M*NOBR,  and  L*NOBR */
/*             rows and columns. The submatrix  R_22  contains the matrix */
/*             of left singular vectors used. Also needed, for */
/*             METH = 'N'  or  JOBCV = 'C',  are the submatrices  R_11, */
/*             R_14 : R_44,  and, for  METH = 'M'  and  JOB <> 'C',  the */
/*             submatrices  R_31  and  R_12,  containing the processed */
/*             matrices  R_1c  and  R_2c,  respectively, as returned by */
/*             SLICOT Library routines IB01AD or IB01ND. */
/*             Moreover, if  METH = 'N'  and  JOB = 'A' or 'C',  the */
/*             block-row  R_41 : R_43  must contain the transpose of the */
/*             block-column  R_14 : R_34  as returned by SLICOT Library */
/*             routines IB01AD or IB01ND. */
/*             The remaining part of  R  is used as workspace. */
/*             On exit, part of this array is overwritten. Specifically, */
/*             if  METH = 'M',  R_22  and  R_31  are overwritten if */
/*                 JOB = 'B' or 'D',  and  R_12,  R_22,  R_14 : R_34, */
/*                 and possibly  R_11  are overwritten if  JOBCV = 'C'; */
/*             if  METH = 'N',  all needed submatrices are overwritten. */

/*     LDR     INTEGER */
/*             The leading dimension of the array  R. */
/*             LDR >= 2*(M+L)*NOBR. */

/*     A       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDA,N) */
/*             On entry, if  METH = 'N'  and  JOB = 'B' or 'D',  the */
/*             leading N-by-N part of this array must contain the system */
/*             state matrix. */
/*             If  METH = 'M'  or  (METH = 'N'  and JOB = 'A' or 'C'), */
/*             this array need not be set on input. */
/*             On exit, if  JOB = 'A' or 'C'  and  INFO = 0,  the */
/*             leading N-by-N part of this array contains the system */
/*             state matrix. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= N,  if  JOB = 'A' or 'C',  or  METH = 'N'  and */
/*                            JOB = 'B' or 'D'; */
/*             LDA >= 1,  otherwise. */

/*     C       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDC,N) */
/*             On entry, if  METH = 'N'  and  JOB = 'B' or 'D',  the */
/*             leading L-by-N part of this array must contain the system */
/*             output matrix. */
/*             If  METH = 'M'  or  (METH = 'N'  and JOB = 'A' or 'C'), */
/*             this array need not be set on input. */
/*             On exit, if  JOB = 'A' or 'C'  and  INFO = 0,  or */
/*             INFO = 3  (or  INFO >= 0,  for  METH = 'M'),  the leading */
/*             L-by-N part of this array contains the system output */
/*             matrix. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= L,  if  JOB = 'A' or 'C',  or  METH = 'N'  and */
/*                            JOB = 'B' or 'D'; */
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
/*             If JOBCV = 'C', the leading N-by-N part of this array */
/*             contains the positive semidefinite state covariance matrix */
/*             to be used as state weighting matrix when computing the */
/*             Kalman gain. */
/*             This parameter is not referenced if JOBCV = 'N'. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q. */
/*             LDQ >= N,  if JOBCV = 'C'; */
/*             LDQ >= 1,  if JOBCV = 'N'. */

/*     RY      (output) DOUBLE PRECISION array, dimension (LDRY,L) */
/*             If JOBCV = 'C', the leading L-by-L part of this array */
/*             contains the positive (semi)definite output covariance */
/*             matrix to be used as output weighting matrix when */
/*             computing the Kalman gain. */
/*             This parameter is not referenced if JOBCV = 'N'. */

/*     LDRY    INTEGER */
/*             The leading dimension of the array RY. */
/*             LDRY >= L,  if JOBCV = 'C'; */
/*             LDRY >= 1,  if JOBCV = 'N'. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,L) */
/*             If JOBCV = 'C', the leading N-by-L part of this array */
/*             contains the state-output cross-covariance matrix to be */
/*             used as cross-weighting matrix when computing the Kalman */
/*             gain. */
/*             This parameter is not referenced if JOBCV = 'N'. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S. */
/*             LDS >= N,  if JOBCV = 'C'; */
/*             LDS >= 1,  if JOBCV = 'N'. */

/*     O       (output) DOUBLE PRECISION array, dimension ( LDO,N ) */
/*             If  METH = 'M'  and  JOBCV = 'C',  or  METH = 'N', */
/*             the leading  L*NOBR-by-N  part of this array contains */
/*             the estimated extended observability matrix, i.e., the */
/*             first  N  columns of the relevant singular vectors. */
/*             If  METH = 'M'  and  JOBCV = 'N',  this array is not */
/*             referenced. */

/*     LDO     INTEGER */
/*             The leading dimension of the array  O. */
/*             LDO >= L*NOBR,  if  JOBCV = 'C'  or  METH = 'N'; */
/*             LDO >= 1,       otherwise. */

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
/*             LIWORK = N,                   if METH = 'M' and M = 0 */
/*                                        or JOB = 'C' and JOBCV = 'N'; */
/*             LIWORK = M*NOBR+N,            if METH = 'M', JOB = 'C', */
/*                                           and JOBCV = 'C'; */
/*             LIWORK = max(L*NOBR,M*NOBR),  if METH = 'M', JOB <> 'C', */
/*                                           and JOBCV = 'N'; */
/*             LIWORK = max(L*NOBR,M*NOBR+N),  if METH = 'M', JOB <> 'C', */
/*                                             and JOBCV = 'C'; */
/*             LIWORK = max(M*NOBR+N,M*(N+L)), if METH = 'N'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK,  and  DWORK(2),  DWORK(3),  DWORK(4),  and */
/*             DWORK(5)  contain the reciprocal condition numbers of the */
/*             triangular factors of the matrices, defined in the code, */
/*             GaL  (GaL = Un(1:(s-1)*L,1:n)),  R_1c  (if  METH = 'M'), */
/*             M  (if  JOBCV = 'C'  or  METH = 'N'),  and  Q  or  T  (see */
/*             SLICOT Library routines IB01PY or IB01PX),  respectively. */
/*             If  METH = 'N',  DWORK(3)  is set to one without any */
/*             calculations. Similarly, if  METH = 'M'  and  JOBCV = 'N', */
/*             DWORK(4)  is set to one. If  M = 0  or  JOB = 'C', */
/*             DWORK(3)  and  DWORK(5)  are set to one. */
/*             On exit, if  INFO = -30,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( LDW1,LDW2 ), where, if METH = 'M', */
/*             LDW1 >= max( 2*(L*NOBR-L)*N+2*N, (L*NOBR-L)*N+N*N+7*N ), */
/*                     if JOB = 'C' or JOB = 'A' and M = 0; */
/*             LDW1 >= max( 2*(L*NOBR-L)*N+N*N+7*N, */
/*                          (L*NOBR-L)*N+N+6*M*NOBR, (L*NOBR-L)*N+N+ */
/*                          max( L+M*NOBR, L*NOBR + */
/*                                         max( 3*L*NOBR+1, M ) ) ) */
/*                     if M > 0 and JOB = 'A', 'B', or 'D'; */
/*             LDW2 >= 0,                                 if JOBCV = 'N'; */
/*             LDW2 >= max( (L*NOBR-L)*N+Aw+2*N+max(5*N,(2*M+L)*NOBR+L), */
/*                          4*(M*NOBR+N)+1, M*NOBR+2*N+L ), */
/*                                                        if JOBCV = 'C', */
/*             where Aw = N+N*N, if M = 0 or JOB = 'C'; */
/*                   Aw = 0,     otherwise; */
/*             and, if METH = 'N', */
/*             LDW1 >= max( (L*NOBR-L)*N+2*N+(2*M+L)*NOBR+L, */
/*                          2*(L*NOBR-L)*N+N*N+8*N, N+4*(M*NOBR+N)+1, */
/*                          M*NOBR+3*N+L ); */
/*             LDW2 >= 0, if M = 0 or JOB = 'C'; */
/*             LDW2 >= M*NOBR*(N+L)*(M*(N+L)+1)+ */
/*                     max( (N+L)**2, 4*M*(N+L)+1 ), */
/*                     if M > 0 and JOB = 'A', 'B', or 'D'. */
/*             For good performance,  LDWORK  should be larger. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 4:  a least squares problem to be solved has a */
/*                   rank-deficient coefficient matrix; */
/*             = 5:  the computed covariance matrices are too small. */
/*                   The problem seems to be a deterministic one. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 2:  the singular value decomposition (SVD) algorithm did */
/*                   not converge; */
/*             = 3:  a singular upper triangular matrix was found. */

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

/*     REFERENCES */

/*     [1] Verhaegen M., and Dewilde, P. */
/*         Subspace Model Identification. Part 1: The output-error state- */
/*         space model identification class of algorithms. */
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

/*     The implemented method is numerically stable. */

/*     FURTHER COMMENTS */

/*     In some applications, it is useful to compute the system matrices */
/*     using two calls to this routine, the first one with  JOB = 'C', */
/*     and the second one with  JOB = 'B' or 'D'.  This is slightly less */
/*     efficient than using a single call with  JOB = 'A',  because some */
/*     calculations are repeated. If  METH = 'N',  all the calculations */
/*     at the first call are performed again at the second call; */
/*     moreover, it is required to save the needed submatrices of  R */
/*     before the first call and restore them before the second call. */
/*     If the covariance matrices are desired,  JOBCV  should be set */
/*     to  'C'  at the second call. If  B  and  D  are both needed, they */
/*     should be computed at once. */
/*     It is possible to compute the matrices A and C using the MOESP */
/*     algorithm (METH = 'M'), and the matrices B and D using the N4SID */
/*     algorithm (METH = 'N'). This combination could be slightly more */
/*     efficient than N4SID algorithm alone and it could be more accurate */
/*     than MOESP algorithm. No saving/restoring is needed in such a */
/*     combination, provided  JOBCV  is set to  'N'  at the first call. */
/*     Recommended usage:  either one call with  JOB = 'A',  or */
/*        first  call with  METH = 'M',  JOB = 'C',  JOBCV = 'N', */
/*        second call with  METH = 'M',  JOB = 'D',  JOBCV = 'C',  or */
/*        first  call with  METH = 'M',  JOB = 'C',  JOBCV = 'N', */
/*        second call with  METH = 'N',  JOB = 'D',  JOBCV = 'C'. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 1999. */

/*     REVISIONS */

/*     March 2000, Feb. 2001, Sep. 2001, March 2005. */

/*     KEYWORDS */

/*     Identification methods; least squares solutions; multivariable */
/*     systems; QR decomposition; singular value decomposition. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. Local Array .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

/*     Decode the scalar input parameters. */

#line 422 "IB01PD.f"
    /* Parameter adjustments */
#line 422 "IB01PD.f"
    r_dim1 = *ldr;
#line 422 "IB01PD.f"
    r_offset = 1 + r_dim1;
#line 422 "IB01PD.f"
    r__ -= r_offset;
#line 422 "IB01PD.f"
    a_dim1 = *lda;
#line 422 "IB01PD.f"
    a_offset = 1 + a_dim1;
#line 422 "IB01PD.f"
    a -= a_offset;
#line 422 "IB01PD.f"
    c_dim1 = *ldc;
#line 422 "IB01PD.f"
    c_offset = 1 + c_dim1;
#line 422 "IB01PD.f"
    c__ -= c_offset;
#line 422 "IB01PD.f"
    b_dim1 = *ldb;
#line 422 "IB01PD.f"
    b_offset = 1 + b_dim1;
#line 422 "IB01PD.f"
    b -= b_offset;
#line 422 "IB01PD.f"
    d_dim1 = *ldd;
#line 422 "IB01PD.f"
    d_offset = 1 + d_dim1;
#line 422 "IB01PD.f"
    d__ -= d_offset;
#line 422 "IB01PD.f"
    q_dim1 = *ldq;
#line 422 "IB01PD.f"
    q_offset = 1 + q_dim1;
#line 422 "IB01PD.f"
    q -= q_offset;
#line 422 "IB01PD.f"
    ry_dim1 = *ldry;
#line 422 "IB01PD.f"
    ry_offset = 1 + ry_dim1;
#line 422 "IB01PD.f"
    ry -= ry_offset;
#line 422 "IB01PD.f"
    s_dim1 = *lds;
#line 422 "IB01PD.f"
    s_offset = 1 + s_dim1;
#line 422 "IB01PD.f"
    s -= s_offset;
#line 422 "IB01PD.f"
    o_dim1 = *ldo;
#line 422 "IB01PD.f"
    o_offset = 1 + o_dim1;
#line 422 "IB01PD.f"
    o -= o_offset;
#line 422 "IB01PD.f"
    --iwork;
#line 422 "IB01PD.f"
    --dwork;
#line 422 "IB01PD.f"

#line 422 "IB01PD.f"
    /* Function Body */
#line 422 "IB01PD.f"
    moesp = lsame_(meth, "M", (ftnlen)1, (ftnlen)1);
#line 423 "IB01PD.f"
    n4sid = lsame_(meth, "N", (ftnlen)1, (ftnlen)1);
#line 424 "IB01PD.f"
    withal = lsame_(job, "A", (ftnlen)1, (ftnlen)1);
#line 425 "IB01PD.f"
    withc = lsame_(job, "C", (ftnlen)1, (ftnlen)1) || withal;
#line 426 "IB01PD.f"
    withd = lsame_(job, "D", (ftnlen)1, (ftnlen)1) || withal;
#line 427 "IB01PD.f"
    withb = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || withd;
#line 428 "IB01PD.f"
    withco = lsame_(jobcv, "C", (ftnlen)1, (ftnlen)1);
#line 429 "IB01PD.f"
    mnobr = *m * *nobr;
#line 430 "IB01PD.f"
    lnobr = *l * *nobr;
#line 431 "IB01PD.f"
    lmnobr = lnobr + mnobr;
#line 432 "IB01PD.f"
    lmmnob = lnobr + (mnobr << 1);
#line 433 "IB01PD.f"
    mnobrn = mnobr + *n;
#line 434 "IB01PD.f"
    lnobrn = lnobr - *n;
#line 435 "IB01PD.f"
    ldun2 = lnobr - *l;
#line 436 "IB01PD.f"
    ldunn = ldun2 * *n;
#line 437 "IB01PD.f"
    lmmnol = lmmnob + *l;
#line 438 "IB01PD.f"
    nr = lmnobr + lmnobr;
#line 439 "IB01PD.f"
    npl = *n + *l;
#line 440 "IB01PD.f"
    n2 = *n + *n;
#line 441 "IB01PD.f"
    nn = *n * *n;
#line 442 "IB01PD.f"
    minwrk = 1;
#line 443 "IB01PD.f"
    *iwarn = 0;
#line 444 "IB01PD.f"
    *info = 0;

/*     Check the scalar input parameters. */

#line 448 "IB01PD.f"
    if (! (moesp || n4sid)) {
#line 449 "IB01PD.f"
	*info = -1;
#line 450 "IB01PD.f"
    } else if (! (withb || withc)) {
#line 451 "IB01PD.f"
	*info = -2;
#line 452 "IB01PD.f"
    } else if (! (withco || lsame_(jobcv, "N", (ftnlen)1, (ftnlen)1))) {
#line 453 "IB01PD.f"
	*info = -3;
#line 454 "IB01PD.f"
    } else if (*nobr <= 1) {
#line 455 "IB01PD.f"
	*info = -4;
#line 456 "IB01PD.f"
    } else if (*n <= 0 || *n >= *nobr) {
#line 457 "IB01PD.f"
	*info = -5;
#line 458 "IB01PD.f"
    } else if (*m < 0) {
#line 459 "IB01PD.f"
	*info = -6;
#line 460 "IB01PD.f"
    } else if (*l <= 0) {
#line 461 "IB01PD.f"
	*info = -7;
#line 462 "IB01PD.f"
    } else if (withco && *nsmpl < nr) {
#line 463 "IB01PD.f"
	*info = -8;
#line 464 "IB01PD.f"
    } else if (*ldr < nr) {
#line 465 "IB01PD.f"
	*info = -10;
#line 466 "IB01PD.f"
    } else if (*lda < 1 || (withc || withb && n4sid) && *lda < *n) {
#line 468 "IB01PD.f"
	*info = -12;
#line 469 "IB01PD.f"
    } else if (*ldc < 1 || (withc || withb && n4sid) && *ldc < *l) {
#line 471 "IB01PD.f"
	*info = -14;
#line 472 "IB01PD.f"
    } else if (*ldb < 1 || withb && *ldb < *n && *m > 0) {
#line 474 "IB01PD.f"
	*info = -16;
#line 475 "IB01PD.f"
    } else if (*ldd < 1 || withd && *ldd < *l && *m > 0) {
#line 477 "IB01PD.f"
	*info = -18;
#line 478 "IB01PD.f"
    } else if (*ldq < 1 || withco && *ldq < *n) {
#line 479 "IB01PD.f"
	*info = -20;
#line 480 "IB01PD.f"
    } else if (*ldry < 1 || withco && *ldry < *l) {
#line 481 "IB01PD.f"
	*info = -22;
#line 482 "IB01PD.f"
    } else if (*lds < 1 || withco && *lds < *n) {
#line 483 "IB01PD.f"
	*info = -24;
#line 484 "IB01PD.f"
    } else if (*ldo < 1 || (withco || n4sid) && *ldo < lnobr) {
#line 486 "IB01PD.f"
	*info = -26;
#line 487 "IB01PD.f"
    } else {

/*        Compute workspace. */
/*        (Note: Comments in the code beginning "Workspace:" describe the */
/*         minimal amount of workspace needed at that point in the code, */
/*         as well as the preferred amount for good performance. */
/*         NB refers to the optimal block size for the immediately */
/*         following subroutine, as returned by ILAENV.) */

#line 496 "IB01PD.f"
	iaw = 0;
#line 497 "IB01PD.f"
	minwrk = ldunn + (*n << 2);
#line 498 "IB01PD.f"
	maxwrk = ldunn + *n + *n * ilaenv_(&c__1, "DGEQRF", " ", &ldun2, n, &
		c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 500 "IB01PD.f"
	if (moesp) {
#line 501 "IB01PD.f"
	    id = 0;
#line 502 "IB01PD.f"
	    if (withc) {
/* Computing MAX */
#line 503 "IB01PD.f"
		i__1 = minwrk, i__2 = (ldunn << 1) + n2, i__1 = max(i__1,i__2)
			, i__2 = ldunn + nn + *n * 7;
#line 503 "IB01PD.f"
		minwrk = max(i__1,i__2);
/* Computing MAX */
#line 504 "IB01PD.f"
		i__1 = maxwrk, i__2 = (ldunn << 1) + *n + *n * ilaenv_(&c__1, 
			"DORMQR", "LT", &ldun2, n, n, &c_n1, (ftnlen)6, (
			ftnlen)2);
#line 504 "IB01PD.f"
		maxwrk = max(i__1,i__2);
#line 506 "IB01PD.f"
	    }
#line 507 "IB01PD.f"
	} else {
#line 508 "IB01PD.f"
	    id = *n;
#line 509 "IB01PD.f"
	}

#line 511 "IB01PD.f"
	if (*m > 0 && withb || n4sid) {
/* Computing MAX */
#line 512 "IB01PD.f"
	    i__1 = minwrk, i__2 = (ldunn << 1) + nn + id + *n * 7;
#line 512 "IB01PD.f"
	    minwrk = max(i__1,i__2);
#line 513 "IB01PD.f"
	    if (moesp) {
/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 513 "IB01PD.f"
		i__5 = lnobr * 3 + 1;
#line 513 "IB01PD.f"
		i__3 = *l + mnobr, i__4 = lnobr + max(i__5,*m);
#line 513 "IB01PD.f"
		i__1 = minwrk, i__2 = ldunn + *n + mnobr * 6, i__1 = max(i__1,
			i__2), i__2 = ldunn + *n + max(i__3,i__4);
#line 513 "IB01PD.f"
		minwrk = max(i__1,i__2);
#line 513 "IB01PD.f"
	    }
#line 517 "IB01PD.f"
	} else {
#line 518 "IB01PD.f"
	    if (moesp) {
#line 518 "IB01PD.f"
		iaw = *n + nn;
#line 518 "IB01PD.f"
	    }
#line 520 "IB01PD.f"
	}

#line 522 "IB01PD.f"
	if (n4sid || withco) {
/* Computing MAX */
/* Computing MAX */
#line 523 "IB01PD.f"
	    i__3 = *n * 5;
#line 523 "IB01PD.f"
	    i__1 = minwrk, i__2 = ldunn + iaw + n2 + max(i__3,lmmnol), i__1 = 
		    max(i__1,i__2), i__2 = id + (mnobrn << 2) + 1, i__1 = max(
		    i__1,i__2), i__2 = id + mnobrn + npl;
#line 523 "IB01PD.f"
	    minwrk = max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 525 "IB01PD.f"
	    i__3 = *n * ilaenv_(&c__1, "DGEQRF", " ", &lnobr, n, &c_n1, &c_n1,
		     (ftnlen)6, (ftnlen)1), i__4 = lmmnob * ilaenv_(&c__1, 
		    "DORMQR", "LT", &lnobr, &lmmnob, n, &c_n1, (ftnlen)6, (
		    ftnlen)2), i__3 = max(i__3,i__4), i__4 = lmmnol * ilaenv_(
		    &c__1, "DORMQR", "LT", &ldun2, &lmmnol, n, &c_n1, (ftnlen)
		    6, (ftnlen)2);
#line 525 "IB01PD.f"
	    i__1 = maxwrk, i__2 = ldunn + iaw + n2 + max(i__3,i__4), i__1 = 
		    max(i__1,i__2), i__2 = id + *n + *n * ilaenv_(&c__1, 
		    "DGEQRF", " ", &lmnobr, n, &c_n1, &c_n1, (ftnlen)6, (
		    ftnlen)1), i__1 = max(i__1,i__2), i__2 = id + *n + npl * 
		    ilaenv_(&c__1, "DORMQR", "LT", &lmnobr, &npl, n, &c_n1, (
		    ftnlen)6, (ftnlen)2);
#line 525 "IB01PD.f"
	    maxwrk = max(i__1,i__2);
#line 536 "IB01PD.f"
	    if (n4sid && (*m > 0 && withb)) {
/* Computing MAX */
/* Computing MAX */
/* Computing 2nd power */
#line 536 "IB01PD.f"
		i__5 = npl;
#line 536 "IB01PD.f"
		i__3 = i__5 * i__5, i__4 = (*m << 2) * npl + 1;
#line 536 "IB01PD.f"
		i__1 = minwrk, i__2 = mnobr * npl * (*m * npl + 1) + max(i__3,
			i__4);
#line 536 "IB01PD.f"
		minwrk = max(i__1,i__2);
#line 536 "IB01PD.f"
	    }
#line 539 "IB01PD.f"
	}
#line 540 "IB01PD.f"
	maxwrk = max(minwrk,maxwrk);

#line 542 "IB01PD.f"
	if (*ldwork < minwrk) {
#line 543 "IB01PD.f"
	    *info = -30;
#line 544 "IB01PD.f"
	    dwork[1] = (doublereal) minwrk;
#line 545 "IB01PD.f"
	}
#line 546 "IB01PD.f"
    }

/*     Return if there are illegal arguments. */

#line 550 "IB01PD.f"
    if (*info != 0) {
#line 551 "IB01PD.f"
	i__1 = -(*info);
#line 551 "IB01PD.f"
	xerbla_("IB01PD", &i__1, (ftnlen)6);
#line 552 "IB01PD.f"
	return 0;
#line 553 "IB01PD.f"
    }

#line 555 "IB01PD.f"
    nr2 = mnobr + 1;
#line 556 "IB01PD.f"
    nr3 = lmnobr + 1;
#line 557 "IB01PD.f"
    nr4 = lmmnob + 1;

/*     Set the precision parameters. A threshold value  EPS**(2/3)  is */
/*     used for deciding to use pivoting or not, where  EPS  is the */
/*     relative machine precision (see LAPACK Library routine DLAMCH). */

#line 563 "IB01PD.f"
    eps = dlamch_("Precision", (ftnlen)9);
#line 564 "IB01PD.f"
    thresh = pow_dd(&eps, &c_b43);
#line 565 "IB01PD.f"
    svlmax = 0.;
#line 566 "IB01PD.f"
    rcond4 = 1.;

/*     Let  Un  be the matrix of left singular vectors (stored in  R_22). */
/*     Copy  un1 = GaL = Un(1:(s-1)*L,1:n)  in the workspace. */

#line 571 "IB01PD.f"
    igal = 1;
#line 572 "IB01PD.f"
    dlacpy_("Full", &ldun2, n, &r__[nr2 + nr2 * r_dim1], ldr, &dwork[igal], &
	    ldun2, (ftnlen)4);

/*     Factor un1 = Q1*[r1'  0]' (' means transposition). */
/*     Workspace: need   L*(NOBR-1)*N+2*N, */
/*                prefer L*(NOBR-1)*N+N+N*NB. */

#line 579 "IB01PD.f"
    itau1 = igal + ldunn;
#line 580 "IB01PD.f"
    jwork = itau1 + *n;
#line 581 "IB01PD.f"
    ldw = jwork;
#line 582 "IB01PD.f"
    i__1 = *ldwork - jwork + 1;
#line 582 "IB01PD.f"
    dgeqrf_(&ldun2, n, &dwork[igal], &ldun2, &dwork[itau1], &dwork[jwork], &
	    i__1, &ierr);

/*     Compute the reciprocal of the condition number of r1. */
/*     Workspace: need L*(NOBR-1)*N+4*N. */

#line 588 "IB01PD.f"
    dtrcon_("1-norm", "Upper", "NonUnit", n, &dwork[igal], &ldun2, &rcond1, &
	    dwork[jwork], &iwork[1], info, (ftnlen)6, (ftnlen)5, (ftnlen)7);

#line 591 "IB01PD.f"
    toll1 = *tol;
#line 592 "IB01PD.f"
    if (toll1 <= 0.) {
#line 592 "IB01PD.f"
	toll1 = nn * eps;
#line 592 "IB01PD.f"
    }

#line 595 "IB01PD.f"
    if (*m > 0 && withb || n4sid) {
#line 596 "IB01PD.f"
	*(unsigned char *)jobp = 'P';
#line 597 "IB01PD.f"
	if (withal) {
#line 598 "IB01PD.f"
	    *(unsigned char *)jobpy = 'D';
#line 599 "IB01PD.f"
	} else {
#line 600 "IB01PD.f"
	    *(unsigned char *)jobpy = *(unsigned char *)job;
#line 601 "IB01PD.f"
	}
#line 602 "IB01PD.f"
    } else {
#line 603 "IB01PD.f"
	*(unsigned char *)jobp = 'N';
#line 604 "IB01PD.f"
    }

#line 606 "IB01PD.f"
    if (moesp) {
#line 607 "IB01PD.f"
	ncol = 0;
#line 608 "IB01PD.f"
	iun2 = jwork;
#line 609 "IB01PD.f"
	if (withc) {

/*           Set  C = Un(1:L,1:n)  and then compute the system matrix A. */

/*           Set  un2 = Un(L+1:L*s,1:n)  in  DWORK(IUN2). */
/*           Workspace: need   2*L*(NOBR-1)*N+N. */

#line 616 "IB01PD.f"
	    dlacpy_("Full", l, n, &r__[nr2 + nr2 * r_dim1], ldr, &c__[
		    c_offset], ldc, (ftnlen)4);
#line 617 "IB01PD.f"
	    dlacpy_("Full", &ldun2, n, &r__[nr2 + *l + nr2 * r_dim1], ldr, &
		    dwork[iun2], &ldun2, (ftnlen)4);

/*           Note that un1 has already been factored as */
/*           un1 = Q1*[r1'  0]'  and usually (generically, assuming */
/*           observability) has full column rank. */
/*           Update  un2 <-- Q1'*un2  in  DWORK(IUN2)  and save its */
/*           first  n  rows in  A. */
/*           Workspace: need   2*L*(NOBR-1)*N+2*N; */
/*                      prefer 2*L*(NOBR-1)*N+N+N*NB. */

#line 628 "IB01PD.f"
	    jwork = iun2 + ldunn;
#line 629 "IB01PD.f"
	    i__1 = *ldwork - jwork + 1;
#line 629 "IB01PD.f"
	    dormqr_("Left", "Transpose", &ldun2, n, n, &dwork[igal], &ldun2, &
		    dwork[itau1], &dwork[iun2], &ldun2, &dwork[jwork], &i__1, 
		    &ierr, (ftnlen)4, (ftnlen)9);
#line 632 "IB01PD.f"
	    dlacpy_("Full", n, n, &dwork[iun2], &ldun2, &a[a_offset], lda, (
		    ftnlen)4);
#line 633 "IB01PD.f"
	    ncol = *n;
#line 634 "IB01PD.f"
	    jwork = iun2;
#line 635 "IB01PD.f"
	}

#line 637 "IB01PD.f"
	if (rcond1 > max(toll1,thresh)) {

/*           The triangular factor r1 is considered to be of full rank. */
/*           Solve for  A  (if requested),  r1*A = un2(1:n,:)  in  A. */

#line 642 "IB01PD.f"
	    if (withc) {
#line 643 "IB01PD.f"
		dtrtrs_("Upper", "NoTranspose", "NonUnit", n, n, &dwork[igal],
			 &ldun2, &a[a_offset], lda, &ierr, (ftnlen)5, (ftnlen)
			11, (ftnlen)7);
#line 645 "IB01PD.f"
		if (ierr > 0) {
#line 646 "IB01PD.f"
		    *info = 3;
#line 647 "IB01PD.f"
		    return 0;
#line 648 "IB01PD.f"
		}
#line 649 "IB01PD.f"
	    }
#line 650 "IB01PD.f"
	    rank = *n;
#line 651 "IB01PD.f"
	} else {

/*           Rank-deficient triangular factor r1.  Use SVD of r1, */
/*           r1 = U*S*V',  also for computing  A  (if requested) from */
/*           r1*A = un2(1:n,:).  Matrix  U  is computed in  DWORK(IU), */
/*           and  V' overwrites  r1.  If  B  is requested, the */
/*           pseudoinverse of  r1  and then of  GaL  are also computed */
/*           in  R(NR3,NR2). */
/*           Workspace: need   c*L*(NOBR-1)*N+N*N+7*N, */
/*                             where  c = 1  if  B and D  are not needed, */
/*                                    c = 2  if  B and D  are needed; */
/*                      prefer larger. */

#line 664 "IB01PD.f"
	    iu = iun2;
#line 665 "IB01PD.f"
	    isv = iu + nn;
#line 666 "IB01PD.f"
	    jwork = isv + *n;
#line 667 "IB01PD.f"
	    if (*m > 0 && withb) {

/*              Save the elementary reflectors used for computing r1, */
/*              if  B, D  are needed. */
/*              Workspace: need   2*L*(NOBR-1)*N+2*N+N*N. */

#line 673 "IB01PD.f"
		ihous = jwork;
#line 674 "IB01PD.f"
		jwork = ihous + ldunn;
#line 675 "IB01PD.f"
		dlacpy_("Lower", &ldun2, n, &dwork[igal], &ldun2, &dwork[
			ihous], &ldun2, (ftnlen)5);
#line 677 "IB01PD.f"
	    } else {
#line 678 "IB01PD.f"
		ihous = igal;
#line 679 "IB01PD.f"
	    }

#line 681 "IB01PD.f"
	    i__1 = *ldwork - jwork + 1;
#line 681 "IB01PD.f"
	    mb02ud_("Not factored", "Left", "NoTranspose", jobp, n, &ncol, &
		    c_b63, &toll1, &rank, &dwork[igal], &ldun2, &dwork[iu], n,
		     &dwork[isv], &a[a_offset], lda, &r__[nr3 + nr2 * r_dim1],
		     ldr, &dwork[jwork], &i__1, &ierr, (ftnlen)12, (ftnlen)4, 
		    (ftnlen)11, (ftnlen)1);
#line 685 "IB01PD.f"
	    if (ierr != 0) {
#line 686 "IB01PD.f"
		*info = 2;
#line 687 "IB01PD.f"
		return 0;
#line 688 "IB01PD.f"
	    }
/* Computing MAX */
#line 689 "IB01PD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 689 "IB01PD.f"
	    maxwrk = max(i__1,i__2);

#line 691 "IB01PD.f"
	    if (rank == 0) {
#line 692 "IB01PD.f"
		*(unsigned char *)jobp = 'N';
#line 693 "IB01PD.f"
	    } else if (*m > 0 && withb) {

/*              Compute  pinv(GaL)  in  R(NR3,NR2)  if  B, D  are needed. */
/*              Workspace: need   2*L*(NOBR-1)*N+N*N+3*N; */
/*                         prefer 2*L*(NOBR-1)*N+N*N+2*N+N*NB. */

#line 699 "IB01PD.f"
		i__1 = ldun2 - *n;
#line 699 "IB01PD.f"
		dlaset_("Full", n, &i__1, &c_b66, &c_b66, &r__[nr3 + (nr2 + *
			n) * r_dim1], ldr, (ftnlen)4);
#line 701 "IB01PD.f"
		i__1 = *ldwork - jwork + 1;
#line 701 "IB01PD.f"
		dormqr_("Right", "Transpose", n, &ldun2, n, &dwork[ihous], &
			ldun2, &dwork[itau1], &r__[nr3 + nr2 * r_dim1], ldr, &
			dwork[jwork], &i__1, &ierr, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 705 "IB01PD.f"
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 705 "IB01PD.f"
		maxwrk = max(i__1,i__2);
#line 706 "IB01PD.f"
		if (withco) {

/*                 Save  pinv(GaL)  in  DWORK(IGAL). */

#line 710 "IB01PD.f"
		    dlacpy_("Full", n, &ldun2, &r__[nr3 + nr2 * r_dim1], ldr, 
			    &dwork[igal], n, (ftnlen)4);
#line 712 "IB01PD.f"
		}
#line 713 "IB01PD.f"
		jwork = iun2;
#line 714 "IB01PD.f"
	    }
#line 715 "IB01PD.f"
	    ldw = jwork;
#line 716 "IB01PD.f"
	}

#line 718 "IB01PD.f"
	if (*m > 0 && withb) {

/*           Computation of  B  and  D. */

/*           Compute the reciprocal of the condition number of R_1c. */
/*           Workspace: need L*(NOBR-1)*N+N+3*M*NOBR. */

#line 725 "IB01PD.f"
	    dtrcon_("1-norm", "Upper", "NonUnit", &mnobr, &r__[nr3 + r_dim1], 
		    ldr, &rcond2, &dwork[jwork], &iwork[1], &ierr, (ftnlen)6, 
		    (ftnlen)5, (ftnlen)7);

#line 728 "IB01PD.f"
	    toll = *tol;
#line 729 "IB01PD.f"
	    if (toll <= 0.) {
#line 729 "IB01PD.f"
		toll = mnobr * mnobr * eps;
#line 729 "IB01PD.f"
	    }

/*           Compute the right hand side and solve for  K  (in  R_23), */
/*              K*R_1c' = u2'*R_2c', */
/*           where  u2 = Un(:,n+1:L*s),  and  K  is  (Ls-n) x ms. */

#line 736 "IB01PD.f"
	    dgemm_("Transpose", "Transpose", &lnobrn, &mnobr, &lnobr, &c_b63, 
		    &r__[nr2 + (nr2 + *n) * r_dim1], ldr, &r__[nr2 * r_dim1 + 
		    1], ldr, &c_b66, &r__[nr2 + nr3 * r_dim1], ldr, (ftnlen)9,
		     (ftnlen)9);

#line 740 "IB01PD.f"
	    if (rcond2 > max(toll,thresh)) {

/*              The triangular factor R_1c is considered to be of full */
/*              rank. Solve for  K,  K*R_1c' = u2'*R_2c'. */

#line 745 "IB01PD.f"
		dtrsm_("Right", "Upper", "Transpose", "Non-unit", &lnobrn, &
			mnobr, &c_b63, &r__[nr3 + r_dim1], ldr, &r__[nr2 + 
			nr3 * r_dim1], ldr, (ftnlen)5, (ftnlen)5, (ftnlen)9, (
			ftnlen)8);
#line 748 "IB01PD.f"
	    } else {

/*              Rank-deficient triangular factor  R_1c.  Use SVD of  R_1c */
/*              for computing  K  from  K*R_1c' = u2'*R_2c',  where */
/*              R_1c = U1*S1*V1'.  Matrix  U1  is computed in  R_33, */
/*              and  V1'  overwrites  R_1c. */
/*              Workspace: need   L*(NOBR-1)*N+N+6*M*NOBR; */
/*                         prefer larger. */

#line 757 "IB01PD.f"
		isv = ldw;
#line 758 "IB01PD.f"
		jwork = isv + mnobr;
#line 759 "IB01PD.f"
		i__1 = *ldwork - jwork + 1;
#line 759 "IB01PD.f"
		mb02ud_("Not factored", "Right", "Transpose", "No pinv", &
			lnobrn, &mnobr, &c_b63, &toll, &rank11, &r__[nr3 + 
			r_dim1], ldr, &r__[nr3 + nr3 * r_dim1], ldr, &dwork[
			isv], &r__[nr2 + nr3 * r_dim1], ldr, &dwork[jwork], &
			c__1, &dwork[jwork], &i__1, &ierr, (ftnlen)12, (
			ftnlen)5, (ftnlen)9, (ftnlen)7);
#line 764 "IB01PD.f"
		if (ierr != 0) {
#line 765 "IB01PD.f"
		    *info = 2;
#line 766 "IB01PD.f"
		    return 0;
#line 767 "IB01PD.f"
		}
/* Computing MAX */
#line 768 "IB01PD.f"
		i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 768 "IB01PD.f"
		maxwrk = max(i__1,i__2);
#line 769 "IB01PD.f"
		jwork = ldw;
#line 770 "IB01PD.f"
	    }

/*           Compute the triangular factor of the structured matrix  Q */
/*           and apply the transformations to the matrix  Kexpand,  where */
/*           Q  and  Kexpand  are defined in SLICOT Library routine */
/*           IB01PY.  Compute also the matrices  B,  D. */
/*           Workspace: need   L*(NOBR-1)*N+N+max(L+M*NOBR,L*NOBR+ */
/*                                                max(3*L*NOBR+1,M)); */
/*                      prefer larger. */

#line 780 "IB01PD.f"
	    if (withco) {
#line 780 "IB01PD.f"
		dlacpy_("Full", &lnobr, n, &r__[nr2 + nr2 * r_dim1], ldr, &o[
			o_offset], ldo, (ftnlen)4);
#line 780 "IB01PD.f"
	    }
#line 782 "IB01PD.f"
	    i__1 = *ldwork - jwork + 1;
#line 782 "IB01PD.f"
	    ib01py_(meth, jobpy, nobr, n, m, l, &rank, &r__[nr2 + nr2 * 
		    r_dim1], ldr, &dwork[igal], &ldun2, &dwork[itau1], &r__[
		    nr3 + nr2 * r_dim1], ldr, &r__[nr2 + nr3 * r_dim1], ldr, &
		    r__[nr4 + nr2 * r_dim1], ldr, &r__[nr4 + nr3 * r_dim1], 
		    ldr, &b[b_offset], ldb, &d__[d_offset], ldd, tol, &iwork[
		    1], &dwork[jwork], &i__1, iwarn, info, (ftnlen)1, (ftnlen)
		    1);
#line 788 "IB01PD.f"
	    if (*info != 0) {
#line 788 "IB01PD.f"
		return 0;
#line 788 "IB01PD.f"
	    }
/* Computing MAX */
#line 790 "IB01PD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 790 "IB01PD.f"
	    maxwrk = max(i__1,i__2);
#line 791 "IB01PD.f"
	    rcond4 = dwork[jwork + 1];
#line 792 "IB01PD.f"
	    if (withco) {
#line 792 "IB01PD.f"
		dlacpy_("Full", &lnobr, n, &o[o_offset], ldo, &r__[nr2 + 
			r_dim1], ldr, (ftnlen)4);
#line 792 "IB01PD.f"
	    }

#line 795 "IB01PD.f"
	} else {
#line 796 "IB01PD.f"
	    rcond2 = 1.;
#line 797 "IB01PD.f"
	}

#line 799 "IB01PD.f"
	if (! withco) {
#line 800 "IB01PD.f"
	    rcond3 = 1.;
#line 801 "IB01PD.f"
	    goto L30;
#line 802 "IB01PD.f"
	}
#line 803 "IB01PD.f"
    } else {

/*        For N4SID, set  RCOND2  to one. */

#line 807 "IB01PD.f"
	rcond2 = 1.;
#line 808 "IB01PD.f"
    }

/*     If needed, save the first  n  columns, representing  Gam,  of the */
/*     matrix of left singular vectors,  Un,  in  R_21  and in  O. */

#line 813 "IB01PD.f"
    if (n4sid || withc && ! withal) {
#line 814 "IB01PD.f"
	if (*m > 0) {
#line 814 "IB01PD.f"
	    dlacpy_("Full", &lnobr, n, &r__[nr2 + nr2 * r_dim1], ldr, &r__[
		    nr2 + r_dim1], ldr, (ftnlen)4);
#line 814 "IB01PD.f"
	}
#line 817 "IB01PD.f"
	dlacpy_("Full", &lnobr, n, &r__[nr2 + nr2 * r_dim1], ldr, &o[o_offset]
		, ldo, (ftnlen)4);
#line 818 "IB01PD.f"
    }

/*     Computations for covariance matrices, and system matrices (N4SID). */
/*     Solve the least squares problems  Gam*Y = R4(1:L*s,1:(2*m+L)*s), */
/*                                       GaL*X = R4(L+1:L*s,:),  where */
/*     GaL = Gam(1:L*(s-1),:),  Gam  has full column rank, and */
/*     R4 = [ R_14' R_24' R_34' R_44L' ],  R_44L = R_44(1:L,:), as */
/*     returned by SLICOT Library routine  IB01ND. */
/*     First, find the  QR  factorization of  Gam,  Gam = Q*R. */
/*     Workspace: need   L*(NOBR-1)*N+Aw+3*N; */
/*                prefer L*(NOBR-1)*N+Aw+2*N+N*NB, where */
/*                Aw = N+N*N,  if  (M = 0  or  JOB = 'C'),  rank(r1) < N, */
/*                             and  METH = 'M'; */
/*                Aw = 0,      otherwise. */

#line 833 "IB01PD.f"
    itau2 = ldw;
#line 834 "IB01PD.f"
    jwork = itau2 + *n;
#line 835 "IB01PD.f"
    i__1 = *ldwork - jwork + 1;
#line 835 "IB01PD.f"
    dgeqrf_(&lnobr, n, &r__[nr2 + r_dim1], ldr, &dwork[itau2], &dwork[jwork], 
	    &i__1, &ierr);

/*     For METH = 'M' or when JOB = 'B' or 'D', transpose */
/*     [ R_14' R_24' R_34' ]'  in the last block-row of  R, obtaining  Z, */
/*     and for METH = 'N' and JOB = 'A' or 'C', use the matrix  Z */
/*     already available in the last block-row of  R,  and then apply */
/*     the transformations, Z <-- Q'*Z. */
/*     Workspace: need   L*(NOBR-1)*N+Aw+2*N+(2*M+L)*NOBR; */
/*                prefer L*(NOBR-1)*N+Aw+2*N+(2*M+L)*NOBR*NB. */

#line 846 "IB01PD.f"
    if (moesp || withb && ! withal) {
#line 846 "IB01PD.f"
	ma02ad_("Full", &lmmnob, &lnobr, &r__[nr4 * r_dim1 + 1], ldr, &r__[
		nr4 + r_dim1], ldr, (ftnlen)4);
#line 846 "IB01PD.f"
    }
#line 849 "IB01PD.f"
    i__1 = *ldwork - jwork + 1;
#line 849 "IB01PD.f"
    dormqr_("Left", "Transpose", &lnobr, &lmmnob, n, &r__[nr2 + r_dim1], ldr, 
	    &dwork[itau2], &r__[nr4 + r_dim1], ldr, &dwork[jwork], &i__1, &
	    ierr, (ftnlen)4, (ftnlen)9);

/*     Solve for  Y,  RY = Z  in  Z  and save the transpose of the */
/*     solution  Y  in the second block-column of  R. */

#line 856 "IB01PD.f"
    dtrtrs_("Upper", "NoTranspose", "NonUnit", n, &lmmnob, &r__[nr2 + r_dim1],
	     ldr, &r__[nr4 + r_dim1], ldr, &ierr, (ftnlen)5, (ftnlen)11, (
	    ftnlen)7);
#line 858 "IB01PD.f"
    if (ierr > 0) {
#line 859 "IB01PD.f"
	*info = 3;
#line 860 "IB01PD.f"
	return 0;
#line 861 "IB01PD.f"
    }
#line 862 "IB01PD.f"
    ma02ad_("Full", n, &lmmnob, &r__[nr4 + r_dim1], ldr, &r__[nr2 * r_dim1 + 
	    1], ldr, (ftnlen)4);
#line 863 "IB01PD.f"
    nr4mn = nr4 - *n;
#line 864 "IB01PD.f"
    nr4pl = nr4 + *l;
#line 865 "IB01PD.f"
    nrow = lmmnol;

/*     SHIFT is .TRUE. if some columns of  R_14 : R_44L  should be */
/*     shifted to the right, to avoid overwriting useful information. */

#line 870 "IB01PD.f"
    shift = *m == 0 && lnobr < n2;

#line 872 "IB01PD.f"
    if (rcond1 > max(toll1,thresh)) {

/*        The triangular factor  r1  of  GaL  (GaL = Q1*r1)  is */
/*        considered to be of full rank. */

/*        Transpose  [ R_14' R_24' R_34' R_44L' ]'(:,L+1:L*s)  in the */
/*        last block-row of  R  (beginning with the  (L+1)-th  row), */
/*        obtaining  Z1,  and then apply the transformations, */
/*        Z1 <-- Q1'*Z1. */
/*        Workspace: need   L*(NOBR-1)*N+Aw+2*N+ (2*M+L)*NOBR + L; */
/*                   prefer L*(NOBR-1)*N+Aw+2*N+((2*M+L)*NOBR + L)*NB. */

#line 884 "IB01PD.f"
	ma02ad_("Full", &lmmnol, &ldun2, &r__[nr4pl * r_dim1 + 1], ldr, &r__[
		nr4pl + r_dim1], ldr, (ftnlen)4);
#line 886 "IB01PD.f"
	i__1 = *ldwork - jwork + 1;
#line 886 "IB01PD.f"
	dormqr_("Left", "Transpose", &ldun2, &lmmnol, n, &dwork[igal], &ldun2,
		 &dwork[itau1], &r__[nr4pl + r_dim1], ldr, &dwork[jwork], &
		i__1, &ierr, (ftnlen)4, (ftnlen)9);

/*        Solve for  X,  r1*X = Z1  in  Z1,  and copy the transpose of  X */
/*        into the last part of the third block-column of  R. */

#line 893 "IB01PD.f"
	dtrtrs_("Upper", "NoTranspose", "NonUnit", n, &lmmnol, &dwork[igal], &
		ldun2, &r__[nr4pl + r_dim1], ldr, &ierr, (ftnlen)5, (ftnlen)
		11, (ftnlen)7);
#line 895 "IB01PD.f"
	if (ierr > 0) {
#line 896 "IB01PD.f"
	    *info = 3;
#line 897 "IB01PD.f"
	    return 0;
#line 898 "IB01PD.f"
	}

#line 900 "IB01PD.f"
	if (shift) {
#line 901 "IB01PD.f"
	    nr4mn = nr4;

#line 903 "IB01PD.f"
	    for (i__ = *l - 1; i__ >= 0; --i__) {
#line 904 "IB01PD.f"
		dcopy_(&lmmnol, &r__[(nr4 + i__) * r_dim1 + 1], &c__1, &r__[(
			nr4 + *n + i__) * r_dim1 + 1], &c__1);
#line 905 "IB01PD.f"
/* L10: */
#line 905 "IB01PD.f"
	    }

#line 907 "IB01PD.f"
	}
#line 908 "IB01PD.f"
	ma02ad_("Full", n, &lmmnol, &r__[nr4pl + r_dim1], ldr, &r__[nr4mn * 
		r_dim1 + 1], ldr, (ftnlen)4);
#line 910 "IB01PD.f"
	nrow = 0;
#line 911 "IB01PD.f"
    }

#line 913 "IB01PD.f"
    if (n4sid || nrow > 0) {

/*        METH = 'N'  or rank-deficient triangular factor r1. */
/*        For  METH = 'N',  use SVD of  r1,  r1 = U*S*V', for computing */
/*        X'  from  X'*GaL' = Z1',  if  rank(r1) < N.  Matrix  U  is */
/*        computed in  DWORK(IU)  and  V'  overwrites  r1.  Then, the */
/*        pseudoinverse of  GaL  is determined in  R(NR4+L,NR2). */
/*        For METH = 'M', the pseudoinverse of  GaL  is already available */
/*        if  M > 0  and  B  is requested;  otherwise, the SVD of  r1  is */
/*        available in  DWORK(IU),  DWORK(ISV),  and  DWORK(IGAL). */
/*        Workspace for N4SID: need   2*L*(NOBR-1)*N+N*N+8*N; */
/*                             prefer larger. */

#line 926 "IB01PD.f"
	if (moesp) {
#line 927 "IB01PD.f"
	    *(unsigned char *)fact = 'F';
#line 928 "IB01PD.f"
	    if (*m > 0 && withb) {
#line 928 "IB01PD.f"
		dlacpy_("Full", n, &ldun2, &dwork[igal], n, &r__[nr4pl + nr2 *
			 r_dim1], ldr, (ftnlen)4);
#line 928 "IB01PD.f"
	    }
#line 931 "IB01PD.f"
	} else {

/*           Save the elementary reflectors used for computing r1. */

#line 935 "IB01PD.f"
	    ihous = jwork;
#line 936 "IB01PD.f"
	    dlacpy_("Lower", &ldun2, n, &dwork[igal], &ldun2, &dwork[ihous], &
		    ldun2, (ftnlen)5);
#line 938 "IB01PD.f"
	    *(unsigned char *)fact = 'N';
#line 939 "IB01PD.f"
	    iu = ihous + ldunn;
#line 940 "IB01PD.f"
	    isv = iu + nn;
#line 941 "IB01PD.f"
	    jwork = isv + *n;
#line 942 "IB01PD.f"
	}

#line 944 "IB01PD.f"
	i__1 = *ldwork - jwork + 1;
#line 944 "IB01PD.f"
	mb02ud_(fact, "Right", "Transpose", jobp, &nrow, n, &c_b63, &toll1, &
		rank, &dwork[igal], &ldun2, &dwork[iu], n, &dwork[isv], &r__[
		nr4pl * r_dim1 + 1], ldr, &r__[nr4pl + nr2 * r_dim1], ldr, &
		dwork[jwork], &i__1, &ierr, (ftnlen)1, (ftnlen)5, (ftnlen)9, (
		ftnlen)1);
#line 948 "IB01PD.f"
	if (nrow > 0) {
#line 949 "IB01PD.f"
	    if (shift) {
#line 950 "IB01PD.f"
		nr4mn = nr4;
#line 951 "IB01PD.f"
		dlacpy_("Full", &lmmnol, l, &r__[nr4 * r_dim1 + 1], ldr, &r__[
			(nr4 - *l) * r_dim1 + 1], ldr, (ftnlen)4);
#line 953 "IB01PD.f"
		dlacpy_("Full", &lmmnol, n, &r__[nr4pl * r_dim1 + 1], ldr, &
			r__[nr4mn * r_dim1 + 1], ldr, (ftnlen)4);
#line 955 "IB01PD.f"
		dlacpy_("Full", &lmmnol, l, &r__[(nr4 - *l) * r_dim1 + 1], 
			ldr, &r__[(nr4 + *n) * r_dim1 + 1], ldr, (ftnlen)4);
#line 957 "IB01PD.f"
	    } else {
#line 958 "IB01PD.f"
		dlacpy_("Full", &lmmnol, n, &r__[nr4pl * r_dim1 + 1], ldr, &
			r__[nr4mn * r_dim1 + 1], ldr, (ftnlen)4);
#line 960 "IB01PD.f"
	    }
#line 961 "IB01PD.f"
	}

#line 963 "IB01PD.f"
	if (n4sid) {
#line 964 "IB01PD.f"
	    if (ierr != 0) {
#line 965 "IB01PD.f"
		*info = 2;
#line 966 "IB01PD.f"
		return 0;
#line 967 "IB01PD.f"
	    }
/* Computing MAX */
#line 968 "IB01PD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 968 "IB01PD.f"
	    maxwrk = max(i__1,i__2);

/*           Compute  pinv(GaL)  in  R(NR4+L,NR2). */
/*           Workspace: need   2*L*(NOBR-1)*N+3*N; */
/*                      prefer 2*L*(NOBR-1)*N+2*N+N*NB. */

#line 974 "IB01PD.f"
	    jwork = iu;
#line 975 "IB01PD.f"
	    i__1 = ldun2 - *n;
#line 975 "IB01PD.f"
	    dlaset_("Full", n, &i__1, &c_b66, &c_b66, &r__[nr4pl + (nr2 + *n) 
		    * r_dim1], ldr, (ftnlen)4);
#line 977 "IB01PD.f"
	    i__1 = *ldwork - jwork + 1;
#line 977 "IB01PD.f"
	    dormqr_("Right", "Transpose", n, &ldun2, n, &dwork[ihous], &ldun2,
		     &dwork[itau1], &r__[nr4pl + nr2 * r_dim1], ldr, &dwork[
		    jwork], &i__1, &ierr, (ftnlen)5, (ftnlen)9);
/* Computing MAX */
#line 981 "IB01PD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 981 "IB01PD.f"
	    maxwrk = max(i__1,i__2);
#line 982 "IB01PD.f"
	}
#line 983 "IB01PD.f"
    }

/*     For METH = 'N', find part of the solution (corresponding to A */
/*     and C) and, optionally, for both  METH = 'M',  or  METH = 'N', */
/*     find the residual of the least squares problem that gives the */
/*     covariances,  M*V = N,  where */
/*         (     R_11 ) */
/*     M = (  Y'      ),  N = (  X'   R4'(:,1:L) ),  V = V(n+m*s, n+L), */
/*         (  0   0   ) */
/*     with  M((2*m+L)*s+L, n+m*s),  N((2*m+L)*s+L, n+L),  R4'  being */
/*     stored in the last block-column of  R.  The last  L  rows of  M */
/*     are not explicitly considered. Note that, for efficiency, the */
/*     last  m*s  columns of  M  are in the first positions of arrray  R. */
/*     This permutation does not affect the residual, only the */
/*     solution.  (The solution is not needed for METH = 'M'.) */
/*     Note that R_11 corresponds to the future outputs for both */
/*     METH = 'M', or METH = 'N' approaches.  (For  METH = 'N',  the */
/*     first two block-columns have been interchanged.) */
/*     For  METH = 'N',  A and C are obtained as follows: */
/*     [ A'  C' ] = V(m*s+1:m*s+n,:). */

/*     First, find the  QR  factorization of  Y'(m*s+1:(2*m+L)*s,:) */
/*     and apply the transformations to the corresponding part of N. */
/*     Compress the workspace for N4SID by moving the scalar reflectors */
/*     corresponding to  Q. */
/*     Workspace: need   d*N+2*N; */
/*                prefer d*N+N+N*NB; */
/*     where  d = 0,  for  MOESP,  and  d = 1,  for  N4SID. */

#line 1012 "IB01PD.f"
    if (moesp) {
#line 1013 "IB01PD.f"
	itau = 1;
#line 1014 "IB01PD.f"
    } else {
#line 1015 "IB01PD.f"
	dcopy_(n, &dwork[itau2], &c__1, &dwork[1], &c__1);
#line 1016 "IB01PD.f"
	itau = *n + 1;
#line 1017 "IB01PD.f"
    }

#line 1019 "IB01PD.f"
    jwork = itau + *n;
#line 1020 "IB01PD.f"
    i__1 = *ldwork - jwork + 1;
#line 1020 "IB01PD.f"
    dgeqrf_(&lmnobr, n, &r__[nr2 + nr2 * r_dim1], ldr, &dwork[itau], &dwork[
	    jwork], &i__1, &ierr);

/*     Workspace: need   d*N+N+(N+L); */
/*                prefer d*N+N+(N+L)*NB. */

#line 1026 "IB01PD.f"
    i__1 = *ldwork - jwork + 1;
#line 1026 "IB01PD.f"
    dormqr_("Left", "Transpose", &lmnobr, &npl, n, &r__[nr2 + nr2 * r_dim1], 
	    ldr, &dwork[itau], &r__[nr2 + nr4mn * r_dim1], ldr, &dwork[jwork],
	     &i__1, &ierr, (ftnlen)4, (ftnlen)9);

#line 1030 "IB01PD.f"
    i__1 = *l - 1;
#line 1030 "IB01PD.f"
    i__2 = *l - 1;
#line 1030 "IB01PD.f"
    dlaset_("Lower", &i__1, &i__2, &c_b66, &c_b66, &r__[nr4 + 1 + nr4 * 
	    r_dim1], ldr, (ftnlen)5);

/*     Now, matrix  M  with permuted block-columns has been */
/*     triangularized. */
/*     Compute the reciprocal of the condition number of its */
/*     triangular factor in  R(1:m*s+n,1:m*s+n). */
/*     Workspace: need d*N+3*(M*NOBR+N). */

#line 1038 "IB01PD.f"
    jwork = itau;
#line 1039 "IB01PD.f"
    dtrcon_("1-norm", "Upper", "NonUnit", &mnobrn, &r__[r_offset], ldr, &
	    rcond3, &dwork[jwork], &iwork[1], info, (ftnlen)6, (ftnlen)5, (
	    ftnlen)7);

#line 1042 "IB01PD.f"
    toll = *tol;
#line 1043 "IB01PD.f"
    if (toll <= 0.) {
#line 1043 "IB01PD.f"
	toll = mnobrn * mnobrn * eps;
#line 1043 "IB01PD.f"
    }
#line 1045 "IB01PD.f"
    if (rcond3 > max(toll,thresh)) {

/*        The triangular factor is considered to be of full rank. */
/*        Solve for  V(m*s+1:m*s+n,:),  giving  [ A'  C' ]. */

#line 1050 "IB01PD.f"
	fullr = TRUE_;
#line 1051 "IB01PD.f"
	rankm = mnobrn;
#line 1052 "IB01PD.f"
	if (n4sid) {
#line 1052 "IB01PD.f"
	    dtrsm_("Left", "Upper", "NoTranspose", "NonUnit", n, &npl, &c_b63,
		     &r__[nr2 + nr2 * r_dim1], ldr, &r__[nr2 + nr4mn * r_dim1]
		    , ldr, (ftnlen)4, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 1052 "IB01PD.f"
	}
#line 1055 "IB01PD.f"
    } else {
#line 1056 "IB01PD.f"
	fullr = FALSE_;

/*        Use QR factorization (with pivoting). For METH = 'N', save */
/*        (and then restore) information about the QR factorization of */
/*        Gam,  for later use. Note that  R_11  could be modified by */
/*        MB03OD, but the corresponding part of  N  is also modified */
/*        accordingly. */
/*        Workspace: need   d*N+4*(M*NOBR+N)+1; */
/*                   prefer d*N+3*(M*NOBR+N)+(M*NOBR+N+1)*NB. */

#line 1066 "IB01PD.f"
	i__1 = mnobrn;
#line 1066 "IB01PD.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 1067 "IB01PD.f"
	    iwork[i__] = 0;
#line 1068 "IB01PD.f"
/* L20: */
#line 1068 "IB01PD.f"
	}

#line 1070 "IB01PD.f"
	if (n4sid && (*m > 0 && withb)) {
#line 1070 "IB01PD.f"
	    dlacpy_("Full", &lnobr, n, &r__[nr2 + r_dim1], ldr, &r__[nr4 + 
		    r_dim1], ldr, (ftnlen)4);
#line 1070 "IB01PD.f"
	}
#line 1073 "IB01PD.f"
	jwork = itau + mnobrn;
#line 1074 "IB01PD.f"
	i__1 = mnobrn - 1;
#line 1074 "IB01PD.f"
	dlaset_("Lower", &i__1, &mnobrn, &c_b66, &c_b66, &r__[r_dim1 + 2], 
		ldr, (ftnlen)5);
#line 1076 "IB01PD.f"
	i__1 = *ldwork - jwork + 1;
#line 1076 "IB01PD.f"
	mb03od_("QR", &mnobrn, &mnobrn, &r__[r_offset], ldr, &iwork[1], &toll,
		 &svlmax, &dwork[itau], &rankm, sval, &dwork[jwork], &i__1, &
		ierr, (ftnlen)2);
/* Computing MAX */
#line 1079 "IB01PD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 1079 "IB01PD.f"
	maxwrk = max(i__1,i__2);

/*        Workspace: need   d*N+M*NOBR+N+N+L; */
/*                   prefer d*N+M*NOBR+N+(N+L)*NB. */

#line 1084 "IB01PD.f"
	i__1 = *ldwork - jwork + 1;
#line 1084 "IB01PD.f"
	dormqr_("Left", "Transpose", &mnobrn, &npl, &mnobrn, &r__[r_offset], 
		ldr, &dwork[itau], &r__[nr4mn * r_dim1 + 1], ldr, &dwork[
		jwork], &i__1, &ierr, (ftnlen)4, (ftnlen)9);
/* Computing MAX */
#line 1087 "IB01PD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 1087 "IB01PD.f"
	maxwrk = max(i__1,i__2);
#line 1088 "IB01PD.f"
    }

#line 1090 "IB01PD.f"
    if (withco) {

/*        The residual (transposed) of the least squares solution */
/*        (multiplied by a matrix with orthogonal columns) is stored */
/*        in the rows  RANKM+1:(2*m+L)*s+L  of V,  and it should be */
/*        squared-up for getting the covariance matrices. (Generically, */
/*        RANKM = m*s+n.) */

#line 1098 "IB01PD.f"
	rnrm = 1. / (doublereal) (*nsmpl);
#line 1099 "IB01PD.f"
	if (moesp) {
#line 1100 "IB01PD.f"
	    i__1 = lmmnol - rankm;
#line 1100 "IB01PD.f"
	    dsyrk_("Upper", "Transpose", &npl, &i__1, &rnrm, &r__[rankm + 1 + 
		    nr4mn * r_dim1], ldr, &c_b66, &r__[r_offset], ldr, (
		    ftnlen)5, (ftnlen)9);
#line 1102 "IB01PD.f"
	    dlacpy_("Upper", n, n, &r__[r_offset], ldr, &q[q_offset], ldq, (
		    ftnlen)5);
#line 1103 "IB01PD.f"
	    dlacpy_("Full", n, l, &r__[(*n + 1) * r_dim1 + 1], ldr, &s[
		    s_offset], lds, (ftnlen)4);
#line 1104 "IB01PD.f"
	    dlacpy_("Upper", l, l, &r__[*n + 1 + (*n + 1) * r_dim1], ldr, &ry[
		    ry_offset], ldry, (ftnlen)5);
#line 1105 "IB01PD.f"
	} else {
#line 1106 "IB01PD.f"
	    i__1 = lmmnol - rankm;
#line 1106 "IB01PD.f"
	    dsyrk_("Upper", "Transpose", &npl, &i__1, &rnrm, &r__[rankm + 1 + 
		    nr4mn * r_dim1], ldr, &c_b66, &dwork[jwork], &npl, (
		    ftnlen)5, (ftnlen)9);
#line 1108 "IB01PD.f"
	    dlacpy_("Upper", n, n, &dwork[jwork], &npl, &q[q_offset], ldq, (
		    ftnlen)5);
#line 1109 "IB01PD.f"
	    dlacpy_("Full", n, l, &dwork[jwork + *n * npl], &npl, &s[s_offset]
		    , lds, (ftnlen)4);
#line 1111 "IB01PD.f"
	    dlacpy_("Upper", l, l, &dwork[jwork + *n * (npl + 1)], &npl, &ry[
		    ry_offset], ldry, (ftnlen)5);
#line 1113 "IB01PD.f"
	}
#line 1114 "IB01PD.f"
	ma02ed_("Upper", n, &q[q_offset], ldq, (ftnlen)5);
#line 1115 "IB01PD.f"
	ma02ed_("Upper", l, &ry[ry_offset], ldry, (ftnlen)5);

/*        Check the magnitude of the residual. */

#line 1119 "IB01PD.f"
	i__1 = lmmnol - rankm;
#line 1119 "IB01PD.f"
	rnrm = dlange_("1-norm", &i__1, &npl, &r__[rankm + 1 + nr4mn * r_dim1]
		, ldr, &dwork[jwork], (ftnlen)6);
#line 1121 "IB01PD.f"
	if (rnrm < thresh) {
#line 1121 "IB01PD.f"
	    *iwarn = 5;
#line 1121 "IB01PD.f"
	}
#line 1123 "IB01PD.f"
    }

#line 1125 "IB01PD.f"
    if (n4sid) {
#line 1126 "IB01PD.f"
	if (! fullr) {
#line 1127 "IB01PD.f"
	    *iwarn = 4;

/*           Compute part of the solution of the least squares problem, */
/*           M*V = N,  for the rank-deficient problem. */
/*           Remark: this computation should not be performed before the */
/*           symmetric updating operation above. */
/*           Workspace: need   M*NOBR+3*N+L; */
/*                      prefer larger. */

#line 1136 "IB01PD.f"
	    i__1 = *ldwork - jwork + 1;
#line 1136 "IB01PD.f"
	    mb03od_("No QR", n, n, &r__[nr2 + nr2 * r_dim1], ldr, &iwork[1], &
		    toll1, &svlmax, &dwork[itau], &rankm, sval, &dwork[jwork],
		     &i__1, &ierr, (ftnlen)5);
#line 1139 "IB01PD.f"
	    i__1 = *ldwork - jwork + 1;
#line 1139 "IB01PD.f"
	    mb02qy_(n, n, &npl, &rankm, &r__[nr2 + nr2 * r_dim1], ldr, &iwork[
		    1], &r__[nr2 + nr4mn * r_dim1], ldr, &dwork[itau + mnobr],
		     &dwork[jwork], &i__1, info);
/* Computing MAX */
#line 1142 "IB01PD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 1142 "IB01PD.f"
	    maxwrk = max(i__1,i__2);
#line 1143 "IB01PD.f"
	    jwork = itau;
#line 1144 "IB01PD.f"
	    if (*m > 0 && withb) {
#line 1144 "IB01PD.f"
		dlacpy_("Full", &lnobr, n, &r__[nr4 + r_dim1], ldr, &r__[nr2 
			+ r_dim1], ldr, (ftnlen)4);
#line 1144 "IB01PD.f"
	    }
#line 1147 "IB01PD.f"
	}

#line 1149 "IB01PD.f"
	if (withc) {

/*           Obtain  A  and  C,  noting that block-permutations have been */
/*           implicitly used. */

#line 1154 "IB01PD.f"
	    ma02ad_("Full", n, n, &r__[nr2 + nr4mn * r_dim1], ldr, &a[
		    a_offset], lda, (ftnlen)4);
#line 1155 "IB01PD.f"
	    ma02ad_("Full", n, l, &r__[nr2 + (nr4mn + *n) * r_dim1], ldr, &
		    c__[c_offset], ldc, (ftnlen)4);
#line 1156 "IB01PD.f"
	} else {

/*           Use the given  A  and  C. */

#line 1160 "IB01PD.f"
	    ma02ad_("Full", n, n, &a[a_offset], lda, &r__[nr2 + nr4mn * 
		    r_dim1], ldr, (ftnlen)4);
#line 1161 "IB01PD.f"
	    ma02ad_("Full", l, n, &c__[c_offset], ldc, &r__[nr2 + (nr4mn + *n)
		     * r_dim1], ldr, (ftnlen)4);
#line 1162 "IB01PD.f"
	}

#line 1164 "IB01PD.f"
	if (*m > 0 && withb) {

/*           Obtain  B  and  D. */
/*           First, compute the transpose of the matrix K as */
/*           N(1:m*s,:) - M(1:m*s,m*s+1:m*s+n)*[A'  C'],  in the first */
/*           m*s  rows of  R(1,NR4MN). */

#line 1171 "IB01PD.f"
	    dgemm_("NoTranspose", "NoTranspose", &mnobr, &npl, n, &c_b173, &
		    r__[nr2 * r_dim1 + 1], ldr, &r__[nr2 + nr4mn * r_dim1], 
		    ldr, &c_b63, &r__[nr4mn * r_dim1 + 1], ldr, (ftnlen)11, (
		    ftnlen)11);

/*           Denote   M = pinv(GaL)  and construct */

/*                    [ [ A ]   -1   ]                      [ R ] */
/*           and  L = [ [   ]  R   0 ] Q',  where Gam = Q * [   ]. */
/*                    [ [ C ]        ]                      [ 0 ] */

/*           Then, solve the least squares problem. */

#line 1183 "IB01PD.f"
	    dlacpy_("Full", n, n, &a[a_offset], lda, &r__[nr2 + nr4 * r_dim1],
		     ldr, (ftnlen)4);
#line 1184 "IB01PD.f"
	    dlacpy_("Full", l, n, &c__[c_offset], ldc, &r__[nr2 + *n + nr4 * 
		    r_dim1], ldr, (ftnlen)4);
#line 1185 "IB01PD.f"
	    dtrsm_("Right", "Upper", "NoTranspose", "NonUnit", &npl, n, &
		    c_b63, &r__[nr2 + r_dim1], ldr, &r__[nr2 + nr4 * r_dim1], 
		    ldr, (ftnlen)5, (ftnlen)5, (ftnlen)11, (ftnlen)7);
#line 1187 "IB01PD.f"
	    dlaset_("Full", &npl, &lnobrn, &c_b66, &c_b66, &r__[nr2 + (nr4 + *
		    n) * r_dim1], ldr, (ftnlen)4);

/*           Workspace: need 2*N+L; prefer N + (N+L)*NB. */

#line 1192 "IB01PD.f"
	    i__1 = *ldwork - jwork + 1;
#line 1192 "IB01PD.f"
	    dormqr_("Right", "Transpose", &npl, &lnobr, n, &r__[nr2 + r_dim1],
		     ldr, &dwork[1], &r__[nr2 + nr4 * r_dim1], ldr, &dwork[
		    jwork], &i__1, &ierr, (ftnlen)5, (ftnlen)9);

/*           Obtain the matrix  K  by transposition, and find  B  and  D. */
/*           Workspace: need   NOBR*(M*(N+L))**2+M*NOBR*(N+L)+ */
/*                             max((N+L)**2,4*M*(N+L)+1); */
/*                      prefer larger. */

#line 1201 "IB01PD.f"
	    ma02ad_("Full", &mnobr, &npl, &r__[nr4mn * r_dim1 + 1], ldr, &r__[
		    nr2 + nr3 * r_dim1], ldr, (ftnlen)4);
/* Computing 2nd power */
#line 1203 "IB01PD.f"
	    i__1 = npl;
#line 1203 "IB01PD.f"
	    ix = mnobr * (i__1 * i__1) * *m + 1;
#line 1204 "IB01PD.f"
	    jwork = ix + mnobr * npl;
#line 1205 "IB01PD.f"
	    i__1 = mnobr * npl;
#line 1205 "IB01PD.f"
	    i__2 = *ldwork - jwork + 1;
#line 1205 "IB01PD.f"
	    ib01px_(jobpy, nobr, n, m, l, &r__[r_offset], ldr, &o[o_offset], 
		    ldo, &r__[nr2 + nr4 * r_dim1], ldr, &r__[nr4pl + nr2 * 
		    r_dim1], ldr, &r__[nr2 + nr3 * r_dim1], ldr, &dwork[1], &
		    i__1, &dwork[ix], &b[b_offset], ldb, &d__[d_offset], ldd, 
		    tol, &iwork[1], &dwork[jwork], &i__2, &iwarnl, info, (
		    ftnlen)1);
#line 1210 "IB01PD.f"
	    if (*info != 0) {
#line 1210 "IB01PD.f"
		return 0;
#line 1210 "IB01PD.f"
	    }
#line 1212 "IB01PD.f"
	    *iwarn = max(*iwarn,iwarnl);
/* Computing MAX */
#line 1213 "IB01PD.f"
	    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 1213 "IB01PD.f"
	    maxwrk = max(i__1,i__2);
#line 1214 "IB01PD.f"
	    rcond4 = dwork[jwork + 1];

#line 1216 "IB01PD.f"
	}
#line 1217 "IB01PD.f"
    }

#line 1219 "IB01PD.f"
L30:

/*     Return optimal workspace in  DWORK(1)  and reciprocal condition */
/*     numbers in the next locations. */

#line 1224 "IB01PD.f"
    dwork[1] = (doublereal) maxwrk;
#line 1225 "IB01PD.f"
    dwork[2] = rcond1;
#line 1226 "IB01PD.f"
    dwork[3] = rcond2;
#line 1227 "IB01PD.f"
    dwork[4] = rcond3;
#line 1228 "IB01PD.f"
    dwork[5] = rcond4;
#line 1229 "IB01PD.f"
    return 0;

/* *** Last line of IB01PD *** */
} /* ib01pd_ */

