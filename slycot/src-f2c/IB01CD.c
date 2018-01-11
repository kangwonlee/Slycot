#line 1 "IB01CD.f"
/* IB01CD.f -- translated by f2c (version 20100827).
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

#line 1 "IB01CD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__0 = 0;
static doublereal c_b57 = 1.;
static doublereal c_b58 = 0.;

/* Subroutine */ int ib01cd_(char *jobx0, char *comuse, char *job, integer *n,
	 integer *m, integer *l, integer *nsmp, doublereal *a, integer *lda, 
	doublereal *b, integer *ldb, doublereal *c__, integer *ldc, 
	doublereal *d__, integer *ldd, doublereal *u, integer *ldu, 
	doublereal *y, integer *ldy, doublereal *x0, doublereal *v, integer *
	ldv, doublereal *tol, integer *iwork, doublereal *dwork, integer *
	ldwork, integer *iwarn, integer *info, ftnlen jobx0_len, ftnlen 
	comuse_len, ftnlen job_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, d_dim1, 
	    d_offset, u_dim1, u_offset, v_dim1, v_offset, y_dim1, y_offset, 
	    i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer i__, ia, ib, ic, lm, iq, ln, nm, nn, n2m, ldw;
    static doublereal dum[1];
    static integer iwi, iwr, ncp1, ldw2, ldw3;
    static char jobd[1];
    static integer ncol, ierr, itau, mtmp;
    extern /* Subroutine */ int ib01qd_(char *, char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen, ftnlen), ib01rd_(char *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *, 
	    integer *, ftnlen), dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical usebd;
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), tb01wd_(integer *,
	     integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);
    static doublereal rcond;
    static logical withb;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static logical withd;
    static integer isize, nsmpl, jwork;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static logical withx0, maxdia, compbd;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static logical maxdim;
    static doublereal rcondu;
    static integer iwarnl, minsmp, minwrk, maxwrk, minwls;


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

/*     To estimate the initial state and, optionally, the system matrices */
/*     B  and  D  of a linear time-invariant (LTI) discrete-time system, */
/*     given the system matrices  (A,B,C,D),  or (when  B  and  D  are */
/*     estimated) only the matrix pair  (A,C),  and the input and output */
/*     trajectories of the system. The model structure is : */

/*           x(k+1) = Ax(k) + Bu(k),   k >= 0, */
/*           y(k)   = Cx(k) + Du(k), */

/*     where  x(k)  is the  n-dimensional state vector (at time k), */
/*            u(k)  is the  m-dimensional input vector, */
/*            y(k)  is the  l-dimensional output vector, */
/*     and  A, B, C, and D  are real matrices of appropriate dimensions. */
/*     The input-output data can internally be processed sequentially. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOBX0   CHARACTER*1 */
/*             Specifies whether or not the initial state should be */
/*             computed, as follows: */
/*             = 'X':  compute the initial state x(0); */
/*             = 'N':  do not compute the initial state (possibly, */
/*                     because x(0) is known to be zero). */

/*     COMUSE  CHARACTER*1 */
/*             Specifies whether the system matrices B and D should be */
/*             computed or used, as follows: */
/*             = 'C':  compute the system matrices B and D, as specified */
/*                     by JOB; */
/*             = 'U':  use the system matrices B and D, as specified by */
/*                     JOB; */
/*             = 'N':  do not compute/use the matrices B and D. */
/*             If  JOBX0 = 'N'  and  COMUSE <> 'N',  then  x(0)  is set */
/*             to zero. */
/*             If  JOBX0 = 'N'  and  COMUSE =  'N',  then  x(0)  is */
/*             neither computed nor set to zero. */

/*     JOB     CHARACTER*1 */
/*             If  COMUSE = 'C'  or  'U',  specifies which of the system */
/*             matrices  B and D  should be computed or used, as follows: */
/*             = 'B':  compute/use the matrix B only (D is known to be */
/*                     zero); */
/*             = 'D':  compute/use the matrices B and D. */
/*             The value of  JOB  is irrelevant if  COMUSE = 'N'  or if */
/*             JOBX0 = 'N'  and  COMUSE = 'U'. */
/*             The combinations of options, the data used, and the */
/*             returned results, are given in the table below, where */
/*             '*'  denotes an irrelevant value. */

/*              JOBX0   COMUSE    JOB     Data used    Returned results */
/*             ---------------------------------------------------------- */
/*                X       C        B       A,C,u,y          x,B */
/*                X       C        D       A,C,u,y          x,B,D */
/*                N       C        B       A,C,u,y          x=0,B */
/*                N       C        D       A,C,u,y          x=0,B,D */
/*             ---------------------------------------------------------- */
/*                X       U        B      A,B,C,u,y            x */
/*                X       U        D      A,B,C,D,u,y          x */
/*                N       U        *          -               x=0 */
/*             ---------------------------------------------------------- */
/*                X       N        *        A,C,y              x */
/*                N       N        *          -                - */
/*             ---------------------------------------------------------- */

/*             For  JOBX0 = 'N'  and  COMUSE = 'N',  the routine just */
/*             sets  DWORK(1)  to 2 and  DWORK(2)  to 1, and returns */
/*             (see the parameter DWORK). */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the system.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     L       (input) INTEGER */
/*             The number of system outputs.  L > 0. */

/*     NSMP    (input) INTEGER */
/*             The number of rows of matrices  U  and  Y  (number of */
/*             samples,  t). */
/*             NSMP >= 0,            if  JOBX0 = 'N'  and  COMUSE <> 'C'; */
/*             NSMP >= N,            if  JOBX0 = 'X'  and  COMUSE <> 'C'; */
/*             NSMP >= N*M + a + e,  if  COMUSE = 'C', */
/*             where   a = 0,  if  JOBX0 = 'N'; */
/*                     a = N,  if  JOBX0 = 'X'; */
/*                     e = 0,  if  JOBX0 = 'X'  and  JOB = 'B'; */
/*                     e = 1,  if  JOBX0 = 'N'  and  JOB = 'B'; */
/*                     e = M,  if  JOB   = 'D'. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             If  JOBX0 = 'X'  or  COMUSE = 'C',  the leading N-by-N */
/*             part of this array must contain the system state matrix A. */
/*             If  N = 0,  or  JOBX0 = 'N'  and  COMUSE <> 'C',  this */
/*             array is not referenced. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. */
/*             LDA >= MAX(1,N),  if  JOBX0 = 'X'  or   COMUSE =  'C'; */
/*             LDA >= 1,         if  JOBX0 = 'N'  and  COMUSE <> 'C'. */

/*     B       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDB,M) */
/*             If  JOBX0 = 'X'  and  COMUSE = 'U',  B  is an input */
/*             parameter and, on entry, the leading N-by-M part of this */
/*             array must contain the system input matrix  B. */
/*             If  COMUSE = 'C',  B  is an output parameter and, on exit, */
/*             if  INFO = 0,  the leading N-by-M part of this array */
/*             contains the estimated system input matrix  B. */
/*             If  min(N,M) = 0,  or  JOBX0 = 'N'  and  COMUSE = 'U', */
/*             or  COMUSE = 'N',  this array is not referenced. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. */
/*             LDB >= MAX(1,N),  if  M > 0,  COMUSE = 'U',  JOBX0 = 'X', */
/*                               or  M > 0,  COMUSE = 'C'; */
/*             LDB >= 1,         if  min(N,M) = 0,  or  COMUSE = 'N', */
/*                               or  JOBX0  = 'N'  and  COMUSE = 'U'. */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             If  JOBX0 = 'X'  or  COMUSE = 'C',  the leading L-by-N */
/*             part of this array must contain the system output */
/*             matrix  C. */
/*             If  N = 0,  or  JOBX0 = 'N'  and  COMUSE <> 'C',  this */
/*             array is not referenced. */

/*     LDC     INTEGER */
/*             The leading dimension of the array C. */
/*             LDC >= L,  if  N > 0, and  JOBX0 = 'X'  or  COMUSE = 'C'; */
/*             LDC >= 1,  if  N = 0, or  JOBX0 = 'N'  and  COMUSE <> 'C'. */

/*     D       (input or output) DOUBLE PRECISION array, dimension */
/*             (LDD,M) */
/*             If  JOBX0 = 'X',  COMUSE = 'U',  and  JOB = 'D',  D  is an */
/*             input parameter and, on entry, the leading L-by-M part of */
/*             this array must contain the system input-output matrix  D. */
/*             If  COMUSE = 'C'  and  JOB = 'D',  D  is an output */
/*             parameter and, on exit, if  INFO = 0,  the leading */
/*             L-by-M part of this array contains the estimated system */
/*             input-output matrix  D. */
/*             If  M = 0,  or  JOBX0 = 'N'  and  COMUSE = 'U',  or */
/*             COMUSE = 'N',  or  JOB = 'B',  this array is not */
/*             referenced. */

/*     LDD     INTEGER */
/*             The leading dimension of the array D. */
/*             LDD >= L,  if  M > 0,   JOBX0 = 'X',  COMUSE = 'U',  and */
/*                                                   JOB = 'D',  or */
/*                        if  M > 0,  COMUSE = 'C',  and  JOB = 'D'; */
/*             LDD >= 1,  if  M = 0,  or  JOBX0 = 'N'  and  COMUSE = 'U', */
/*                        or  COMUSE = 'N',  or  JOB = 'B'. */

/*     U       (input or input/output) DOUBLE PRECISION array, dimension */
/*             (LDU,M) */
/*             On entry, if  COMUSE = 'C',  or  JOBX0 = 'X'  and */
/*             COMUSE = 'U',  the leading NSMP-by-M part of this array */
/*             must contain the t-by-m input-data sequence matrix  U, */
/*             U = [u_1 u_2 ... u_m].  Column  j  of  U  contains the */
/*             NSMP  values of the j-th input component for consecutive */
/*             time increments. */
/*             On exit, if  COMUSE = 'C'  and  JOB = 'D',  the leading */
/*             NSMP-by-M part of this array contains details of the */
/*             QR factorization of the t-by-m matrix  U,  possibly */
/*             computed sequentially (see METHOD). */
/*             If  COMUSE = 'C'  and  JOB = 'B',  or  COMUSE = 'U',  this */
/*             array is unchanged on exit. */
/*             If  M = 0,  or  JOBX0 = 'N'  and  COMUSE = 'U',  or */
/*             COMUSE = 'N',  this array is not referenced. */

/*     LDU     INTEGER */
/*             The leading dimension of the array U. */
/*             LDU >= MAX(1,NSMP),  if  M > 0    and  COMUSE = 'C'  or */
/*                                  JOBX0 = 'X'  and  COMUSE = 'U; */
/*             LDU >= 1,            if  M = 0,   or   COMUSE = 'N',  or */
/*                                  JOBX0 = 'N'  and  COMUSE = 'U'. */

/*     Y       (input) DOUBLE PRECISION array, dimension (LDY,L) */
/*             On entry, if  JOBX0 = 'X'  or  COMUSE = 'C',  the leading */
/*             NSMP-by-L part of this array must contain the t-by-l */
/*             output-data sequence matrix  Y,  Y = [y_1 y_2 ... y_l]. */
/*             Column  j  of  Y  contains the  NSMP  values of the j-th */
/*             output component for consecutive time increments. */
/*             If  JOBX0 = 'N'  and  COMUSE <> 'C',  this array is not */
/*             referenced. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y. */
/*             LDY >= MAX(1,NSMP),  if  JOBX0 = 'X'  or   COMUSE = 'C; */
/*             LDY >= 1,            if  JOBX0 = 'N'  and  COMUSE <> 'C'. */

/*     X0      (output) DOUBLE PRECISION array, dimension (N) */
/*             If  INFO = 0  and  JOBX0 = 'X',  this array contains the */
/*             estimated initial state of the system,  x(0). */
/*             If  JOBX0 = 'N'  and  COMUSE = 'C',  this array is used as */
/*             workspace and finally it is set to zero. */
/*             If  JOBX0 = 'N'  and  COMUSE = 'U',  then  x(0)  is set to */
/*             zero without any calculations. */
/*             If  JOBX0 = 'N'  and  COMUSE = 'N',  this array is not */
/*             referenced. */

/*     V       (output) DOUBLE PRECISION array, dimension (LDV,N) */
/*             On exit, if  INFO = 0  or 2,  JOBX0 = 'X'  or */
/*             COMUSE = 'C',  the leading N-by-N part of this array */
/*             contains the orthogonal matrix V of a real Schur */
/*             factorization of the matrix  A. */
/*             If  JOBX0 = 'N'  and  COMUSE <> 'C',  this array is not */
/*             referenced. */

/*     LDV     INTEGER */
/*             The leading dimension of the array V. */
/*             LDV >= MAX(1,N),  if  JOBX0 = 'X'  or   COMUSE =  'C; */
/*             LDV >= 1,         if  JOBX0 = 'N'  and  COMUSE <> 'C'. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             The tolerance to be used for estimating the rank of */
/*             matrices. If the user sets  TOL > 0,  then the given value */
/*             of  TOL  is used as a lower bound for the reciprocal */
/*             condition number;  a matrix whose estimated condition */
/*             number is less than  1/TOL  is considered to be of full */
/*             rank.  If the user sets  TOL <= 0,  then  EPS  is used */
/*             instead, where  EPS  is the relative machine precision */
/*             (see LAPACK Library routine DLAMCH).  TOL <= 1. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (LIWORK), where */
/*             LIWORK >= 0,          if  JOBX0 = 'N'  and  COMUSE <> 'C'; */
/*             LIWORK >= N,          if  JOBX0 = 'X'  and  COMUSE <> 'C'; */
/*             LIWORK >= N*M + a,        if COMUSE = 'C' and JOB = 'B', */
/*             LIWORK >= max(N*M + a,M), if COMUSE = 'C' and JOB = 'D', */
/*             with  a = 0,  if  JOBX0 = 'N'; */
/*                   a = N,  if  JOBX0 = 'X'. */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if  INFO = 0,  DWORK(1) returns the optimal value */
/*             of LDWORK;  DWORK(2)  contains the reciprocal condition */
/*             number of the triangular factor of the QR factorization of */
/*             the matrix  W2,  if  COMUSE = 'C',  or of the matrix */
/*             Gamma,  if  COMUSE = 'U'  (see METHOD); if  JOBX0 = 'N' */
/*             and  COMUSE <> 'C',  DWORK(2)  is set to one; */
/*             if  COMUSE = 'C',  M > 0,  and  JOB = 'D',   DWORK(3) */
/*             contains the reciprocal condition number of the triangular */
/*             factor of the QR factorization of  U;  denoting */
/*                g = 2,  if  JOBX0  = 'X'  and  COMUSE <> 'C'  or */
/*                            COMUSE = 'C'  and  M = 0  or   JOB = 'B', */
/*                g = 3,  if  COMUSE = 'C'  and  M > 0  and  JOB = 'D', */
/*             then  DWORK(i), i = g+1:g+N*N, */
/*                   DWORK(j), j = g+1+N*N:g+N*N+L*N,  and */
/*                   DWORK(k), k = g+1+N*N+L*N:g+N*N+L*N+N*M, */
/*             contain the transformed system matrices  At, Ct, and Bt, */
/*             respectively, corresponding to the real Schur form of the */
/*             given system state matrix  A,  i.e., */
/*                At = V'*A*V,  Bt = V'*B,  Ct = C*V. */
/*             The matrices  At, Ct, Bt  are not computed if  JOBX0 = 'N' */
/*             and  COMUSE <> 'C'. */
/*             On exit, if  INFO = -26,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 2,  if  JOBX0 = 'N'  and  COMUSE <> 'C',  or */
/*                           if  max( N, M ) = 0. */
/*             Otherwise, */
/*             LDWORK >= LDW1 + N*( N + M + L ) + */
/*                              max( 5*N, LDW1, min( LDW2, LDW3 ) ), */
/*             where, if  COMUSE = 'C',  then */
/*             LDW1 = 2,          if  M = 0  or   JOB = 'B', */
/*             LDW1 = 3,          if  M > 0  and  JOB = 'D', */
/*             LDWa = t*L*(r + 1) + max( N + max( d, f ), 6*r ), */
/*             LDW2 = LDWa,       if  M = 0  or  JOB = 'B', */
/*             LDW2 = max( LDWa, t*L*(r + 1) + 2*M*M + 6*M ), */
/*                                if  M > 0  and JOB = 'D', */
/*             LDWb = (b + r)*(r + 1) + */
/*                     max( q*(r + 1) + N*N*M + c + max( d, f ), 6*r ), */
/*             LDW3 = LDWb,       if  M = 0  or  JOB = 'B', */
/*             LDW3 = max( LDWb, (b + r)*(r + 1) + 2*M*M + 6*M ), */
/*                                if  M > 0  and JOB = 'D', */
/*                r = N*M + a, */
/*                a = 0,                  if  JOBX0 = 'N', */
/*                a = N,                  if  JOBX0 = 'X'; */
/*                b = 0,                  if  JOB   = 'B', */
/*                b = L*M,                if  JOB   = 'D'; */
/*                c = 0,                  if  JOBX0 = 'N', */
/*                c = L*N,                if  JOBX0 = 'X'; */
/*                d = 0,                  if  JOBX0 = 'N', */
/*                d = 2*N*N + N,          if  JOBX0 = 'X'; */
/*                f = 2*r,                if  JOB   = 'B'   or  M = 0, */
/*                f = M + max( 2*r, M ),  if  JOB   = 'D'  and  M > 0; */
/*                q = b + r*L; */
/*             and, if  JOBX0 = 'X'  and  COMUSE <> 'C',  then */
/*             LDW1 = 2, */
/*             LDW2 = t*L*(N + 1) + 2*N + max( 2*N*N, 4*N ), */
/*             LDW3 = N*(N + 1) + 2*N + max( q*(N + 1) + 2*N*N + L*N, */
/*                                           4*N ), */
/*                q = N*L. */
/*             For good performance,  LDWORK  should be larger. */
/*             If  LDWORK >= LDW2,  or if  COMUSE = 'C'  and */
/*                 LDWORK >= t*L*(r + 1) + (b + r)*(r + 1) + N*N*M + c + */
/*                           max( d, f ), */
/*             then standard QR factorizations of the matrices  U  and/or */
/*             W2,  if  COMUSE = 'C',  or of the matrix  Gamma,  if */
/*             JOBX0 = 'X'  and  COMUSE <> 'C'  (see METHOD), are used. */
/*             Otherwise, the QR factorizations are computed sequentially */
/*             by performing  NCYCLE  cycles, each cycle (except possibly */
/*             the last one) processing  s < t  samples, where  s  is */
/*             chosen by equating  LDWORK  to the first term of  LDWb, */
/*             if  COMUSE = 'C',  or of  LDW3,  if  COMUSE <> 'C',  for */
/*             q  replaced by  s*L.  (s  is larger than or equal to the */
/*             minimum value of  NSMP.)  The computational effort may */
/*             increase and the accuracy may slightly decrease with the */
/*             decrease of  s.  Recommended value is  LDWORK = LDW2, */
/*             assuming a large enough cache size, to also accommodate */
/*             A,  (B,)  C,  (D,)  U,  and  Y. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             = 0:  no warning; */
/*             = 4:  the least squares problem to be solved has a */
/*                   rank-deficient coefficient matrix; */
/*             = 6:  the matrix  A  is unstable;  the estimated  x(0) */
/*                   and/or  B and D  could be inaccurate. */
/*             NOTE: the value 4 of  IWARN  has no significance for the */
/*                   identification problem. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  if the QR algorithm failed to compute all the */
/*                   eigenvalues of the matrix A (see LAPACK Library */
/*                   routine DGEES); the locations  DWORK(i),  for */
/*                   i = g+1:g+N*N,  contain the partially converged */
/*                   Schur form; */
/*             = 2:  the singular value decomposition (SVD) algorithm did */
/*                   not converge. */

/*     METHOD */

/*     Matrix  A  is initially reduced to a real Schur form, A = V*At*V', */
/*     and the given system matrices are transformed accordingly. For the */
/*     reduced system, an extension and refinement of the method in [1,2] */
/*     is used. Specifically, for  JOBX0 = 'X',  COMUSE = 'C',  and */
/*     JOB = 'D',  denoting */

/*           X = [ vec(D')' vec(B)' x0' ]', */

/*     where  vec(M)  is the vector obtained by stacking the columns of */
/*     the matrix  M,  then  X  is the least squares solution of the */
/*     system  S*X = vec(Y),  with the matrix  S = [ diag(U)  W ], */
/*     defined by */

/*           ( U         |     | ... |     |     | ... |     |         ) */
/*           (   U       |  11 | ... |  n1 |  12 | ... |  nm |         ) */
/*       S = (     :     | y   | ... | y   | y   | ... | y   | P*Gamma ), */
/*           (       :   |     | ... |     |     | ... |     |         ) */
/*           (         U |     | ... |     |     | ... |     |         ) */
/*                                                                     ij */
/*     diag(U)  having  L  block rows and columns.  In this formula,  y */
/*     are the outputs of the system for zero initial state computed */
/*     using the following model, for j = 1:m, and for i = 1:n, */
/*            ij          ij                    ij */
/*           x  (k+1) = Ax  (k) + e_i u_j(k),  x  (0) = 0, */

/*            ij          ij */
/*           y  (k)   = Cx  (k), */

/*     where  e_i  is the i-th n-dimensional unit vector,  Gamma  is */
/*     given by */

/*                (     C     ) */
/*                (    C*A    ) */
/*        Gamma = (   C*A^2   ), */
/*                (     :     ) */
/*                ( C*A^(t-1) ) */

/*     and  P  is a permutation matrix that groups together the rows of */
/*     Gamma  depending on the same row of  C,  namely */
/*     [ c_j;  c_j*A;  c_j*A^2; ...  c_j*A^(t-1) ],  for j = 1:L. */
/*     The first block column,  diag(U),  is not explicitly constructed, */
/*     but its structure is exploited. The last block column is evaluated */
/*     using powers of A with exponents 2^k. No interchanges are applied. */
/*     A special QR decomposition of the matrix  S  is computed. Let */
/*     U = q*[ r' 0 ]'  be the QR decomposition of  U,  if  M > 0,  where */
/*     r  is  M-by-M.   Then,  diag(q')  is applied to  W  and  vec(Y). */
/*     The block-rows of  S  and  vec(Y)  are implicitly permuted so that */
/*     matrix  S  becomes */

/*        ( diag(r)  W1 ) */
/*        (    0     W2 ), */

/*     where  W1  has L*M rows. Then, the QR decomposition of  W2 is */
/*     computed (sequentially, if  M > 0) and used to obtain  B  and  x0. */
/*     The intermediate results and the QR decomposition of  U  are */
/*     needed to find  D.  If a triangular factor is too ill conditioned, */
/*     then singular value decomposition (SVD) is employed. SVD is not */
/*     generally needed if the input sequence is sufficiently */
/*     persistently exciting and  NSMP  is large enough. */
/*     If the matrix  W  cannot be stored in the workspace (i.e., */
/*     LDWORK < LDW2),  the QR decompositions of  W2  and  U  are */
/*     computed sequentially. */
/*     For  JOBX0 = 'N'  and  COMUSE = 'C',  or  JOB = 'B',  a simpler */
/*     problem is solved efficiently. */

/*     For  JOBX0 = 'X'  and  COMUSE <> 'C',  a simpler method is used. */
/*     Specifically, the output y0(k) of the system for zero initial */
/*     state is computed for k = 0, 1, ...,  t-1 using the given model. */
/*     Then the following least squares problem is solved for x(0) */

/*                         (   y(0) - y0(0)   ) */
/*                         (   y(1) - y0(1)   ) */
/*        Gamma * x(0)  =  (        :         ). */
/*                         (        :         ) */
/*                         ( y(t-1) - y0(t-1) ) */

/*     The coefficient matrix  Gamma  is evaluated using powers of A with */
/*     exponents 2^k. The QR decomposition of this matrix is computed. */
/*     If its triangular factor  R  is too ill conditioned, then singular */
/*     value decomposition of  R  is used. */
/*     If the coefficient matrix cannot be stored in the workspace (i.e., */
/*     LDWORK < LDW2),  the QR decomposition is computed sequentially. */


/*     REFERENCES */

/*     [1] Verhaegen M., and Varga, A. */
/*         Some Experience with the MOESP Class of Subspace Model */
/*         Identification Methods in Identifying the BO105 Helicopter. */
/*         Report TR R165-94, DLR Oberpfaffenhofen, 1994. */

/*     [2] Sima, V., and Varga, A. */
/*         RASP-IDENT : Subspace Model Identification Programs. */
/*         Deutsche Forschungsanstalt fur Luft- und Raumfahrt e. V., */
/*         Report TR R888-94, DLR Oberpfaffenhofen, Oct. 1994. */

/*     NUMERICAL ASPECTS */

/*     The implemented method is numerically stable. */

/*     FURTHER COMMENTS */

/*     The algorithm for computing the system matrices  B  and  D  is */
/*     less efficient than the MOESP or N4SID algorithms implemented in */
/*     SLICOT Library routines IB01BD/IB01PD, because a large least */
/*     squares problem has to be solved, but the accuracy is better, as */
/*     the computed matrices  B  and  D  are fitted to the input and */
/*     output trajectories. However, if matrix  A  is unstable, the */
/*     computed matrices  B  and  D  could be inaccurate. */

/*     CONTRIBUTOR */

/*     V. Sima, Research Institute for Informatics, Bucharest, Apr. 2000. */

/*     REVISIONS */

/*     - */

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

/*     Check the input parameters. */

#line 535 "IB01CD.f"
    /* Parameter adjustments */
#line 535 "IB01CD.f"
    a_dim1 = *lda;
#line 535 "IB01CD.f"
    a_offset = 1 + a_dim1;
#line 535 "IB01CD.f"
    a -= a_offset;
#line 535 "IB01CD.f"
    b_dim1 = *ldb;
#line 535 "IB01CD.f"
    b_offset = 1 + b_dim1;
#line 535 "IB01CD.f"
    b -= b_offset;
#line 535 "IB01CD.f"
    c_dim1 = *ldc;
#line 535 "IB01CD.f"
    c_offset = 1 + c_dim1;
#line 535 "IB01CD.f"
    c__ -= c_offset;
#line 535 "IB01CD.f"
    d_dim1 = *ldd;
#line 535 "IB01CD.f"
    d_offset = 1 + d_dim1;
#line 535 "IB01CD.f"
    d__ -= d_offset;
#line 535 "IB01CD.f"
    u_dim1 = *ldu;
#line 535 "IB01CD.f"
    u_offset = 1 + u_dim1;
#line 535 "IB01CD.f"
    u -= u_offset;
#line 535 "IB01CD.f"
    y_dim1 = *ldy;
#line 535 "IB01CD.f"
    y_offset = 1 + y_dim1;
#line 535 "IB01CD.f"
    y -= y_offset;
#line 535 "IB01CD.f"
    --x0;
#line 535 "IB01CD.f"
    v_dim1 = *ldv;
#line 535 "IB01CD.f"
    v_offset = 1 + v_dim1;
#line 535 "IB01CD.f"
    v -= v_offset;
#line 535 "IB01CD.f"
    --iwork;
#line 535 "IB01CD.f"
    --dwork;
#line 535 "IB01CD.f"

#line 535 "IB01CD.f"
    /* Function Body */
#line 535 "IB01CD.f"
    withx0 = lsame_(jobx0, "X", (ftnlen)1, (ftnlen)1);
#line 536 "IB01CD.f"
    compbd = lsame_(comuse, "C", (ftnlen)1, (ftnlen)1);
#line 537 "IB01CD.f"
    usebd = lsame_(comuse, "U", (ftnlen)1, (ftnlen)1);
#line 538 "IB01CD.f"
    withd = lsame_(job, "D", (ftnlen)1, (ftnlen)1);
#line 539 "IB01CD.f"
    withb = lsame_(job, "B", (ftnlen)1, (ftnlen)1) || withd;
#line 540 "IB01CD.f"
    maxdim = withx0 && usebd || compbd;
#line 541 "IB01CD.f"
    maxdia = withx0 || compbd;

#line 543 "IB01CD.f"
    *iwarn = 0;
#line 544 "IB01CD.f"
    *info = 0;
#line 545 "IB01CD.f"
    ldw = max(1,*n);
#line 546 "IB01CD.f"
    lm = *l * *m;
#line 547 "IB01CD.f"
    ln = *l * *n;
#line 548 "IB01CD.f"
    nn = *n * *n;
#line 549 "IB01CD.f"
    nm = *n * *m;
#line 550 "IB01CD.f"
    n2m = *n * nm;
#line 551 "IB01CD.f"
    if (compbd) {
#line 552 "IB01CD.f"
	ncol = nm;
#line 553 "IB01CD.f"
	if (withx0) {
#line 553 "IB01CD.f"
	    ncol += *n;
#line 553 "IB01CD.f"
	}
#line 555 "IB01CD.f"
	minsmp = ncol;
#line 556 "IB01CD.f"
	if (withd) {
#line 557 "IB01CD.f"
	    minsmp += *m;
#line 558 "IB01CD.f"
	    iq = minsmp;
#line 559 "IB01CD.f"
	} else if (! withx0) {
#line 560 "IB01CD.f"
	    iq = minsmp;
#line 561 "IB01CD.f"
	    ++minsmp;
#line 562 "IB01CD.f"
	} else {
#line 563 "IB01CD.f"
	    iq = minsmp;
#line 564 "IB01CD.f"
	}
#line 565 "IB01CD.f"
    } else {
#line 566 "IB01CD.f"
	ncol = *n;
#line 567 "IB01CD.f"
	if (withx0) {
#line 568 "IB01CD.f"
	    minsmp = *n;
#line 569 "IB01CD.f"
	} else {
#line 570 "IB01CD.f"
	    minsmp = 0;
#line 571 "IB01CD.f"
	}
#line 572 "IB01CD.f"
	iq = minsmp;
#line 573 "IB01CD.f"
    }

#line 575 "IB01CD.f"
    if (! (withx0 || lsame_(jobx0, "N", (ftnlen)1, (ftnlen)1))) {
#line 576 "IB01CD.f"
	*info = -1;
#line 577 "IB01CD.f"
    } else if (! (compbd || usebd || lsame_(comuse, "N", (ftnlen)1, (ftnlen)1)
	    )) {
#line 579 "IB01CD.f"
	*info = -2;
#line 580 "IB01CD.f"
    } else if (! withb) {
#line 581 "IB01CD.f"
	*info = -3;
#line 582 "IB01CD.f"
    } else if (*n < 0) {
#line 583 "IB01CD.f"
	*info = -4;
#line 584 "IB01CD.f"
    } else if (*m < 0) {
#line 585 "IB01CD.f"
	*info = -5;
#line 586 "IB01CD.f"
    } else if (*l <= 0) {
#line 587 "IB01CD.f"
	*info = -6;
#line 588 "IB01CD.f"
    } else if (*nsmp < minsmp) {
#line 589 "IB01CD.f"
	*info = -7;
#line 590 "IB01CD.f"
    } else if (*lda < 1 || maxdia && *lda < ldw) {
#line 591 "IB01CD.f"
	*info = -9;
#line 592 "IB01CD.f"
    } else if (*ldb < 1 || *m > 0 && maxdim && *ldb < ldw) {
#line 594 "IB01CD.f"
	*info = -11;
#line 595 "IB01CD.f"
    } else if (*ldc < 1 || *n > 0 && maxdia && *ldc < *l) {
#line 597 "IB01CD.f"
	*info = -13;
#line 598 "IB01CD.f"
    } else if (*ldd < 1 || *m > 0 && maxdim && withd && *ldd < *l) {
#line 600 "IB01CD.f"
	*info = -15;
#line 601 "IB01CD.f"
    } else if (*ldu < 1 || *m > 0 && maxdim && *ldu < *nsmp) {
#line 603 "IB01CD.f"
	*info = -17;
#line 604 "IB01CD.f"
    } else if (*ldy < 1 || maxdia && *ldy < *nsmp) {
#line 605 "IB01CD.f"
	*info = -19;
#line 606 "IB01CD.f"
    } else if (*ldv < 1 || maxdia && *ldv < ldw) {
#line 607 "IB01CD.f"
	*info = -22;
#line 608 "IB01CD.f"
    } else if (*tol > 1.) {
#line 609 "IB01CD.f"
	*info = -23;
#line 610 "IB01CD.f"
    }

/*     Compute workspace. */
/*      (Note: Comments in the code beginning "Workspace:" describe the */
/*       minimal amount of workspace needed at that point in the code, */
/*       as well as the preferred amount for good performance. */
/*       NB refers to the optimal block size for the immediately */
/*       following subroutine, as returned by ILAENV.) */

#line 619 "IB01CD.f"
    if (! maxdia || max(*n,*m) == 0) {
#line 620 "IB01CD.f"
	minwrk = 2;
#line 621 "IB01CD.f"
    } else {
#line 622 "IB01CD.f"
	nsmpl = *nsmp * *l;
#line 623 "IB01CD.f"
	iq *= *l;
#line 624 "IB01CD.f"
	ncp1 = ncol + 1;
#line 625 "IB01CD.f"
	isize = nsmpl * ncp1;
#line 626 "IB01CD.f"
	if (compbd) {
#line 627 "IB01CD.f"
	    if (*n > 0 && withx0) {
#line 628 "IB01CD.f"
		ic = (nn << 1) + *n;
#line 629 "IB01CD.f"
	    } else {
#line 630 "IB01CD.f"
		ic = 0;
#line 631 "IB01CD.f"
	    }
#line 632 "IB01CD.f"
	} else {
#line 633 "IB01CD.f"
	    ic = nn << 1;
#line 634 "IB01CD.f"
	}
#line 635 "IB01CD.f"
	minwls = ncol * ncp1;
#line 636 "IB01CD.f"
	if (compbd) {
#line 637 "IB01CD.f"
	    if (withd) {
#line 637 "IB01CD.f"
		minwls += lm * ncp1;
#line 637 "IB01CD.f"
	    }
#line 639 "IB01CD.f"
	    if (*m > 0 && withd) {
/* Computing MAX */
#line 640 "IB01CD.f"
		i__1 = ncol << 1;
#line 640 "IB01CD.f"
		ia = *m + max(i__1,*m);
#line 641 "IB01CD.f"
	    } else {
#line 642 "IB01CD.f"
		ia = ncol << 1;
#line 643 "IB01CD.f"
	    }
#line 644 "IB01CD.f"
	    itau = n2m + max(ic,ia);
#line 645 "IB01CD.f"
	    if (withx0) {
#line 645 "IB01CD.f"
		itau += ln;
#line 645 "IB01CD.f"
	    }
/* Computing MAX */
#line 647 "IB01CD.f"
	    i__1 = *n + max(ic,ia), i__2 = ncol * 6;
#line 647 "IB01CD.f"
	    ldw2 = isize + max(i__1,i__2);
/* Computing MAX */
#line 648 "IB01CD.f"
	    i__1 = iq * ncp1 + itau, i__2 = ncol * 6;
#line 648 "IB01CD.f"
	    ldw3 = minwls + max(i__1,i__2);
#line 649 "IB01CD.f"
	    if (*m > 0 && withd) {
/* Computing MAX */
#line 650 "IB01CD.f"
		i__1 = ldw2, i__2 = isize + (*m << 1) * *m + *m * 6;
#line 650 "IB01CD.f"
		ldw2 = max(i__1,i__2);
/* Computing MAX */
#line 651 "IB01CD.f"
		i__1 = ldw3, i__2 = minwls + (*m << 1) * *m + *m * 6;
#line 651 "IB01CD.f"
		ldw3 = max(i__1,i__2);
#line 652 "IB01CD.f"
		ia = 3;
#line 653 "IB01CD.f"
	    } else {
#line 654 "IB01CD.f"
		ia = 2;
#line 655 "IB01CD.f"
	    }
#line 656 "IB01CD.f"
	} else {
#line 657 "IB01CD.f"
	    itau = ic + ln;
/* Computing MAX */
#line 658 "IB01CD.f"
	    i__1 = ic, i__2 = *n << 2;
#line 658 "IB01CD.f"
	    ldw2 = isize + (*n << 1) + max(i__1,i__2);
/* Computing MAX */
#line 659 "IB01CD.f"
	    i__1 = iq * ncp1 + itau, i__2 = *n << 2;
#line 659 "IB01CD.f"
	    ldw3 = minwls + (*n << 1) + max(i__1,i__2);
#line 660 "IB01CD.f"
	    ia = 2;
#line 661 "IB01CD.f"
	}
/* Computing MAX */
#line 662 "IB01CD.f"
	i__1 = *n * 5, i__1 = max(i__1,ia), i__2 = min(ldw2,ldw3);
#line 662 "IB01CD.f"
	minwrk = ia + nn + nm + ln + max(i__1,i__2);

#line 664 "IB01CD.f"
	if (*info == 0 && *ldwork >= minwrk) {
/* Computing MAX */
#line 665 "IB01CD.f"
	    i__1 = *n * 5;
#line 665 "IB01CD.f"
	    maxwrk = max(i__1,ia);
#line 666 "IB01CD.f"
	    if (compbd) {
#line 667 "IB01CD.f"
		if (*m > 0 && withd) {
/* Computing MAX */
/* Computing MAX */
#line 668 "IB01CD.f"
		    i__5 = *nsmp - *m;
#line 668 "IB01CD.f"
		    i__3 = *m * ilaenv_(&c__1, "DGEQRF", " ", nsmp, m, &c_n1, 
			    &c_n1, (ftnlen)6, (ftnlen)1), i__4 = ncol + ncol *
			     ilaenv_(&c__1, "DGEQRF", " ", &i__5, &ncol, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1);
#line 668 "IB01CD.f"
		    i__1 = maxwrk, i__2 = isize + *n + *m + max(i__3,i__4);
#line 668 "IB01CD.f"
		    maxwrk = max(i__1,i__2);
/* Computing MAX */
/* Computing MAX */
#line 673 "IB01CD.f"
		    i__5 = *nsmp - *m;
#line 673 "IB01CD.f"
		    i__3 = ncp1 * ilaenv_(&c__1, "DORMQR", "LT", nsmp, &ncp1, 
			    m, &c_n1, (ftnlen)6, (ftnlen)2), i__4 = ncol + 
			    ilaenv_(&c__1, "DORMQR", "LT", &i__5, &c__1, &
			    ncol, &c_n1, (ftnlen)6, (ftnlen)2);
#line 673 "IB01CD.f"
		    i__1 = maxwrk, i__2 = isize + *n + *m + max(i__3,i__4);
#line 673 "IB01CD.f"
		    maxwrk = max(i__1,i__2);
#line 678 "IB01CD.f"
		} else {
/* Computing MAX */
/* Computing MAX */
#line 679 "IB01CD.f"
		    i__3 = ncol * ilaenv_(&c__1, "DGEQRF", " ", &nsmpl, &ncol,
			     &c_n1, &c_n1, (ftnlen)6, (ftnlen)1), i__4 = 
			    ilaenv_(&c__1, "DORMQR", "LT", &nsmpl, &c__1, &
			    ncol, &c_n1, (ftnlen)6, (ftnlen)2);
#line 679 "IB01CD.f"
		    i__1 = maxwrk, i__2 = isize + *n + ncol + max(i__3,i__4);
#line 679 "IB01CD.f"
		    maxwrk = max(i__1,i__2);
#line 684 "IB01CD.f"
		}
#line 685 "IB01CD.f"
	    } else {
/* Computing MAX */
/* Computing MAX */
#line 686 "IB01CD.f"
		i__3 = *n * ilaenv_(&c__1, "DGEQRF", " ", &nsmpl, n, &c_n1, &
			c_n1, (ftnlen)6, (ftnlen)1), i__4 = ilaenv_(&c__1, 
			"DORMQR", "LT", &nsmpl, &c__1, n, &c_n1, (ftnlen)6, (
			ftnlen)2);
#line 686 "IB01CD.f"
		i__1 = maxwrk, i__2 = isize + (*n << 1) + max(i__3,i__4);
#line 686 "IB01CD.f"
		maxwrk = max(i__1,i__2);
#line 691 "IB01CD.f"
	    }
#line 692 "IB01CD.f"
	    maxwrk = ia + nn + nm + ln + maxwrk;
#line 693 "IB01CD.f"
	    maxwrk = max(maxwrk,minwrk);
#line 694 "IB01CD.f"
	}
#line 695 "IB01CD.f"
    }

#line 697 "IB01CD.f"
    if (*info == 0 && *ldwork < minwrk) {
#line 698 "IB01CD.f"
	*info = -26;
#line 699 "IB01CD.f"
	dwork[1] = (doublereal) minwrk;
#line 700 "IB01CD.f"
    }

/*     Return if there are illegal arguments. */

#line 704 "IB01CD.f"
    if (*info != 0) {
#line 705 "IB01CD.f"
	i__1 = -(*info);
#line 705 "IB01CD.f"
	xerbla_("IB01CD", &i__1, (ftnlen)6);
#line 706 "IB01CD.f"
	return 0;
#line 707 "IB01CD.f"
    }

/*     Quick return if possible. */

#line 711 "IB01CD.f"
    if (! maxdia || max(*n,*m) == 0) {
#line 712 "IB01CD.f"
	dwork[2] = 1.;
#line 713 "IB01CD.f"
	if (compbd && *m > 0 && withd) {
#line 714 "IB01CD.f"
	    dwork[1] = 3.;
#line 715 "IB01CD.f"
	    dwork[3] = 1.;
#line 716 "IB01CD.f"
	} else {
#line 717 "IB01CD.f"
	    dwork[1] = 2.;
#line 718 "IB01CD.f"
	}
#line 719 "IB01CD.f"
	if (*n > 0 && usebd) {
#line 720 "IB01CD.f"
	    dum[0] = 0.;
#line 721 "IB01CD.f"
	    dcopy_(n, dum, &c__0, &x0[1], &c__1);
#line 722 "IB01CD.f"
	}
#line 723 "IB01CD.f"
	return 0;
#line 724 "IB01CD.f"
    }

/*     Compute the Schur factorization of  A  and transform the other */
/*     given system matrices accordingly. */
/*     Workspace:  need   g + N*N + L*N + N*M + 5*N,  where */
/*                        g = 2,  if  M = 0, COMUSE = 'C', or  JOB = 'B', */
/*                        g = 3,  if  M > 0, COMUSE = 'C', and JOB = 'D', */
/*                        g = 2,  if  JOBX0 = 'X'  and  COMUSE <> 'C'; */
/*                 prefer larger. */

#line 734 "IB01CD.f"
    ++ia;
#line 735 "IB01CD.f"
    ic = ia + nn;
#line 736 "IB01CD.f"
    ib = ic + ln;
#line 737 "IB01CD.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[ia], &ldw, (ftnlen)4);
#line 738 "IB01CD.f"
    dlacpy_("Full", l, n, &c__[c_offset], ldc, &dwork[ic], l, (ftnlen)4);

#line 740 "IB01CD.f"
    if (usebd) {
#line 741 "IB01CD.f"
	mtmp = *m;
#line 742 "IB01CD.f"
	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ib], &ldw, (ftnlen)4);
#line 743 "IB01CD.f"
    } else {
#line 744 "IB01CD.f"
	mtmp = 0;
#line 745 "IB01CD.f"
    }
#line 746 "IB01CD.f"
    iwr = ib + nm;
#line 747 "IB01CD.f"
    iwi = iwr + *n;
#line 748 "IB01CD.f"
    jwork = iwi + *n;

#line 750 "IB01CD.f"
    i__1 = *ldwork - jwork + 1;
#line 750 "IB01CD.f"
    tb01wd_(n, &mtmp, l, &dwork[ia], &ldw, &dwork[ib], &ldw, &dwork[ic], l, &
	    v[v_offset], ldv, &dwork[iwr], &dwork[iwi], &dwork[jwork], &i__1, 
	    &ierr);
#line 753 "IB01CD.f"
    if (ierr > 0) {
#line 754 "IB01CD.f"
	*info = 1;
#line 755 "IB01CD.f"
	return 0;
#line 756 "IB01CD.f"
    }
/* Computing MAX */
#line 757 "IB01CD.f"
    i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 757 "IB01CD.f"
    maxwrk = max(i__1,i__2);

#line 759 "IB01CD.f"
    i__1 = iwi - 1;
#line 759 "IB01CD.f"
    for (i__ = iwr; i__ <= i__1; ++i__) {
#line 760 "IB01CD.f"
	if (dlapy2_(&dwork[i__], &dwork[i__ + *n]) >= 1.) {
#line 760 "IB01CD.f"
	    *iwarn = 6;
#line 760 "IB01CD.f"
	}
#line 762 "IB01CD.f"
/* L10: */
#line 762 "IB01CD.f"
    }

#line 764 "IB01CD.f"
    jwork = iwr;

/*     Estimate  x(0)  and/or the system matrices  B and D. */
/*     Workspace: need   g + N*N + L*N + N*M + */
/*                           max( g, min( LDW2, LDW3 ) ) (see LDWORK); */
/*                prefer larger. */

#line 771 "IB01CD.f"
    if (compbd) {
#line 772 "IB01CD.f"
	i__1 = *ldwork - jwork + 1;
#line 772 "IB01CD.f"
	ib01qd_(jobx0, job, n, m, l, nsmp, &dwork[ia], &ldw, &dwork[ic], l, &
		u[u_offset], ldu, &y[y_offset], ldy, &x0[1], &dwork[ib], &ldw,
		 &d__[d_offset], ldd, tol, &iwork[1], &dwork[jwork], &i__1, &
		iwarnl, info, (ftnlen)1, (ftnlen)1);

#line 777 "IB01CD.f"
	if (*info == 0) {
#line 778 "IB01CD.f"
	    if (*m > 0 && withd) {
#line 778 "IB01CD.f"
		rcondu = dwork[jwork + 2];
#line 778 "IB01CD.f"
	    }

/*           Compute the system input matrix  B  corresponding to the */
/*           original system. */

#line 784 "IB01CD.f"
	    dgemm_("NoTranspose", "NoTranspose", n, m, n, &c_b57, &v[v_offset]
		    , ldv, &dwork[ib], &ldw, &c_b58, &b[b_offset], ldb, (
		    ftnlen)11, (ftnlen)11);
#line 786 "IB01CD.f"
	}
#line 787 "IB01CD.f"
    } else {
#line 788 "IB01CD.f"
	if (withd) {
#line 789 "IB01CD.f"
	    *(unsigned char *)jobd = 'N';
#line 790 "IB01CD.f"
	} else {
#line 791 "IB01CD.f"
	    *(unsigned char *)jobd = 'Z';
#line 792 "IB01CD.f"
	}

#line 794 "IB01CD.f"
	i__1 = *ldwork - jwork + 1;
#line 794 "IB01CD.f"
	ib01rd_(jobd, n, &mtmp, l, nsmp, &dwork[ia], &ldw, &dwork[ib], &ldw, &
		dwork[ic], l, &d__[d_offset], ldd, &u[u_offset], ldu, &y[
		y_offset], ldy, &x0[1], tol, &iwork[1], &dwork[jwork], &i__1, 
		&iwarnl, info, (ftnlen)1);
#line 798 "IB01CD.f"
    }
#line 799 "IB01CD.f"
    *iwarn = max(*iwarn,iwarnl);

#line 801 "IB01CD.f"
    if (*info == 0) {
#line 802 "IB01CD.f"
	rcond = dwork[jwork + 1];
/* Computing MAX */
#line 803 "IB01CD.f"
	i__1 = maxwrk, i__2 = (integer) dwork[jwork] + jwork - 1;
#line 803 "IB01CD.f"
	maxwrk = max(i__1,i__2);
#line 804 "IB01CD.f"
	if (withx0) {

/*           Transform the initial state estimate to obtain the initial */
/*           state corresponding to the original system. */
/*           Workspace: need g + N*N + L*N + N*M + N. */

#line 810 "IB01CD.f"
	    dgemv_("NoTranspose", n, n, &c_b57, &v[v_offset], ldv, &x0[1], &
		    c__1, &c_b58, &dwork[jwork], &c__1, (ftnlen)11);
#line 812 "IB01CD.f"
	    dcopy_(n, &dwork[jwork], &c__1, &x0[1], &c__1);
#line 813 "IB01CD.f"
	}

#line 815 "IB01CD.f"
	dwork[1] = (doublereal) maxwrk;
#line 816 "IB01CD.f"
	dwork[2] = rcond;
#line 817 "IB01CD.f"
	if (compbd && *m > 0 && withd) {
#line 817 "IB01CD.f"
	    dwork[3] = rcondu;
#line 817 "IB01CD.f"
	}
#line 819 "IB01CD.f"
    }
#line 820 "IB01CD.f"
    return 0;

/* *** End of IB01CD *** */
} /* ib01cd_ */

