#line 1 "MD03BD.f"
/* MD03BD.f -- translated by f2c (version 20100827).
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

#line 1 "MD03BD.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int md03bd_(char *xinit, char *scale, char *cond, S_fp fcn, 
	S_fp qrfact, S_fp lmparm, integer *m, integer *n, integer *itmax, 
	doublereal *factor, integer *nprint, integer *ipar, integer *lipar, 
	doublereal *dpar1, integer *ldpar1, doublereal *dpar2, integer *
	ldpar2, doublereal *x, doublereal *diag, integer *nfev, integer *njev,
	 doublereal *ftol, doublereal *xtol, doublereal *gtol, doublereal *
	tol, integer *iwork, doublereal *dwork, integer *ldwork, integer *
	iwarn, integer *info, ftnlen xinit_len, ftnlen scale_len, ftnlen 
	cond_len)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7, i__8;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer e, j, l, nc, iw1, iw2, iw3, jw1, jw2, jac, ldj;
    static doublereal par;
    static integer seed[4];
    static logical init;
    static integer iter, llmp, lqrf;
    static doublereal temp;
    static integer lfcn1, lfcn2;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static doublereal temp1, temp2;
    static integer iflag;
    static doublereal ftdef, delta, gtdef;
    static logical iscal;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical sscal;
    static integer infol, nfevl;
    static doublereal xtdef, ratio;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal fnorm, gnorm;
    static integer sizej, jwork;
    static doublereal pnorm, xnorm, fnorm1;
    static logical badscl;
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal actred, dirder, toldef, epsmch, prered;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer ldjsav;
    extern /* Subroutine */ int dlarnv_(integer *, integer *, integer *, 
	    doublereal *);
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

/*     To minimize the sum of the squares of m nonlinear functions, e, in */
/*     n variables, x, by a modification of the Levenberg-Marquardt */
/*     algorithm. The user must provide a subroutine FCN which calculates */
/*     the functions and the Jacobian (possibly by finite differences). */
/*     In addition, specialized subroutines QRFACT, for QR factorization */
/*     with pivoting of the Jacobian, and LMPARM, for the computation of */
/*     Levenberg-Marquardt parameter, exploiting the possible structure */
/*     of the Jacobian matrix, should be provided. Template */
/*     implementations of these routines are included in SLICOT Library. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     XINIT   CHARACTER*1 */
/*             Specifies how the variables x are initialized, as follows: */
/*             = 'R' :  the array X is initialized to random values; the */
/*                      entries DWORK(1:4) are used to initialize the */
/*                      random number generator: the first three values */
/*                      are converted to integers between 0 and 4095, and */
/*                      the last one is converted to an odd integer */
/*                      between 1 and 4095; */
/*             = 'G' :  the given entries of X are used as initial values */
/*                      of variables. */

/*     SCALE   CHARACTER*1 */
/*             Specifies how the variables will be scaled, as follows: */
/*             = 'I' :  use internal scaling; */
/*             = 'S' :  use specified scaling factors, given in DIAG. */

/*     COND    CHARACTER*1 */
/*             Specifies whether the condition of the linear systems */
/*             involved should be estimated, as follows: */
/*             = 'E' :  use incremental condition estimation to find the */
/*                      numerical rank; */
/*             = 'N' :  do not use condition estimation, but check the */
/*                      diagonal entries of matrices for zero values. */

/*     Function Parameters */

/*     FCN     EXTERNAL */
/*             Subroutine which evaluates the functions and the Jacobian. */
/*             FCN must be declared in an external statement in the user */
/*             calling program, and must have the following interface: */

/*             SUBROUTINE FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, */
/*            $                DPAR2, LDPAR2, X, NFEVL, E, J, LDJ, DWORK, */
/*            $                LDWORK, INFO ) */

/*             where */

/*             IFLAG   (input/output) INTEGER */
/*                     On entry, this parameter must contain a value */
/*                     defining the computations to be performed: */
/*                     = 0 :  Optionally, print the current iterate X, */
/*                            function values E, and Jacobian matrix J, */
/*                            or other results defined in terms of these */
/*                            values. See the argument NPRINT of MD03BD. */
/*                            Do not alter E and J. */
/*                     = 1 :  Calculate the functions at X and return */
/*                            this vector in E. Do not alter J. */
/*                     = 2 :  Calculate the Jacobian at X and return */
/*                            this matrix in J. Also return NFEVL */
/*                            (see below). Do not alter E. */
/*                     = 3 :  Do not compute neither the functions nor */
/*                            the Jacobian, but return in LDJ and */
/*                            IPAR/DPAR1,DPAR2 (some of) the integer/real */
/*                            parameters needed. */
/*                     On exit, the value of this parameter should not be */
/*                     changed by FCN unless the user wants to terminate */
/*                     execution of MD03BD, in which case IFLAG must be */
/*                     set to a negative integer. */

/*             M       (input) INTEGER */
/*                     The number of functions.  M >= 0. */

/*             N       (input) INTEGER */
/*                     The number of variables.  M >= N >= 0. */

/*             IPAR    (input/output) INTEGER array, dimension (LIPAR) */
/*                     The integer parameters describing the structure of */
/*                     the Jacobian matrix or needed for problem solving. */
/*                     IPAR is an input parameter, except for IFLAG = 3 */
/*                     on entry, when it is also an output parameter. */
/*                     On exit, if IFLAG = 3, IPAR(1) contains the length */
/*                     of the array J, for storing the Jacobian matrix, */
/*                     and the entries IPAR(2:5) contain the workspace */
/*                     required by FCN for IFLAG = 1, FCN for IFLAG = 2, */
/*                     QRFACT, and LMPARM, respectively. */

/*             LIPAR   (input) INTEGER */
/*                     The length of the array IPAR.  LIPAR >= 5. */

/*             DPAR1   (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDPAR1,*) or (LDPAR1) */
/*                     A first set of real parameters needed for */
/*                     describing or solving the problem. */
/*                     DPAR1 can also be used as an additional array for */
/*                     intermediate results when computing the functions */
/*                     or the Jacobian. For control problems, DPAR1 could */
/*                     store the input trajectory of a system. */

/*             LDPAR1  (input) INTEGER */
/*                     The leading dimension or the length of the array */
/*                     DPAR1, as convenient.  LDPAR1 >= 0.  (LDPAR1 >= 1, */
/*                     if leading dimension.) */

/*             DPAR2   (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDPAR2,*) or (LDPAR2) */
/*                     A second set of real parameters needed for */
/*                     describing or solving the problem. */
/*                     DPAR2 can also be used as an additional array for */
/*                     intermediate results when computing the functions */
/*                     or the Jacobian. For control problems, DPAR2 could */
/*                     store the output trajectory of a system. */

/*             LDPAR2  (input) INTEGER */
/*                     The leading dimension or the length of the array */
/*                     DPAR2, as convenient.  LDPAR2 >= 0.  (LDPAR2 >= 1, */
/*                     if leading dimension.) */

/*             X       (input) DOUBLE PRECISION array, dimension (N) */
/*                     This array must contain the value of the */
/*                     variables x where the functions or the Jacobian */
/*                     must be evaluated. */

/*             NFEVL   (input/output) INTEGER */
/*                     The number of function evaluations needed to */
/*                     compute the Jacobian by a finite difference */
/*                     approximation. */
/*                     NFEVL is an input parameter if IFLAG = 0, or an */
/*                     output parameter if IFLAG = 2. If the Jacobian is */
/*                     computed analytically, NFEVL should be set to a */
/*                     non-positive value. */

/*             E       (input/output) DOUBLE PRECISION array, */
/*                     dimension (M) */
/*                     This array contains the value of the (error) */
/*                     functions e evaluated at X. */
/*                     E is an input parameter if IFLAG = 0 or 2, or an */
/*                     output parameter if IFLAG = 1. */

/*             J       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDJ,NC), where NC is the number of columns */
/*                     needed. */
/*                     This array contains a possibly compressed */
/*                     representation of the Jacobian matrix evaluated */
/*                     at X. If full Jacobian is stored, then NC = N. */
/*                     J is an input parameter if IFLAG = 0, or an output */
/*                     parameter if IFLAG = 2. */

/*             LDJ     (input/output) INTEGER */
/*                     The leading dimension of array J.  LDJ >= 1. */
/*                     LDJ is essentially used inside the routines FCN, */
/*                     QRFACT and LMPARM. */
/*                     LDJ is an input parameter, except for IFLAG = 3 */
/*                     on entry, when it is an output parameter. */
/*                     It is assumed in MD03BD that LDJ is not larger */
/*                     than needed. */

/*             DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*                     The workspace array for subroutine FCN. */
/*                     On exit, if INFO = 0, DWORK(1) returns the optimal */
/*                     value of LDWORK. */

/*             LDWORK  (input) INTEGER */
/*                     The size of the array DWORK (as large as needed */
/*                     in the subroutine FCN).  LDWORK >= 1. */

/*             INFO    INTEGER */
/*                     Error indicator, set to a negative value if an */
/*                     input (scalar) argument is erroneous, and to */
/*                     positive values for other possible errors in the */
/*                     subroutine FCN. The LAPACK Library routine XERBLA */
/*                     should be used in conjunction with negative INFO. */
/*                     INFO must be zero if the subroutine finished */
/*                     successfully. */

/*             Parameters marked with "(input)" must not be changed. */

/*     QRFACT  EXTERNAL */
/*             Subroutine which computes the QR factorization with */
/*             (block) column pivoting of the Jacobian matrix, J*P = Q*R. */
/*             QRFACT must be declared in an external statement in the */
/*             calling program, and must have the following interface: */

/*             SUBROUTINE QRFACT( N, IPAR, LIPAR, FNORM, J, LDJ, E, */
/*            $                   JNORMS, GNORM, IPVT, DWORK, LDWORK, */
/*            $                   INFO ) */

/*             where */

/*             N       (input) INTEGER */
/*                     The number of columns of the Jacobian matrix J. */
/*                     N >= 0. */

/*             IPAR    (input) INTEGER array, dimension (LIPAR) */
/*                     The integer parameters describing the structure of */
/*                     the Jacobian matrix. */

/*             LIPAR   (input) INTEGER */
/*                     The length of the array IPAR.  LIPAR >= 0. */

/*             FNORM   (input) DOUBLE PRECISION */
/*                     The Euclidean norm of the vector e.  FNORM >= 0. */

/*             J       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDJ, NC), where NC is the number of columns. */
/*                     On entry, the leading NR-by-NC part of this array */
/*                     must contain the (compressed) representation */
/*                     of the Jacobian matrix J, where NR is the number */
/*                     of rows of J (function of IPAR entries). */
/*                     On exit, the leading N-by-NC part of this array */
/*                     contains a (compressed) representation of the */
/*                     upper triangular factor R of the Jacobian matrix. */
/*                     For efficiency of the later calculations, the */
/*                     matrix R is delivered with the leading dimension */
/*                     MAX(1,N), possibly much smaller than the value */
/*                     of LDJ on entry. */

/*             LDJ     (input/output) INTEGER */
/*                     The leading dimension of array J. */
/*                     On entry, LDJ >= MAX(1,NR). */
/*                     On exit,  LDJ >= MAX(1,N). */

/*             E       (input/output) DOUBLE PRECISION array, dimension */
/*                     (NR) */
/*                     On entry, this array contains the error vector e. */
/*                     On exit, this array contains the updated vector */
/*                     Z*Q'*e, where Z is a block row permutation matrix */
/*                     (possibly identity) used in the QR factorization */
/*                     of J. (See, for example, the SLICOT Library */
/*                     routine NF01BS, Section METHOD.) */

/*             JNORMS  (output) DOUBLE PRECISION array, dimension (N) */
/*                     This array contains the Euclidean norms of the */
/*                     columns of the Jacobian matrix (in the original */
/*                     order). */

/*             GNORM   (output) DOUBLE PRECISION */
/*                     If FNORM > 0, the 1-norm of the scaled vector */
/*                     J'*e/FNORM, with each element i further divided */
/*                     by JNORMS(i) (if JNORMS(i) is nonzero). */
/*                     If FNORM = 0, the returned value of GNORM is 0. */

/*             IPVT    (output) INTEGER array, dimension (N) */
/*                     This array defines the permutation matrix P such */
/*                     that J*P = Q*R. Column j of P is column IPVT(j) of */
/*                     the identity matrix. */

/*             DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*                     The workspace array for subroutine QRFACT. */
/*                     On exit, if INFO = 0, DWORK(1) returns the optimal */
/*                     value of LDWORK. */

/*             LDWORK  (input) INTEGER */
/*                     The size of the array DWORK (as large as needed */
/*                     in the subroutine QRFACT).  LDWORK >= 1. */

/*             INFO    INTEGER */
/*                     Error indicator, set to a negative value if an */
/*                     input (scalar) argument is erroneous, and to */
/*                     positive values for other possible errors in the */
/*                     subroutine QRFACT. The LAPACK Library routine */
/*                     XERBLA should be used in conjunction with negative */
/*                     INFO. INFO must be zero if the subroutine finished */
/*                     successfully. */

/*             Parameters marked with "(input)" must not be changed. */

/*     LMPARM  EXTERNAL */
/*             Subroutine which determines a value for the Levenberg- */
/*             Marquardt parameter PAR such that if x solves the system */

/*                   J*x = b ,     sqrt(PAR)*D*x = 0 , */

/*             in the least squares sense, where J is an m-by-n matrix, */
/*             D is an n-by-n nonsingular diagonal matrix, and b is an */
/*             m-vector, and if DELTA is a positive number, DXNORM is */
/*             the Euclidean norm of D*x, then either PAR is zero and */

/*                   ( DXNORM - DELTA ) .LE. 0.1*DELTA , */

/*             or PAR is positive and */

/*                   ABS( DXNORM - DELTA ) .LE. 0.1*DELTA . */

/*             It is assumed that a block QR factorization, with column */
/*             pivoting, of J is available, that is, J*P = Q*R, where P */
/*             is a permutation matrix, Q has orthogonal columns, and */
/*             R is an upper triangular matrix (possibly stored in a */
/*             compressed form), with diagonal elements of nonincreasing */
/*             magnitude for each block. On output, LMPARM also provides */
/*             a (compressed) representation of an upper triangular */
/*             matrix S, such that */

/*                   P'*(J'*J + PAR*D*D)*P = S'*S . */

/*             LMPARM must be declared in an external statement in the */
/*             calling program, and must have the following interface: */

/*             SUBROUTINE LMPARM( COND, N, IPAR, LIPAR, R, LDR, IPVT, */
/*            $                   DIAG, QTB, DELTA, PAR, RANKS, X, RX, */
/*            $                   TOL, DWORK, LDWORK, INFO ) */

/*             where */

/*             COND    CHARACTER*1 */
/*                     Specifies whether the condition of the linear */
/*                     systems involved should be estimated, as follows: */
/*                     = 'E' :  use incremental condition estimation */
/*                              to find the numerical rank; */
/*                     = 'N' :  do not use condition estimation, but */
/*                              check the diagonal entries for zero */
/*                              values; */
/*                     = 'U' :  use the ranks already stored in RANKS */
/*                              (for R). */

/*             N       (input) INTEGER */
/*                     The order of the matrix R.  N >= 0. */

/*             IPAR    (input) INTEGER array, dimension (LIPAR) */
/*                     The integer parameters describing the structure of */
/*                     the Jacobian matrix. */

/*             LIPAR   (input) INTEGER */
/*                     The length of the array IPAR.  LIPAR >= 0. */

/*             R       (input/output) DOUBLE PRECISION array, dimension */
/*                     (LDR, NC), where NC is the number of columns. */
/*                     On entry, the leading N-by-NC part of this array */
/*                     must contain the (compressed) representation (Rc) */
/*                     of the upper triangular matrix R. */
/*                     On exit, the full upper triangular part of R */
/*                     (in representation Rc), is unaltered, and the */
/*                     remaining part contains (part of) the (compressed) */
/*                     representation of the transpose of the upper */
/*                     triangular matrix S. */

/*             LDR     (input) INTEGER */
/*                     The leading dimension of array R. */
/*                     LDR >= MAX(1,N). */

/*             IPVT    (input) INTEGER array, dimension (N) */
/*                     This array must define the permutation matrix P */
/*                     such that J*P = Q*R. Column j of P is column */
/*                     IPVT(j) of the identity matrix. */

/*             DIAG    (input) DOUBLE PRECISION array, dimension (N) */
/*                     This array must contain the diagonal elements of */
/*                     the matrix D.  DIAG(I) <> 0, I = 1,...,N. */

/*             QTB     (input) DOUBLE PRECISION array, dimension (N) */
/*                     This array must contain the first n elements of */
/*                     the vector Q'*b. */

/*             DELTA   (input) DOUBLE PRECISION */
/*                     An upper bound on the Euclidean norm of D*x. */
/*                     DELTA > 0. */

/*             PAR     (input/output) DOUBLE PRECISION */
/*                     On entry, PAR must contain an initial estimate of */
/*                     the Levenberg-Marquardt parameter.  PAR >= 0. */
/*                     On exit, it contains the final estimate of this */
/*                     parameter. */

/*             RANKS   (input or output) INTEGER array, dimension (r), */
/*                     where r is the number of diagonal blocks R_k in R, */
/*                     corresponding to the block column structure of J. */
/*                     On entry, if COND = 'U' and N > 0, this array must */
/*                     contain the numerical ranks of the submatrices */
/*                     R_k, k = 1:r. The number r is defined in terms of */
/*                     the entries of IPAR. */
/*                     On exit, if N > 0, this array contains the */
/*                     numerical ranks of the submatrices S_k, k = 1:r. */

/*             X       (output) DOUBLE PRECISION array, dimension (N) */
/*                     This array contains the least squares solution of */
/*                     the system J*x = b, sqrt(PAR)*D*x = 0. */

/*             RX      (output) DOUBLE PRECISION array, dimension (N) */
/*                     This array contains the matrix-vector product */
/*                     -R*P'*x. */

/*             TOL     (input) DOUBLE PRECISION */
/*                     If COND = 'E', the tolerance to be used for */
/*                     finding the ranks of the submatrices R_k and S_k. */
/*                     If the user sets TOL > 0, then the given value of */
/*                     TOL is used as a lower bound for the reciprocal */
/*                     condition number;  a (sub)matrix whose estimated */
/*                     condition number is less than 1/TOL is considered */
/*                     to be of full rank.  If the user sets TOL <= 0, */
/*                     then an implicitly computed, default tolerance, */
/*                     defined by TOLDEF = N*EPS,  is used instead, */
/*                     where EPS is the machine precision (see LAPACK */
/*                     Library routine DLAMCH). */
/*                     This parameter is not relevant if COND = 'U' */
/*                     or 'N'. */

/*             DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*                     The workspace array for subroutine LMPARM. */
/*                     On exit, if INFO = 0, DWORK(1) returns the optimal */
/*                     value of LDWORK. */

/*             LDWORK  (input) INTEGER */
/*                     The size of the array DWORK (as large as needed */
/*                     in the subroutine LMPARM).  LDWORK >= 1. */

/*             INFO    INTEGER */
/*                     Error indicator, set to a negative value if an */
/*                     input (scalar) argument is erroneous, and to */
/*                     positive values for other possible errors in the */
/*                     subroutine LMPARM. The LAPACK Library routine */
/*                     XERBLA should be used in conjunction with negative */
/*                     INFO. INFO must be zero if the subroutine finished */
/*                     successfully. */

/*             Parameters marked with "(input)" must not be changed. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of functions.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of variables.  M >= N >= 0. */

/*     ITMAX   (input) INTEGER */
/*             The maximum number of iterations.  ITMAX >= 0. */

/*     FACTOR  (input) DOUBLE PRECISION */
/*             The value used in determining the initial step bound. This */
/*             bound is set to the product of FACTOR and the Euclidean */
/*             norm of DIAG*X if nonzero, or else to FACTOR itself. */
/*             In most cases FACTOR should lie in the interval (.1,100). */
/*             A generally recommended value is 100.  FACTOR > 0. */

/*     NPRINT  (input) INTEGER */
/*             This parameter enables controlled printing of iterates if */
/*             it is positive. In this case, FCN is called with IFLAG = 0 */
/*             at the beginning of the first iteration and every NPRINT */
/*             iterations thereafter and immediately prior to return, */
/*             with X, E, and J available for printing. Note that when */
/*             called immediately prior to return, J normally contains */
/*             the result returned by QRFACT and LMPARM (the compressed */
/*             R and S factors). If NPRINT is not positive, no special */
/*             calls of FCN with IFLAG = 0 are made. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters needed, for instance, for */
/*             describing the structure of the Jacobian matrix, which */
/*             are handed over to the routines FCN, QRFACT and LMPARM. */
/*             The first five entries of this array are modified */
/*             internally by a call to FCN (with IFLAG = 3), but are */
/*             restored on exit. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 5. */

/*     DPAR1   (input/output) DOUBLE PRECISION array, dimension */
/*             (LDPAR1,*) or (LDPAR1) */
/*             A first set of real parameters needed for describing or */
/*             solving the problem. This argument is not used by MD03BD */
/*             routine, but it is passed to the routine FCN. */

/*     LDPAR1  (input) INTEGER */
/*             The leading dimension or the length of the array DPAR1, as */
/*             convenient.  LDPAR1 >= 0.  (LDPAR1 >= 1, if leading */
/*             dimension.) */

/*     DPAR2   (input/output) DOUBLE PRECISION array, dimension */
/*             (LDPAR2,*) or (LDPAR2) */
/*             A second set of real parameters needed for describing or */
/*             solving the problem. This argument is not used by MD03BD */
/*             routine, but it is passed to the routine FCN. */

/*     LDPAR2  (input) INTEGER */
/*             The leading dimension or the length of the array DPAR2, as */
/*             convenient.  LDPAR2 >= 0.  (LDPAR2 >= 1, if leading */
/*             dimension.) */

/*     X       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, if XINIT = 'G', this array must contain the */
/*             vector of initial variables x to be optimized. */
/*             If XINIT = 'R', this array need not be set before entry, */
/*             and random values will be used to initialize x. */
/*             On exit, if INFO = 0, this array contains the vector of */
/*             values that (approximately) minimize the sum of squares of */
/*             error functions. The values returned in IWARN and */
/*             DWORK(1:4) give details on the iterative process. */

/*     DIAG    (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, if SCALE = 'S', this array must contain some */
/*             positive entries that serve as multiplicative scale */
/*             factors for the variables x.  DIAG(I) > 0, I = 1,...,N. */
/*             If SCALE = 'I', DIAG is internally set. */
/*             On exit, this array contains the scale factors used */
/*             (or finally used, if SCALE = 'I'). */

/*     NFEV    (output) INTEGER */
/*             The number of calls to FCN with IFLAG = 1. If FCN is */
/*             properly implemented, this includes the function */
/*             evaluations needed for finite difference approximation */
/*             of the Jacobian. */

/*     NJEV    (output) INTEGER */
/*             The number of calls to FCN with IFLAG = 2. */

/*     Tolerances */

/*     FTOL    DOUBLE PRECISION */
/*             If FTOL >= 0, the tolerance which measures the relative */
/*             error desired in the sum of squares. Termination occurs */
/*             when both the actual and predicted relative reductions in */
/*             the sum of squares are at most FTOL. If the user sets */
/*             FTOL < 0,  then  SQRT(EPS)  is used instead FTOL, where */
/*             EPS is the machine precision (see LAPACK Library routine */
/*             DLAMCH). */

/*     XTOL    DOUBLE PRECISION */
/*             If XTOL >= 0, the tolerance which measures the relative */
/*             error desired in the approximate solution. Termination */
/*             occurs when the relative error between two consecutive */
/*             iterates is at most XTOL. If the user sets  XTOL < 0, */
/*             then  SQRT(EPS)  is used instead XTOL. */

/*     GTOL    DOUBLE PRECISION */
/*             If GTOL >= 0, the tolerance which measures the */
/*             orthogonality desired between the function vector e and */
/*             the columns of the Jacobian J. Termination occurs when */
/*             the cosine of the angle between e and any column of the */
/*             Jacobian J is at most GTOL in absolute value. If the user */
/*             sets  GTOL < 0,  then  EPS  is used instead GTOL. */

/*     TOL     DOUBLE PRECISION */
/*             If COND = 'E', the tolerance to be used for finding the */
/*             ranks of the matrices of linear systems to be solved. If */
/*             the user sets TOL > 0, then the given value of TOL is used */
/*             as a lower bound for the reciprocal condition number;  a */
/*             (sub)matrix whose estimated condition number is less than */
/*             1/TOL is considered to be of full rank.  If the user sets */
/*             TOL <= 0, then an implicitly computed, default tolerance, */
/*             defined by  TOLDEF = N*EPS,  is used instead. */
/*             This parameter is not relevant if COND = 'N'. */

/*     Workspace */

/*     IWORK   INTEGER array, dimension (N+r), where r is the number */
/*             of diagonal blocks R_k in R (see description of LMPARM). */
/*             On output, if INFO = 0, the first N entries of this array */
/*             define a permutation matrix P such that J*P = Q*R, where */
/*             J is the final calculated Jacobian, Q is an orthogonal */
/*             matrix (not stored), and R is upper triangular with */
/*             diagonal elements of nonincreasing magnitude (possibly */
/*             for each block column of J). Column j of P is column */
/*             IWORK(j) of the identity matrix. If INFO = 0, the entries */
/*             N+1:N+r of this array contain the ranks of the final */
/*             submatrices S_k (see description of LMPARM). */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, DWORK(2) returns the residual error norm (the */
/*             sum of squares), DWORK(3) returns the number of iterations */
/*             performed, and DWORK(4) returns the final Levenberg */
/*             factor. If INFO = 0, N > 0, and IWARN >= 0, the elements */
/*             DWORK(5) to DWORK(4+M) contain the final matrix-vector */
/*             product Z*Q'*e, and the elements DWORK(5+M) to */
/*             DWORK(4+M+N*NC) contain the (compressed) representation of */
/*             final upper triangular matrices R and S (if IWARN <> 4). */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( 4, M + max( size(J) + */
/*                                        max( DW( FCN|IFLAG = 1 ), */
/*                                             DW( FCN|IFLAG = 2 ), */
/*                                             DW( QRFACT ) + N ), */
/*                                        N*NC + N + */
/*                                        max( M + DW( FCN|IFLAG = 1 ), */
/*                                             N + DW( LMPARM ) ) ) ), */
/*             where size(J) is the size of the Jacobian (provided by FCN */
/*             in IPAR(1), for IFLAG = 3), and DW( f ) is the workspace */
/*             needed by the routine f, where f is FCN, QRFACT, or LMPARM */
/*             (provided by FCN in IPAR(2:5), for IFLAG = 3). */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             < 0:  the user set IFLAG = IWARN in the subroutine FCN; */
/*             = 1:  both actual and predicted relative reductions in */
/*                   the sum of squares are at most FTOL; */
/*             = 2:  relative error between two consecutive iterates is */
/*                   at most XTOL; */
/*             = 3:  conditions for IWARN = 1 and IWARN = 2 both hold; */
/*             = 4:  the cosine of the angle between e and any column of */
/*                   the Jacobian is at most GTOL in absolute value; */
/*             = 5:  the number of iterations has reached ITMAX without */
/*                   satisfying any convergence condition; */
/*             = 6:  FTOL is too small: no further reduction in the sum */
/*                   of squares is possible; */
/*             = 7:  XTOL is too small: no further improvement in the */
/*                   approximate solution x is possible; */
/*             = 8:  GTOL is too small: e is orthogonal to the columns of */
/*                   the Jacobian to machine precision. */
/*             In all these cases, DWORK(1:4) are set as described above. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  user-defined routine FCN returned with INFO <> 0 */
/*                   for IFLAG = 1; */
/*             = 2:  user-defined routine FCN returned with INFO <> 0 */
/*                   for IFLAG = 2; */
/*             = 3:  user-defined routine QRFACT returned with INFO <> 0; */
/*             = 4:  user-defined routine LMPARM returned with INFO <> 0. */

/*     METHOD */

/*     If XINIT = 'R', the initial value for x is set to a vector of */
/*     pseudo-random values uniformly distributed in (-1,1). */

/*     The Levenberg-Marquardt algorithm (described in [1,3]) is used for */
/*     optimizing the variables x. This algorithm needs the Jacobian */
/*     matrix J, which is provided by the subroutine FCN. A trust region */
/*     method is used. The algorithm tries to update x by the formula */

/*         x = x - p, */

/*     using an approximate solution of the system of linear equations */

/*         (J'*J + PAR*D*D)*p = J'*e, */

/*     with e the error function vector, and D a diagonal nonsingular */
/*     matrix, where either PAR = 0 and */

/*         ( norm( D*x ) - DELTA ) <= 0.1*DELTA , */

/*     or PAR > 0 and */

/*         ABS( norm( D*x ) - DELTA ) <= 0.1*DELTA . */

/*     DELTA is the radius of the trust region. If the Gauss-Newton */
/*     direction is not acceptable, then an iterative algorithm obtains */
/*     improved lower and upper bounds for the Levenberg-Marquardt */
/*     parameter PAR. Only a few iterations are generally needed for */
/*     convergence of the algorithm. The trust region radius DELTA */
/*     and the Levenberg factor PAR are updated based on the ratio */
/*     between the actual and predicted reduction in the sum of squares. */

/*     REFERENCES */

/*     [1] More, J.J., Garbow, B.S, and Hillstrom, K.E. */
/*         User's Guide for MINPACK-1. */
/*         Applied Math. Division, Argonne National Laboratory, Argonne, */
/*         Illinois, Report ANL-80-74, 1980. */

/*     [2] Golub, G.H. and van Loan, C.F. */
/*         Matrix Computations. Third Edition. */
/*         M. D. Johns Hopkins University Press, Baltimore, pp. 520-528, */
/*         1996. */

/*     [3] More, J.J. */
/*         The Levenberg-Marquardt algorithm: implementation and theory. */
/*         In Watson, G.A. (Ed.), Numerical Analysis, Lecture Notes in */
/*         Mathematics, vol. 630, Springer-Verlag, Berlin, Heidelberg */
/*         and New York, pp. 105-116, 1978. */

/*     NUMERICAL ASPECTS */

/*     The Levenberg-Marquardt algorithm described in [3] is scaling */
/*     invariant and globally convergent to (maybe local) minima. */
/*     The convergence rate near a local minimum is quadratic, if the */
/*     Jacobian is computed analytically, and linear, if the Jacobian */
/*     is computed numerically. */

/*     FURTHER COMMENTS */

/*     This routine is a more general version of the subroutines LMDER */
/*     and LMDER1 from the MINPACK package [1], which enables to exploit */
/*     the structure of the problem, and optionally use condition */
/*     estimation. Unstructured problems could be solved as well. */

/*     Template SLICOT Library implementations for FCN, QRFACT and */
/*     LMPARM routines are: */
/*     MD03BF, MD03BA, and MD03BB, respectively, for standard problems; */
/*     NF01BF, NF01BS, and NF01BP, respectively, for optimizing the */
/*     parameters of Wiener systems (structured problems). */

/*     CONTRIBUTORS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Feb. 15, 2004. */

/*     KEYWORDS */

/*     Least-squares approximation,  Levenberg-Marquardt algorithm, */
/*     matrix operations, optimization. */

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

/*     Check the scalar input parameters. */

#line 768 "MD03BD.f"
    /* Parameter adjustments */
#line 768 "MD03BD.f"
    --dwork;
#line 768 "MD03BD.f"
    --iwork;
#line 768 "MD03BD.f"
    --diag;
#line 768 "MD03BD.f"
    --x;
#line 768 "MD03BD.f"
    --dpar2;
#line 768 "MD03BD.f"
    --dpar1;
#line 768 "MD03BD.f"
    --ipar;
#line 768 "MD03BD.f"

#line 768 "MD03BD.f"
    /* Function Body */
#line 768 "MD03BD.f"
    init = lsame_(xinit, "R", (ftnlen)1, (ftnlen)1);
#line 769 "MD03BD.f"
    iscal = lsame_(scale, "I", (ftnlen)1, (ftnlen)1);
#line 770 "MD03BD.f"
    sscal = lsame_(scale, "S", (ftnlen)1, (ftnlen)1);
#line 771 "MD03BD.f"
    *info = 0;
#line 772 "MD03BD.f"
    *iwarn = 0;
#line 773 "MD03BD.f"
    if (! (init || lsame_(xinit, "G", (ftnlen)1, (ftnlen)1))) {
#line 774 "MD03BD.f"
	*info = -1;
#line 775 "MD03BD.f"
    } else if (! (iscal || sscal)) {
#line 776 "MD03BD.f"
	*info = -2;
#line 777 "MD03BD.f"
    } else if (! (lsame_(cond, "E", (ftnlen)1, (ftnlen)1) || lsame_(cond, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 778 "MD03BD.f"
	*info = -3;
#line 779 "MD03BD.f"
    } else if (*m < 0) {
#line 780 "MD03BD.f"
	*info = -7;
#line 781 "MD03BD.f"
    } else if (*n < 0 || *n > *m) {
#line 782 "MD03BD.f"
	*info = -8;
#line 783 "MD03BD.f"
    } else if (*itmax < 0) {
#line 784 "MD03BD.f"
	*info = -9;
#line 785 "MD03BD.f"
    } else if (*factor <= 0.) {
#line 786 "MD03BD.f"
	*info = -10;
#line 787 "MD03BD.f"
    } else if (*lipar < 5) {
#line 788 "MD03BD.f"
	*info = -13;
#line 789 "MD03BD.f"
    } else if (*ldpar1 < 0) {
#line 790 "MD03BD.f"
	*info = -15;
#line 791 "MD03BD.f"
    } else if (*ldpar2 < 0) {
#line 792 "MD03BD.f"
	*info = -17;
#line 793 "MD03BD.f"
    } else if (*ldwork < 4) {
#line 794 "MD03BD.f"
	*info = -28;
#line 795 "MD03BD.f"
    } else if (sscal) {
#line 796 "MD03BD.f"
	badscl = FALSE_;

#line 798 "MD03BD.f"
	i__1 = *n;
#line 798 "MD03BD.f"
	for (j = 1; j <= i__1; ++j) {
#line 799 "MD03BD.f"
	    badscl = badscl || diag[j] <= 0.;
#line 800 "MD03BD.f"
/* L10: */
#line 800 "MD03BD.f"
	}

#line 802 "MD03BD.f"
	if (badscl) {
#line 802 "MD03BD.f"
	    *info = -19;
#line 802 "MD03BD.f"
	}
#line 804 "MD03BD.f"
    }

/*     Return if there are illegal arguments. */

#line 808 "MD03BD.f"
    if (*info != 0) {
#line 809 "MD03BD.f"
	i__1 = -(*info);
#line 809 "MD03BD.f"
	xerbla_("MD03BD", &i__1, (ftnlen)6);
#line 810 "MD03BD.f"
	return 0;
#line 811 "MD03BD.f"
    }

/*     Quick return if possible. */

#line 815 "MD03BD.f"
    *nfev = 0;
#line 816 "MD03BD.f"
    *njev = 0;
#line 817 "MD03BD.f"
    if (*n == 0) {
#line 818 "MD03BD.f"
	dwork[1] = 4.;
#line 819 "MD03BD.f"
	dwork[2] = 0.;
#line 820 "MD03BD.f"
	dwork[3] = 0.;
#line 821 "MD03BD.f"
	dwork[4] = 0.;
#line 822 "MD03BD.f"
	return 0;
#line 823 "MD03BD.f"
    }

/*     Call FCN to get the size of the array J, for storing the Jacobian */
/*     matrix, the leading dimension LDJ and the workspace required */
/*     by FCN for IFLAG = 1 and IFLAG = 2, QRFACT and LMPARM. The */
/*     entries DWORK(1:4) should not be modified by the special call of */
/*     FCN below, if XINIT = 'R' and the values in DWORK(1:4) are */
/*     explicitly desired for initialization of the random number */
/*     generator. */

#line 833 "MD03BD.f"
    iflag = 3;
#line 834 "MD03BD.f"
    iw1 = ipar[1];
#line 835 "MD03BD.f"
    iw2 = ipar[2];
#line 836 "MD03BD.f"
    iw3 = ipar[3];
#line 837 "MD03BD.f"
    jw1 = ipar[4];
#line 838 "MD03BD.f"
    jw2 = ipar[5];

#line 840 "MD03BD.f"
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
	    ldpar2, &x[1], &nfevl, &dwork[1], &dwork[1], &ldjsav, &dwork[1], 
	    ldwork, &infol);
#line 842 "MD03BD.f"
    sizej = ipar[1];
#line 843 "MD03BD.f"
    lfcn1 = ipar[2];
#line 844 "MD03BD.f"
    lfcn2 = ipar[3];
#line 845 "MD03BD.f"
    lqrf = ipar[4];
#line 846 "MD03BD.f"
    llmp = ipar[5];
#line 847 "MD03BD.f"
    if (ldjsav > 0) {
#line 848 "MD03BD.f"
	nc = sizej / ldjsav;
#line 849 "MD03BD.f"
    } else {
#line 850 "MD03BD.f"
	nc = sizej;
#line 851 "MD03BD.f"
    }

#line 853 "MD03BD.f"
    ipar[1] = iw1;
#line 854 "MD03BD.f"
    ipar[2] = iw2;
#line 855 "MD03BD.f"
    ipar[3] = iw3;
#line 856 "MD03BD.f"
    ipar[4] = jw1;
#line 857 "MD03BD.f"
    ipar[5] = jw2;

/*     Check the workspace length. */

#line 861 "MD03BD.f"
    e = 1;
#line 862 "MD03BD.f"
    jac = e + *m;
#line 863 "MD03BD.f"
    jw1 = jac + sizej;
#line 864 "MD03BD.f"
    jw2 = jw1 + *n;
#line 865 "MD03BD.f"
    iw1 = jac + *n * nc;
#line 866 "MD03BD.f"
    iw2 = iw1 + *n;
#line 867 "MD03BD.f"
    iw3 = iw2 + *n;
#line 868 "MD03BD.f"
    jwork = iw2 + *m;

/* Computing MAX */
/* Computing MAX */
/* Computing MAX */
#line 870 "MD03BD.f"
    i__5 = max(lfcn1,lfcn2), i__6 = *n + lqrf;
/* Computing MAX */
#line 870 "MD03BD.f"
    i__7 = *m + lfcn1, i__8 = *n + llmp;
#line 870 "MD03BD.f"
    i__3 = sizej + max(i__5,i__6), i__4 = *n * nc + *n + max(i__7,i__8);
#line 870 "MD03BD.f"
    i__1 = 4, i__2 = *m + max(i__3,i__4);
#line 870 "MD03BD.f"
    l = max(i__1,i__2);
#line 872 "MD03BD.f"
    if (*ldwork < l) {
#line 873 "MD03BD.f"
	*info = -28;
#line 874 "MD03BD.f"
	i__1 = -(*info);
#line 874 "MD03BD.f"
	xerbla_("MD03BD", &i__1, (ftnlen)6);
#line 875 "MD03BD.f"
	return 0;
#line 876 "MD03BD.f"
    }

/*     Set default tolerances. EPSMCH is the machine precision. */

#line 880 "MD03BD.f"
    epsmch = dlamch_("Epsilon", (ftnlen)7);
#line 881 "MD03BD.f"
    ftdef = *ftol;
#line 882 "MD03BD.f"
    xtdef = *xtol;
#line 883 "MD03BD.f"
    gtdef = *gtol;
#line 884 "MD03BD.f"
    toldef = *tol;
/* Computing MIN */
#line 885 "MD03BD.f"
    d__1 = min(ftdef,xtdef), d__1 = min(d__1,gtdef);
#line 885 "MD03BD.f"
    if (min(d__1,toldef) <= 0.) {
#line 886 "MD03BD.f"
	if (ftdef < 0.) {
#line 886 "MD03BD.f"
	    ftdef = sqrt(epsmch);
#line 886 "MD03BD.f"
	}
#line 888 "MD03BD.f"
	if (xtdef < 0.) {
#line 888 "MD03BD.f"
	    xtdef = sqrt(epsmch);
#line 888 "MD03BD.f"
	}
#line 890 "MD03BD.f"
	if (gtdef < 0.) {
#line 890 "MD03BD.f"
	    gtdef = epsmch;
#line 890 "MD03BD.f"
	}
#line 892 "MD03BD.f"
	if (toldef <= 0.) {
#line 892 "MD03BD.f"
	    toldef = (doublereal) (*n) * epsmch;
#line 892 "MD03BD.f"
	}
#line 894 "MD03BD.f"
    }
#line 895 "MD03BD.f"
    wrkopt = 1;

/*     Initialization. */

#line 899 "MD03BD.f"
    if (init) {

/*        SEED is the initial state of the random number generator. */
/*        SEED(4) must be odd. */

#line 904 "MD03BD.f"
	seed[0] = (integer) dwork[1] % 4096;
#line 905 "MD03BD.f"
	seed[1] = (integer) dwork[2] % 4096;
#line 906 "MD03BD.f"
	seed[2] = (integer) dwork[3] % 4096;
#line 907 "MD03BD.f"
	seed[3] = (((integer) dwork[4] << 1) + 1) % 4096;
#line 908 "MD03BD.f"
	dlarnv_(&c__2, seed, n, &x[1]);
#line 909 "MD03BD.f"
    }

/*     Initialize Levenberg-Marquardt parameter and iteration counter. */

#line 913 "MD03BD.f"
    par = 0.;
#line 914 "MD03BD.f"
    iter = 1;

/*     Evaluate the function at the starting point */
/*     and calculate its norm. */
/*     Workspace: need:    M + SIZEJ + LFCN1; */
/*                prefer:  larger. */

#line 921 "MD03BD.f"
    iflag = 1;
#line 922 "MD03BD.f"
    i__1 = *ldwork - jw1 + 1;
#line 922 "MD03BD.f"
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
	    ldpar2, &x[1], &nfevl, &dwork[e], &dwork[jac], &ldj, &dwork[jw1], 
	    &i__1, &infol);

#line 926 "MD03BD.f"
    if (infol != 0) {
#line 927 "MD03BD.f"
	*info = 1;
#line 928 "MD03BD.f"
	return 0;
#line 929 "MD03BD.f"
    }
/* Computing MAX */
#line 930 "MD03BD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jw1] + jw1 - 1;
#line 930 "MD03BD.f"
    wrkopt = max(i__1,i__2);
#line 931 "MD03BD.f"
    *nfev = 1;
#line 932 "MD03BD.f"
    fnorm = dnrm2_(m, &dwork[e], &c__1);
#line 933 "MD03BD.f"
    if (iflag < 0 || fnorm == 0.) {
#line 933 "MD03BD.f"
	goto L90;
#line 933 "MD03BD.f"
    }

/*     Beginning of the outer loop. */

#line 938 "MD03BD.f"
L20:

/*        Calculate the Jacobian matrix. */
/*        Workspace: need:    M + SIZEJ + LFCN2; */
/*                   prefer:  larger. */

#line 944 "MD03BD.f"
    ldj = ldjsav;
#line 945 "MD03BD.f"
    iflag = 2;
#line 946 "MD03BD.f"
    i__1 = *ldwork - jw1 + 1;
#line 946 "MD03BD.f"
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
	    ldpar2, &x[1], &nfevl, &dwork[e], &dwork[jac], &ldj, &dwork[jw1], 
	    &i__1, &infol);

#line 950 "MD03BD.f"
    if (infol != 0) {
#line 951 "MD03BD.f"
	*info = 2;
#line 952 "MD03BD.f"
	return 0;
#line 953 "MD03BD.f"
    }
#line 954 "MD03BD.f"
    if (iter == 1) {
/* Computing MAX */
#line 954 "MD03BD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jw1] + jw1 - 1;
#line 954 "MD03BD.f"
	wrkopt = max(i__1,i__2);
#line 954 "MD03BD.f"
    }
#line 956 "MD03BD.f"
    if (nfevl > 0) {
#line 956 "MD03BD.f"
	*nfev += nfevl;
#line 956 "MD03BD.f"
    }
#line 958 "MD03BD.f"
    ++(*njev);
#line 959 "MD03BD.f"
    if (iflag < 0) {
#line 959 "MD03BD.f"
	goto L90;
#line 959 "MD03BD.f"
    }

/*        If requested, call FCN to enable printing of iterates. */

#line 964 "MD03BD.f"
    if (*nprint > 0) {
#line 965 "MD03BD.f"
	iflag = 0;
#line 966 "MD03BD.f"
	if ((iter - 1) % *nprint == 0) {
#line 967 "MD03BD.f"
	    i__1 = *ldwork - jw1 + 1;
#line 967 "MD03BD.f"
	    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1]
		    , ldpar2, &x[1], nfev, &dwork[e], &dwork[jac], &ldj, &
		    dwork[jw1], &i__1, &infol);

#line 971 "MD03BD.f"
	    if (iflag < 0) {
#line 971 "MD03BD.f"
		goto L90;
#line 971 "MD03BD.f"
	    }
#line 973 "MD03BD.f"
	}
#line 974 "MD03BD.f"
    }

/*        Compute the QR factorization of the Jacobian. */
/*        Workspace: need:    M + SIZEJ + N + LQRF; */
/*                   prefer:  larger. */

#line 980 "MD03BD.f"
    i__1 = *ldwork - jw2 + 1;
#line 980 "MD03BD.f"
    (*qrfact)(n, &ipar[1], lipar, &fnorm, &dwork[jac], &ldj, &dwork[e], &
	    dwork[jw1], &gnorm, &iwork[1], &dwork[jw2], &i__1, &infol);
#line 983 "MD03BD.f"
    if (infol != 0) {
#line 984 "MD03BD.f"
	*info = 3;
#line 985 "MD03BD.f"
	return 0;
#line 986 "MD03BD.f"
    }

/*        On the first iteration and if SCALE = 'I', scale according */
/*        to the norms of the columns of the initial Jacobian. */

#line 991 "MD03BD.f"
    if (iter == 1) {
/* Computing MAX */
#line 992 "MD03BD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jw2] + jw2 - 1;
#line 992 "MD03BD.f"
	wrkopt = max(i__1,i__2);
#line 993 "MD03BD.f"
	if (iscal) {

#line 995 "MD03BD.f"
	    i__1 = *n;
#line 995 "MD03BD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 996 "MD03BD.f"
		diag[j] = dwork[jw1 + j - 1];
#line 997 "MD03BD.f"
		if (diag[j] == 0.) {
#line 997 "MD03BD.f"
		    diag[j] = 1.;
#line 997 "MD03BD.f"
		}
#line 999 "MD03BD.f"
/* L30: */
#line 999 "MD03BD.f"
	    }

#line 1001 "MD03BD.f"
	}

/*           On the first iteration, calculate the norm of the scaled */
/*           x and initialize the step bound DELTA. */

#line 1006 "MD03BD.f"
	i__1 = *n;
#line 1006 "MD03BD.f"
	for (j = 1; j <= i__1; ++j) {
#line 1007 "MD03BD.f"
	    dwork[iw1 + j - 1] = diag[j] * x[j];
#line 1008 "MD03BD.f"
/* L40: */
#line 1008 "MD03BD.f"
	}

#line 1010 "MD03BD.f"
	xnorm = dnrm2_(n, &dwork[iw1], &c__1);
#line 1011 "MD03BD.f"
	delta = *factor * xnorm;
#line 1012 "MD03BD.f"
	if (delta == 0.) {
#line 1012 "MD03BD.f"
	    delta = *factor;
#line 1012 "MD03BD.f"
	}
#line 1014 "MD03BD.f"
    } else {

/*           Rescale if necessary. */

#line 1018 "MD03BD.f"
	if (iscal) {

#line 1020 "MD03BD.f"
	    i__1 = *n;
#line 1020 "MD03BD.f"
	    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
#line 1021 "MD03BD.f"
		d__1 = diag[j], d__2 = dwork[jw1 + j - 1];
#line 1021 "MD03BD.f"
		diag[j] = max(d__1,d__2);
#line 1022 "MD03BD.f"
/* L50: */
#line 1022 "MD03BD.f"
	    }

#line 1024 "MD03BD.f"
	}
#line 1025 "MD03BD.f"
    }

/*        Test for convergence of the gradient norm. */

#line 1029 "MD03BD.f"
    if (gnorm <= gtdef) {
#line 1029 "MD03BD.f"
	*iwarn = 4;
#line 1029 "MD03BD.f"
    }
#line 1031 "MD03BD.f"
    if (*iwarn != 0) {
#line 1031 "MD03BD.f"
	goto L90;
#line 1031 "MD03BD.f"
    }

/*        Beginning of the inner loop. */

#line 1036 "MD03BD.f"
L60:

/*           Determine the Levenberg-Marquardt parameter and the */
/*           direction p, and compute -R*P'*p. */
/*           Workspace:  need:    M + N*NC + 2*N + LLMP; */
/*                       prefer:  larger. */

#line 1043 "MD03BD.f"
    i__1 = *ldwork - iw3 + 1;
#line 1043 "MD03BD.f"
    (*lmparm)(cond, n, &ipar[1], lipar, &dwork[jac], &ldj, &iwork[1], &diag[1]
	    , &dwork[e], &delta, &par, &iwork[*n + 1], &dwork[iw1], &dwork[
	    iw2], &toldef, &dwork[iw3], &i__1, &infol, (ftnlen)1);
#line 1047 "MD03BD.f"
    if (infol != 0) {
#line 1048 "MD03BD.f"
	*info = 4;
#line 1049 "MD03BD.f"
	return 0;
#line 1050 "MD03BD.f"
    }
#line 1051 "MD03BD.f"
    if (iter == 1) {
/* Computing MAX */
#line 1051 "MD03BD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[iw3] + iw3 - 1;
#line 1051 "MD03BD.f"
	wrkopt = max(i__1,i__2);
#line 1051 "MD03BD.f"
    }

#line 1054 "MD03BD.f"
    temp1 = dnrm2_(n, &dwork[iw2], &c__1) / fnorm;

/*           Store the direction p and x - p. */

#line 1058 "MD03BD.f"
    i__1 = *n - 1;
#line 1058 "MD03BD.f"
    for (j = 0; j <= i__1; ++j) {
#line 1059 "MD03BD.f"
	dwork[iw2 + j] = diag[j + 1] * dwork[iw1 + j];
#line 1060 "MD03BD.f"
	dwork[iw1 + j] = x[j + 1] - dwork[iw1 + j];
#line 1061 "MD03BD.f"
/* L70: */
#line 1061 "MD03BD.f"
    }

/*           Compute the norm of scaled p and the scaled predicted */
/*           reduction and the scaled directional derivative. */

#line 1066 "MD03BD.f"
    pnorm = dnrm2_(n, &dwork[iw2], &c__1);
#line 1067 "MD03BD.f"
    temp2 = sqrt(par) * pnorm / fnorm;
/* Computing 2nd power */
#line 1068 "MD03BD.f"
    d__1 = temp1;
/* Computing 2nd power */
#line 1068 "MD03BD.f"
    d__2 = temp2;
#line 1068 "MD03BD.f"
    prered = d__1 * d__1 + d__2 * d__2 / .5;
/* Computing 2nd power */
#line 1069 "MD03BD.f"
    d__1 = temp1;
/* Computing 2nd power */
#line 1069 "MD03BD.f"
    d__2 = temp2;
#line 1069 "MD03BD.f"
    dirder = -(d__1 * d__1 + d__2 * d__2);

/*           On the first iteration, adjust the initial step bound. */

#line 1073 "MD03BD.f"
    if (iter == 1) {
#line 1073 "MD03BD.f"
	delta = min(delta,pnorm);
#line 1073 "MD03BD.f"
    }

/*           Evaluate the function at x - p and calculate its norm. */
/*           Workspace:  need:    2*M + N*NC + N + LFCN1; */
/*                       prefer:  larger. */

#line 1080 "MD03BD.f"
    iflag = 1;
#line 1081 "MD03BD.f"
    i__1 = *ldwork - jwork + 1;
#line 1081 "MD03BD.f"
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
	    ldpar2, &dwork[iw1], &nfevl, &dwork[iw2], &dwork[jac], &ldj, &
	    dwork[jwork], &i__1, &infol);
#line 1084 "MD03BD.f"
    if (infol != 0) {
#line 1085 "MD03BD.f"
	*info = 1;
#line 1086 "MD03BD.f"
	return 0;
#line 1087 "MD03BD.f"
    }

#line 1089 "MD03BD.f"
    ++(*nfev);
#line 1090 "MD03BD.f"
    if (iflag < 0) {
#line 1090 "MD03BD.f"
	goto L90;
#line 1090 "MD03BD.f"
    }
#line 1092 "MD03BD.f"
    fnorm1 = dnrm2_(m, &dwork[iw2], &c__1);

/*           Compute the scaled actual reduction. */

#line 1096 "MD03BD.f"
    actred = -1.;
#line 1097 "MD03BD.f"
    if (fnorm1 * .1 < fnorm) {
/* Computing 2nd power */
#line 1097 "MD03BD.f"
	d__1 = fnorm1 / fnorm;
#line 1097 "MD03BD.f"
	actred = 1. - d__1 * d__1;
#line 1097 "MD03BD.f"
    }

/*           Compute the ratio of the actual to the predicted reduction. */

#line 1102 "MD03BD.f"
    ratio = 0.;
#line 1103 "MD03BD.f"
    if (prered != 0.) {
#line 1103 "MD03BD.f"
	ratio = actred / prered;
#line 1103 "MD03BD.f"
    }

/*           Update the step bound. */

#line 1108 "MD03BD.f"
    if (ratio <= .25) {
#line 1109 "MD03BD.f"
	if (actred >= 0.) {
#line 1110 "MD03BD.f"
	    temp = .5;
#line 1111 "MD03BD.f"
	} else {
#line 1112 "MD03BD.f"
	    temp = dirder * .5 / (dirder + actred * .5);
#line 1113 "MD03BD.f"
	}
#line 1114 "MD03BD.f"
	if (fnorm1 * .1 >= fnorm || temp < .1) {
#line 1114 "MD03BD.f"
	    temp = .1;
#line 1114 "MD03BD.f"
	}
/* Computing MIN */
#line 1116 "MD03BD.f"
	d__1 = delta, d__2 = pnorm / .1;
#line 1116 "MD03BD.f"
	delta = temp * min(d__1,d__2);
#line 1117 "MD03BD.f"
	par /= temp;
#line 1118 "MD03BD.f"
    } else {
#line 1119 "MD03BD.f"
	if (par == 0. || ratio >= .75) {
#line 1120 "MD03BD.f"
	    delta = pnorm / .5;
#line 1121 "MD03BD.f"
	    par *= .5;
#line 1122 "MD03BD.f"
	}
#line 1123 "MD03BD.f"
    }

/*           Test for successful iteration. */

#line 1127 "MD03BD.f"
    if (ratio >= 1e-4) {

/*              Successful iteration. Update x, e, and their norms. */

#line 1131 "MD03BD.f"
	i__1 = *n;
#line 1131 "MD03BD.f"
	for (j = 1; j <= i__1; ++j) {
#line 1132 "MD03BD.f"
	    x[j] = dwork[iw1 + j - 1];
#line 1133 "MD03BD.f"
	    dwork[iw1 + j - 1] = diag[j] * x[j];
#line 1134 "MD03BD.f"
/* L80: */
#line 1134 "MD03BD.f"
	}

#line 1136 "MD03BD.f"
	dcopy_(m, &dwork[iw2], &c__1, &dwork[e], &c__1);
#line 1137 "MD03BD.f"
	xnorm = dnrm2_(n, &dwork[iw1], &c__1);
#line 1138 "MD03BD.f"
	fnorm = fnorm1;
#line 1139 "MD03BD.f"
	++iter;
#line 1140 "MD03BD.f"
    }

/*           Tests for convergence. */

#line 1144 "MD03BD.f"
    if (abs(actred) <= ftdef && prered <= ftdef && ratio * .5 <= 1.) {
#line 1144 "MD03BD.f"
	*iwarn = 1;
#line 1144 "MD03BD.f"
    }
#line 1147 "MD03BD.f"
    if (delta <= xtdef * xnorm) {
#line 1147 "MD03BD.f"
	*iwarn = 2;
#line 1147 "MD03BD.f"
    }
#line 1149 "MD03BD.f"
    if (abs(actred) <= ftdef && prered <= ftdef && ratio * .5 <= 1. && *iwarn 
	    == 2) {
#line 1149 "MD03BD.f"
	*iwarn = 3;
#line 1149 "MD03BD.f"
    }
#line 1152 "MD03BD.f"
    if (*iwarn != 0) {
#line 1152 "MD03BD.f"
	goto L90;
#line 1152 "MD03BD.f"
    }

/*           Tests for termination and stringent tolerances. */

#line 1157 "MD03BD.f"
    if (iter >= *itmax) {
#line 1157 "MD03BD.f"
	*iwarn = 5;
#line 1157 "MD03BD.f"
    }
#line 1159 "MD03BD.f"
    if (abs(actred) <= epsmch && prered <= epsmch && ratio * .5 <= 1.) {
#line 1159 "MD03BD.f"
	*iwarn = 6;
#line 1159 "MD03BD.f"
    }
#line 1162 "MD03BD.f"
    if (delta <= epsmch * xnorm) {
#line 1162 "MD03BD.f"
	*iwarn = 7;
#line 1162 "MD03BD.f"
    }
#line 1164 "MD03BD.f"
    if (gnorm <= epsmch) {
#line 1164 "MD03BD.f"
	*iwarn = 8;
#line 1164 "MD03BD.f"
    }
#line 1166 "MD03BD.f"
    if (*iwarn != 0) {
#line 1166 "MD03BD.f"
	goto L90;
#line 1166 "MD03BD.f"
    }

/*           End of the inner loop. Repeat if unsuccessful iteration. */

#line 1171 "MD03BD.f"
    if (ratio < 1e-4) {
#line 1171 "MD03BD.f"
	goto L60;
#line 1171 "MD03BD.f"
    }

/*        End of the outer loop. */

#line 1175 "MD03BD.f"
    goto L20;

#line 1177 "MD03BD.f"
L90:

/*     Termination, either normal or user imposed. */
/*     Note that DWORK(JAC) normally contains the results returned by */
/*     QRFACT and LMPARM (the compressed R and S factors). */

#line 1183 "MD03BD.f"
    if (iflag < 0) {
#line 1183 "MD03BD.f"
	*iwarn = iflag;
#line 1183 "MD03BD.f"
    }
#line 1185 "MD03BD.f"
    if (*nprint > 0) {
#line 1186 "MD03BD.f"
	iflag = 0;
#line 1187 "MD03BD.f"
	i__1 = *ldwork - jwork + 1;
#line 1187 "MD03BD.f"
	(*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[1], ldpar1, &dpar2[1], 
		ldpar2, &x[1], nfev, &dwork[e], &dwork[jac], &ldj, &dwork[
		jwork], &i__1, &infol);
#line 1190 "MD03BD.f"
	if (iflag < 0) {
#line 1190 "MD03BD.f"
	    *iwarn = iflag;
#line 1190 "MD03BD.f"
	}
#line 1192 "MD03BD.f"
    }

#line 1194 "MD03BD.f"
    if (*iwarn >= 0) {
#line 1195 "MD03BD.f"
	for (j = *m + *n * nc; j >= 1; --j) {
#line 1196 "MD03BD.f"
	    dwork[j + 4] = dwork[j];
#line 1197 "MD03BD.f"
/* L100: */
#line 1197 "MD03BD.f"
	}
#line 1198 "MD03BD.f"
    }
#line 1199 "MD03BD.f"
    dwork[1] = (doublereal) wrkopt;
#line 1200 "MD03BD.f"
    dwork[2] = fnorm;
#line 1201 "MD03BD.f"
    dwork[3] = (doublereal) iter;
#line 1202 "MD03BD.f"
    dwork[4] = par;

#line 1204 "MD03BD.f"
    return 0;
/* *** Last line of MD03BD *** */
} /* md03bd_ */

