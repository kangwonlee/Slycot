#line 1 "MD03AD.f"
/* MD03AD.f -- translated by f2c (version 20100827).
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

#line 1 "MD03AD.f"
/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;
static integer c__0 = 0;

/* Subroutine */ int md03ad_(char *xinit, char *alg, char *stor, char *uplo, 
	S_fp fcn, U_fp jpj, integer *m, integer *n, integer *itmax, integer *
	nprint, integer *ipar, integer *lipar, doublereal *dpar1, integer *
	ldpar1, doublereal *dpar2, integer *ldpar2, doublereal *x, integer *
	nfev, integer *njev, doublereal *tol, doublereal *cgtol, doublereal *
	dwork, integer *ldwork, integer *iwarn, integer *info, ftnlen 
	xinit_len, ftnlen alg_len, ftnlen stor_len, ftnlen uplo_len)
{
    /* System generated locals */
    integer dpar1_dim1, dpar1_offset, dpar2_dim1, dpar2_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer e, i__, iw1, iw2, jw1, jw2, jac, ldj, jte;
    static doublereal par;
    static integer ldw, seed[4];
    static logical chol;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical full, init;
    static integer iter, ljtj, lfcn1, lfcn2;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    static integer iflag;
    extern /* Subroutine */ int mb02wd_(char *, U_fp, integer *, integer *, 
	    integer *, doublereal *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    ftnlen), mb02xd_(char *, char *, char *, U_fp, integer *, integer 
	    *, integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, ftnlen, ftnlen, 
	    ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static integer infol, ljtjd, nfevl, ljtji;
    static doublereal gsmin;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal fnorm;
    static integer dwjtj, sizej;
    static doublereal gnorm;
    static logical upper;
    static integer jwork;
    static doublereal fnorm1;
    extern /* Subroutine */ int dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal cgtdef, actred;
    static integer itercg;
    static doublereal epsmch, toldef, bignum;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen), dlarnv_(
	    integer *, integer *, integer *, doublereal *);
    static integer iwarnl;
    static doublereal smlnum, sqreps;
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
/*     algorithm, using either a Cholesky-based or a conjugate gradients */
/*     solver. The user must provide a subroutine FCN which calculates */
/*     the functions and the Jacobian J (possibly by finite differences), */
/*     and another subroutine JPJ, which computes either J'*J + par*I */
/*     (if ALG = 'D'), or (J'*J + par*I)*x (if ALG = 'I'), where par is */
/*     the Levenberg factor, exploiting the possible structure of the */
/*     Jacobian matrix. Template implementations of these routines are */
/*     included in the SLICOT Library. */

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

/*     ALG     CHARACTER*1 */
/*             Specifies the algorithm used for solving the linear */
/*             systems involving a Jacobian matrix J, as follows: */
/*             = 'D' :  a direct algorithm, which computes the Cholesky */
/*                      factor of the matrix J'*J + par*I is used; */
/*             = 'I' :  an iterative Conjugate Gradients algorithm, which */
/*                      only needs the matrix J, is used. */
/*             In both cases, matrix J is stored in a compressed form. */

/*     STOR    CHARACTER*1 */
/*             If ALG = 'D', specifies the storage scheme for the */
/*             symmetric matrix J'*J, as follows: */
/*             = 'F' :  full storage is used; */
/*             = 'P' :  packed storage is used. */
/*             The option STOR = 'F' usually ensures a faster execution. */
/*             This parameter is not relevant if ALG = 'I'. */

/*     UPLO    CHARACTER*1 */
/*             If ALG = 'D', specifies which part of the matrix J'*J */
/*             is stored, as follows: */
/*             = 'U' :  the upper triagular part is stored; */
/*             = 'L' :  the lower triagular part is stored. */
/*             The option UPLO = 'U' usually ensures a faster execution. */
/*             This parameter is not relevant if ALG = 'I'. */

/*     Function Parameters */

/*     FCN     EXTERNAL */
/*             Subroutine which evaluates the functions and the Jacobian. */
/*             FCN must be declared in an external statement in the user */
/*             calling program, and must have the following interface: */

/*             SUBROUTINE FCN( IFLAG, M, N, IPAR, LIPAR, DPAR1, LDPAR1, */
/*            $                DPAR2, LDPAR2, X, NFEVL, E, J, LDJ, JTE, */
/*            $                DWORK, LDWORK, INFO ) */

/*             where */

/*             IFLAG   (input/output) INTEGER */
/*                     On entry, this parameter must contain a value */
/*                     defining the computations to be performed: */
/*                     = 0 :  Optionally, print the current iterate X, */
/*                            function values E, and Jacobian matrix J, */
/*                            or other results defined in terms of these */
/*                            values. See the argument NPRINT of MD03AD. */
/*                            Do not alter E and J. */
/*                     = 1 :  Calculate the functions at X and return */
/*                            this vector in E. Do not alter J. */
/*                     = 2 :  Calculate the Jacobian at X and return */
/*                            this matrix in J. Also return J'*e in JTE */
/*                            and NFEVL (see below). Do not alter E. */
/*                     = 3 :  Do not compute neither the functions nor */
/*                            the Jacobian, but return in LDJ and */
/*                            IPAR/DPAR1,DPAR2 (some of) the integer/real */
/*                            parameters needed. */
/*                     On exit, the value of this parameter should not be */
/*                     changed by FCN unless the user wants to terminate */
/*                     execution of MD03AD, in which case IFLAG must be */
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
/*                     JPJ for ALG = 'D', and JPJ for ALG = 'I', */
/*                     respectively. */

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
/*                     LDJ is essentially used inside the routines FCN */
/*                     and JPJ. */
/*                     LDJ is an input parameter, except for IFLAG = 3 */
/*                     on entry, when it is an output parameter. */
/*                     It is assumed in MD03AD that LDJ is not larger */
/*                     than needed. */

/*             JTE     (output) DOUBLE PRECISION array, dimension (N) */
/*                     If IFLAG = 2, the matrix-vector product J'*e. */

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

/*     JPJ     EXTERNAL */
/*             Subroutine which computes J'*J + par*I, if ALG = 'D', and */
/*             J'*J*x + par*x, if ALG = 'I', where J is the Jacobian as */
/*             described above. */

/*             JPJ must have the following interface: */

/*             SUBROUTINE JPJ( STOR, UPLO, N, IPAR, LIPAR, DPAR, LDPAR, */
/*            $                J, LDJ, JTJ, LDJTJ, DWORK, LDWORK, INFO ) */

/*             if ALG = 'D', and */

/*             SUBROUTINE JPJ( N, IPAR, LIPAR, DPAR, LDPAR, J, LDJ, X, */
/*            $                INCX, DWORK, LDWORK, INFO ) */

/*             if ALG = 'I', where */

/*             STOR    (input) CHARACTER*1 */
/*                     Specifies the storage scheme for the symmetric */
/*                     matrix J'*J, as follows: */
/*                     = 'F' :  full storage is used; */
/*                     = 'P' :  packed storage is used. */

/*             UPLO    (input) CHARACTER*1 */
/*                     Specifies which part of the matrix J'*J is stored, */
/*                     as follows: */
/*                     = 'U' :  the upper triagular part is stored; */
/*                     = 'L' :  the lower triagular part is stored. */

/*             N       (input) INTEGER */
/*                     The number of columns of the matrix J.  N >= 0. */

/*             IPAR    (input) INTEGER array, dimension (LIPAR) */
/*                     The integer parameters describing the structure of */
/*                     the Jacobian matrix. */

/*             LIPAR   (input) INTEGER */
/*                     The length of the array IPAR.  LIPAR >= 0. */

/*             DPAR    (input) DOUBLE PRECISION array, dimension (LDPAR) */
/*                     DPAR(1) must contain an initial estimate of the */
/*                     Levenberg-Marquardt parameter, par.  DPAR(1) >= 0. */

/*             LDPAR   (input) INTEGER */
/*                     The length of the array DPAR.  LDPAR >= 1. */

/*             J       (input) DOUBLE PRECISION array, dimension */
/*                     (LDJ, NC), where NC is the number of columns. */
/*                     The leading NR-by-NC part of this array must */
/*                     contain the (compressed) representation of the */
/*                     Jacobian matrix J, where NR is the number of rows */
/*                     of J (function of IPAR entries). */

/*             LDJ     (input) INTEGER */
/*                     The leading dimension of array J. */
/*                     LDJ >= MAX(1,NR). */

/*             JTJ     (output) DOUBLE PRECISION array, */
/*                              dimension (LDJTJ,N),    if STOR = 'F', */
/*                              dimension (N*(N+1)/2),  if STOR = 'P'. */
/*                     The leading N-by-N (if STOR = 'F'), or N*(N+1)/2 */
/*                     (if STOR = 'P') part of this array contains the */
/*                     upper or lower triangle of the matrix J'*J+par*I, */
/*                     depending on UPLO = 'U', or UPLO = 'L', */
/*                     respectively, stored either as a two-dimensional, */
/*                     or one-dimensional array, depending on STOR. */

/*             LDJTJ   (input) INTEGER */
/*                     The leading dimension of the array JTJ. */
/*                     LDJTJ >= MAX(1,N), if STOR = 'F'. */
/*                     LDJTJ >= 1,        if STOR = 'P'. */

/*             DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*                     The workspace array for subroutine JPJ. */

/*             LDWORK  (input) INTEGER */
/*                     The size of the array DWORK (as large as needed */
/*                     in the subroutine JPJ). */

/*             INFO    INTEGER */
/*                     Error indicator, set to a negative value if an */
/*                     input (scalar) argument is erroneous, and to */
/*                     positive values for other possible errors in the */
/*                     subroutine JPJ. The LAPACK Library routine XERBLA */
/*                     should be used in conjunction with negative INFO */
/*                     values. INFO must be zero if the subroutine */
/*                     finished successfully. */

/*             If ALG = 'I', the parameters in common with those for */
/*             ALG = 'D', have the same meaning, and the additional */
/*             parameters are: */

/*             X       (input/output) DOUBLE PRECISION array, dimension */
/*                     (1+(N-1)*INCX) */
/*                     On entry, this incremented array must contain the */
/*                     vector x. */
/*                     On exit, this incremented array contains the value */
/*                     of the matrix-vector product (J'*J + par)*x. */

/*             INCX    (input) INTEGER */
/*                     The increment for the elements of X.  INCX > 0. */

/*             Parameters marked with "(input)" must not be changed. */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of functions.  M >= 0. */

/*     N       (input) INTEGER */
/*             The number of variables.  M >= N >= 0. */

/*     ITMAX   (input) INTEGER */
/*             The maximum number of iterations.  ITMAX >= 0. */

/*     NPRINT  (input) INTEGER */
/*             This parameter enables controlled printing of iterates if */
/*             it is positive. In this case, FCN is called with IFLAG = 0 */
/*             at the beginning of the first iteration and every NPRINT */
/*             iterations thereafter and immediately prior to return, */
/*             with X, E, and J available for printing. If NPRINT is not */
/*             positive, no special calls of FCN with IFLAG = 0 are made. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters needed, for instance, for */
/*             describing the structure of the Jacobian matrix, which */
/*             are handed over to the routines FCN and JPJ. */
/*             The first five entries of this array are modified */
/*             internally by a call to FCN (with IFLAG = 3), but are */
/*             restored on exit. */

/*     LIPAR   (input) INTEGER */
/*             The length of the array IPAR.  LIPAR >= 5. */

/*     DPAR1   (input/output) DOUBLE PRECISION array, dimension */
/*             (LDPAR1,*) or (LDPAR1) */
/*             A first set of real parameters needed for describing or */
/*             solving the problem. This argument is not used by MD03AD */
/*             routine, but it is passed to the routine FCN. */

/*     LDPAR1  (input) INTEGER */
/*             The leading dimension or the length of the array DPAR1, as */
/*             convenient.  LDPAR1 >= 0.  (LDPAR1 >= 1, if leading */
/*             dimension.) */

/*     DPAR2   (input/output) DOUBLE PRECISION array, dimension */
/*             (LDPAR2,*) or (LDPAR2) */
/*             A second set of real parameters needed for describing or */
/*             solving the problem. This argument is not used by MD03AD */
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
/*             DWORK(1:5) give details on the iterative process. */

/*     NFEV    (output) INTEGER */
/*             The number of calls to FCN with IFLAG = 1. If FCN is */
/*             properly implemented, this includes the function */
/*             evaluations needed for finite difference approximation */
/*             of the Jacobian. */

/*     NJEV    (output) INTEGER */
/*             The number of calls to FCN with IFLAG = 2. */

/*     Tolerances */

/*     TOL     DOUBLE PRECISION */
/*             If TOL >= 0, the tolerance which measures the relative */
/*             error desired in the sum of squares. Termination occurs */
/*             when the actual relative reduction in the sum of squares */
/*             is at most TOL. If the user sets  TOL < 0, then  SQRT(EPS) */
/*             is used instead TOL, where EPS is the machine precision */
/*             (see LAPACK Library routine DLAMCH). */

/*     CGTOL   DOUBLE PRECISION */
/*             If ALG = 'I' and CGTOL > 0, the tolerance which measures */
/*             the relative residual of the solutions computed by the */
/*             conjugate gradients (CG) algorithm. Termination of a */
/*             CG process occurs when the relative residual is at */
/*             most CGTOL. If the user sets  CGTOL <= 0, then  SQRT(EPS) */
/*             is used instead CGTOL. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK, DWORK(2) returns the residual error norm (the */
/*             sum of squares), DWORK(3) returns the number of iterations */
/*             performed, DWORK(4) returns the total number of conjugate */
/*             gradients iterations performed (zero, if ALG = 'D'), and */
/*             DWORK(5) returns the final Levenberg factor. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= max( 5, M + 2*N + size(J) + */
/*                            max( DW( FCN|IFLAG = 1 ) + N, */
/*                                 DW( FCN|IFLAG = 2 ), */
/*                                 DW( sol ) ) ), */
/*             where size(J) is the size of the Jacobian (provided by FCN */
/*             in IPAR(1), for IFLAG = 3), DW( f ) is the workspace */
/*             needed by the routine f, where f is FCN or JPJ (provided */
/*             by FCN in IPAR(2:5), for IFLAG = 3), and DW( sol ) is the */
/*             workspace needed for solving linear systems, */
/*             DW( sol ) = N*N + DW( JPJ ),  if ALG = 'D', STOR = 'F'; */
/*             DW( sol ) = N*(N+1)/2 + DW( JPJ ), */
/*                                           if ALG = 'D', STOR = 'P'; */
/*             DW( sol ) = 3*N + DW( JPJ ),  if ALG = 'I'. */

/*     Warning Indicator */

/*     IWARN   INTEGER */
/*             < 0:  the user set IFLAG = IWARN in the subroutine FCN; */
/*             = 0:  no warning; */
/*             = 1:  if the iterative process did not converge in ITMAX */
/*                   iterations with tolerance TOL; */
/*             = 2:  if ALG = 'I', and in one or more iterations of the */
/*                   Levenberg-Marquardt algorithm, the conjugate */
/*                   gradient algorithm did not finish after 3*N */
/*                   iterations, with the accuracy required in the */
/*                   call; */
/*             = 3:  the cosine of the angle between e and any column of */
/*                   the Jacobian is at most FACTOR*EPS in absolute */
/*                   value, where FACTOR = 100 is defined in a PARAMETER */
/*                   statement; */
/*             = 4:  TOL is too small: no further reduction in the sum */
/*                   of squares is possible. */
/*                   In all these cases, DWORK(1:5) are set as described */
/*                   above. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  user-defined routine FCN returned with INFO <> 0 */
/*                   for IFLAG = 1; */
/*             = 2:  user-defined routine FCN returned with INFO <> 0 */
/*                   for IFLAG = 2; */
/*             = 3:  SLICOT Library routine MB02XD, if ALG = 'D', or */
/*                   SLICOT Library routine MB02WD, if ALG = 'I' (or */
/*                   user-defined routine JPJ), returned with INFO <> 0. */

/*     METHOD */

/*     If XINIT = 'R', the initial value for X is set to a vector of */
/*     pseudo-random values uniformly distributed in [-1,1]. */

/*     The Levenberg-Marquardt algorithm (described in [1]) is used for */
/*     optimizing the parameters. This algorithm needs the Jacobian */
/*     matrix J, which is provided by the subroutine FCN. The algorithm */
/*     tries to update x by the formula */

/*         x = x - p, */

/*     using the solution of the system of linear equations */

/*         (J'*J + PAR*I)*p = J'*e, */

/*     where I is the identity matrix, and e the error function vector. */
/*     The Levenberg factor PAR is decreased after each successfull step */
/*     and increased in the other case. */

/*     If ALG = 'D', a direct method, which evaluates the matrix product */
/*     J'*J + par*I and then factors it using Cholesky algorithm, */
/*     implemented in the SLICOT Libray routine MB02XD, is used for */
/*     solving the linear system above. */

/*     If ALG = 'I', the Conjugate Gradients method, described in [2], */
/*     and implemented in the SLICOT Libray routine MB02WD, is used for */
/*     solving the linear system above. The main advantage of this method */
/*     is that in most cases the solution of the system can be computed */
/*     in less time than the time needed to compute the matrix J'*J */
/*     This is, however, problem dependent. */

/*     REFERENCES */

/*     [1] Kelley, C.T. */
/*         Iterative Methods for Optimization. */
/*         Society for Industrial and Applied Mathematics (SIAM), */
/*         Philadelphia (Pa.), 1999. */

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
/*     According to [1], the convergence rate near a local minimum is */
/*     quadratic, if the Jacobian is computed analytically, and linear, */
/*     if the Jacobian is computed numerically. */

/*     Whether or not the direct algorithm is faster than the iterative */
/*     Conjugate Gradients algorithm for solving the linear systems */
/*     involved depends on several factors, including the conditioning */
/*     of the Jacobian matrix, and the ratio between its dimensions. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Dec. 2001. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     Mar. 2002. */

/*     KEYWORDS */

/*     Conjugate gradients, least-squares approximation, */
/*     Levenberg-Marquardt algorithm, matrix operations, optimization. */

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

#line 597 "MD03AD.f"
    /* Parameter adjustments */
#line 597 "MD03AD.f"
    --ipar;
#line 597 "MD03AD.f"
    dpar1_dim1 = *ldpar1;
#line 597 "MD03AD.f"
    dpar1_offset = 1 + dpar1_dim1;
#line 597 "MD03AD.f"
    dpar1 -= dpar1_offset;
#line 597 "MD03AD.f"
    dpar2_dim1 = *ldpar2;
#line 597 "MD03AD.f"
    dpar2_offset = 1 + dpar2_dim1;
#line 597 "MD03AD.f"
    dpar2 -= dpar2_offset;
#line 597 "MD03AD.f"
    --x;
#line 597 "MD03AD.f"
    --dwork;
#line 597 "MD03AD.f"

#line 597 "MD03AD.f"
    /* Function Body */
#line 597 "MD03AD.f"
    init = lsame_(xinit, "R", (ftnlen)1, (ftnlen)1);
#line 598 "MD03AD.f"
    chol = lsame_(alg, "D", (ftnlen)1, (ftnlen)1);
#line 599 "MD03AD.f"
    full = lsame_(stor, "F", (ftnlen)1, (ftnlen)1);
#line 600 "MD03AD.f"
    upper = lsame_(uplo, "U", (ftnlen)1, (ftnlen)1);

/*     Check the scalar input parameters. */

#line 604 "MD03AD.f"
    *iwarn = 0;
#line 605 "MD03AD.f"
    *info = 0;
#line 606 "MD03AD.f"
    if (! (init || lsame_(xinit, "G", (ftnlen)1, (ftnlen)1))) {
#line 607 "MD03AD.f"
	*info = -1;
#line 608 "MD03AD.f"
    } else if (! (chol || lsame_(alg, "I", (ftnlen)1, (ftnlen)1))) {
#line 609 "MD03AD.f"
	*info = -2;
#line 610 "MD03AD.f"
    } else if (chol && ! (full || lsame_(stor, "P", (ftnlen)1, (ftnlen)1))) {
#line 611 "MD03AD.f"
	*info = -3;
#line 612 "MD03AD.f"
    } else if (chol && ! (upper || lsame_(uplo, "L", (ftnlen)1, (ftnlen)1))) {
#line 613 "MD03AD.f"
	*info = -4;
#line 614 "MD03AD.f"
    } else if (*m < 0) {
#line 615 "MD03AD.f"
	*info = -7;
#line 616 "MD03AD.f"
    } else if (*n < 0 || *n > *m) {
#line 617 "MD03AD.f"
	*info = -8;
#line 618 "MD03AD.f"
    } else if (*itmax < 0) {
#line 619 "MD03AD.f"
	*info = -9;
#line 620 "MD03AD.f"
    } else if (*lipar < 5) {
#line 621 "MD03AD.f"
	*info = -12;
#line 622 "MD03AD.f"
    } else if (*ldpar1 < 0) {
#line 623 "MD03AD.f"
	*info = -14;
#line 624 "MD03AD.f"
    } else if (*ldpar2 < 0) {
#line 625 "MD03AD.f"
	*info = -16;
#line 626 "MD03AD.f"
    } else if (*ldwork < 5) {
#line 627 "MD03AD.f"
	*info = -23;
#line 628 "MD03AD.f"
    }

/*     Return if there are illegal arguments. */

#line 632 "MD03AD.f"
    if (*info != 0) {
#line 633 "MD03AD.f"
	i__1 = -(*info);
#line 633 "MD03AD.f"
	xerbla_("MD03AD", &i__1, (ftnlen)6);
#line 634 "MD03AD.f"
	return 0;
#line 635 "MD03AD.f"
    }

/*     Quick return if possible. */

#line 639 "MD03AD.f"
    *nfev = 0;
#line 640 "MD03AD.f"
    *njev = 0;
#line 641 "MD03AD.f"
    if (min(*n,*itmax) == 0) {
#line 642 "MD03AD.f"
	dwork[1] = 5.;
#line 643 "MD03AD.f"
	dwork[2] = 0.;
#line 644 "MD03AD.f"
	dwork[3] = 0.;
#line 645 "MD03AD.f"
	dwork[4] = 0.;
#line 646 "MD03AD.f"
	dwork[5] = 0.;
#line 647 "MD03AD.f"
	return 0;
#line 648 "MD03AD.f"
    }

/*     Call FCN to get the size of the array J, for storing the Jacobian */
/*     matrix, the leading dimension LDJ and the workspace required */
/*     by FCN for IFLAG = 1 and IFLAG = 2, and JPJ. The entries */
/*     DWORK(1:4) should not be modified by the special call of FCN */
/*     below, if XINIT = 'R' and the values in DWORK(1:4) are explicitly */
/*     desired for initialization of the random number generator. */

#line 657 "MD03AD.f"
    iflag = 3;
#line 658 "MD03AD.f"
    iw1 = ipar[1];
#line 659 "MD03AD.f"
    iw2 = ipar[2];
#line 660 "MD03AD.f"
    jw1 = ipar[3];
#line 661 "MD03AD.f"
    jw2 = ipar[4];
#line 662 "MD03AD.f"
    ljtj = ipar[5];

#line 664 "MD03AD.f"
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[dpar1_offset], ldpar1, &
	    dpar2[dpar2_offset], ldpar2, &x[1], &nfevl, &dwork[1], &dwork[1], 
	    &ldj, &dwork[1], &dwork[1], ldwork, &infol);

#line 668 "MD03AD.f"
    sizej = ipar[1];
#line 669 "MD03AD.f"
    lfcn1 = ipar[2];
#line 670 "MD03AD.f"
    lfcn2 = ipar[3];
#line 671 "MD03AD.f"
    ljtjd = ipar[4];
#line 672 "MD03AD.f"
    ljtji = ipar[5];

#line 674 "MD03AD.f"
    ipar[1] = iw1;
#line 675 "MD03AD.f"
    ipar[2] = iw2;
#line 676 "MD03AD.f"
    ipar[3] = jw1;
#line 677 "MD03AD.f"
    ipar[4] = jw2;
#line 678 "MD03AD.f"
    ipar[5] = ljtj;

/*     Define pointers to the array variables stored in DWORK. */

#line 682 "MD03AD.f"
    jac = 1;
#line 683 "MD03AD.f"
    e = jac + sizej;
#line 684 "MD03AD.f"
    jte = e + *m;
#line 685 "MD03AD.f"
    iw1 = jte + *n;
#line 686 "MD03AD.f"
    iw2 = iw1 + *n;
#line 687 "MD03AD.f"
    jw1 = iw2;
#line 688 "MD03AD.f"
    jw2 = iw2 + *n;

/*     Check the workspace length. */

#line 692 "MD03AD.f"
    jwork = jw1;
#line 693 "MD03AD.f"
    if (chol) {
#line 694 "MD03AD.f"
	if (full) {
#line 695 "MD03AD.f"
	    ldw = *n * *n;
#line 696 "MD03AD.f"
	} else {
#line 697 "MD03AD.f"
	    ldw = *n * (*n + 1) / 2;
#line 698 "MD03AD.f"
	}
#line 699 "MD03AD.f"
	dwjtj = jwork;
#line 700 "MD03AD.f"
	jwork = dwjtj + ldw;
#line 701 "MD03AD.f"
	ljtj = ljtjd;
#line 702 "MD03AD.f"
    } else {
#line 703 "MD03AD.f"
	ldw = *n * 3;
#line 704 "MD03AD.f"
	ljtj = ljtji;
#line 705 "MD03AD.f"
    }
/* Computing MAX */
/* Computing MAX */
#line 706 "MD03AD.f"
    i__3 = lfcn1 + *n, i__3 = max(i__3,lfcn2), i__4 = ldw + ljtj;
#line 706 "MD03AD.f"
    i__1 = 5, i__2 = sizej + *m + (*n << 1) + max(i__3,i__4);
#line 706 "MD03AD.f"
    if (*ldwork < max(i__1,i__2)) {
#line 709 "MD03AD.f"
	*info = -23;
#line 710 "MD03AD.f"
    }
#line 711 "MD03AD.f"
    if (*info != 0) {
#line 712 "MD03AD.f"
	i__1 = -(*info);
#line 712 "MD03AD.f"
	xerbla_("MD03AD", &i__1, (ftnlen)6);
#line 713 "MD03AD.f"
	return 0;
#line 714 "MD03AD.f"
    }

/*     Set default tolerances. SQREPS is the square root of the machine */
/*     precision, and GSMIN is used in the tests of the gradient norm. */

#line 719 "MD03AD.f"
    epsmch = dlamch_("Epsilon", (ftnlen)7);
#line 720 "MD03AD.f"
    sqreps = sqrt(epsmch);
#line 721 "MD03AD.f"
    toldef = *tol;
#line 722 "MD03AD.f"
    if (toldef < 0.) {
#line 722 "MD03AD.f"
	toldef = sqreps;
#line 722 "MD03AD.f"
    }
#line 724 "MD03AD.f"
    cgtdef = *cgtol;
#line 725 "MD03AD.f"
    if (cgtdef <= 0.) {
#line 725 "MD03AD.f"
	cgtdef = sqreps;
#line 725 "MD03AD.f"
    }
#line 727 "MD03AD.f"
    gsmin = epsmch * 100.;
#line 728 "MD03AD.f"
    wrkopt = 5;

#line 730 "MD03AD.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 731 "MD03AD.f"
    bignum = 1. / smlnum;
#line 732 "MD03AD.f"
    dlabad_(&smlnum, &bignum);

/*     Initialization. */

#line 736 "MD03AD.f"
    if (init) {

/*        SEED is the initial state of the random number generator. */
/*        SEED(4) must be odd. */

#line 741 "MD03AD.f"
	seed[0] = (integer) dwork[1] % 4096;
#line 742 "MD03AD.f"
	seed[1] = (integer) dwork[2] % 4096;
#line 743 "MD03AD.f"
	seed[2] = (integer) dwork[3] % 4096;
#line 744 "MD03AD.f"
	seed[3] = (((integer) dwork[4] << 1) + 1) % 4096;
#line 745 "MD03AD.f"
	dlarnv_(&c__2, seed, n, &x[1]);
#line 746 "MD03AD.f"
    }

/*     Evaluate the function at the starting point and calculate */
/*     its norm. */
/*     Workspace: need:    SIZEJ + M + 2*N + LFCN1; */
/*                prefer:  larger. */

#line 753 "MD03AD.f"
    iflag = 1;
#line 754 "MD03AD.f"
    i__1 = *ldwork - jw1 + 1;
#line 754 "MD03AD.f"
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[dpar1_offset], ldpar1, &
	    dpar2[dpar2_offset], ldpar2, &x[1], &nfevl, &dwork[e], &dwork[jac]
	    , &ldj, &dwork[jte], &dwork[jw1], &i__1, &infol);

#line 758 "MD03AD.f"
    if (infol != 0) {
#line 759 "MD03AD.f"
	*info = 1;
#line 760 "MD03AD.f"
	return 0;
#line 761 "MD03AD.f"
    }
/* Computing MAX */
#line 762 "MD03AD.f"
    i__1 = wrkopt, i__2 = (integer) dwork[jw1] + jw1 - 1;
#line 762 "MD03AD.f"
    wrkopt = max(i__1,i__2);
#line 763 "MD03AD.f"
    *nfev = 1;
#line 764 "MD03AD.f"
    fnorm = dnrm2_(m, &dwork[e], &c__1);
#line 765 "MD03AD.f"
    actred = 0.;
#line 766 "MD03AD.f"
    itercg = 0;
#line 767 "MD03AD.f"
    iter = 0;
#line 768 "MD03AD.f"
    iwarnl = 0;
#line 769 "MD03AD.f"
    par = 0.;
#line 770 "MD03AD.f"
    if (iflag < 0 || fnorm == 0.) {
#line 770 "MD03AD.f"
	goto L40;
#line 770 "MD03AD.f"
    }

/*     Set the initial vector for the conjugate gradients algorithm. */

#line 775 "MD03AD.f"
    dwork[iw1] = 0.;
#line 776 "MD03AD.f"
    dcopy_(n, &dwork[iw1], &c__0, &dwork[iw1], &c__1);

/*     WHILE ( nonconvergence and ITER < ITMAX ) DO */

/*     Beginning of the outer loop. */

#line 782 "MD03AD.f"
L10:

/*        Calculate the Jacobian matrix. */
/*        Workspace: need:    SIZEJ + M + 2*N + LFCN2; */
/*                   prefer:  larger. */

#line 788 "MD03AD.f"
    ++iter;
#line 789 "MD03AD.f"
    iflag = 2;
#line 790 "MD03AD.f"
    i__1 = *ldwork - jw1 + 1;
#line 790 "MD03AD.f"
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[dpar1_offset], ldpar1, &
	    dpar2[dpar2_offset], ldpar2, &x[1], &nfevl, &dwork[e], &dwork[jac]
	    , &ldj, &dwork[jte], &dwork[jw1], &i__1, &infol);

#line 794 "MD03AD.f"
    if (infol != 0) {
#line 795 "MD03AD.f"
	*info = 2;
#line 796 "MD03AD.f"
	return 0;
#line 797 "MD03AD.f"
    }

/*        Compute the gradient norm. */

#line 801 "MD03AD.f"
    gnorm = dnrm2_(n, &dwork[jte], &c__1);
#line 802 "MD03AD.f"
    if (nfevl > 0) {
#line 802 "MD03AD.f"
	*nfev += nfevl;
#line 802 "MD03AD.f"
    }
#line 804 "MD03AD.f"
    ++(*njev);
#line 805 "MD03AD.f"
    if (gnorm <= gsmin) {
#line 805 "MD03AD.f"
	*iwarn = 3;
#line 805 "MD03AD.f"
    }
#line 807 "MD03AD.f"
    if (*iwarn != 0) {
#line 807 "MD03AD.f"
	goto L40;
#line 807 "MD03AD.f"
    }
#line 809 "MD03AD.f"
    if (iter == 1) {
/* Computing MAX */
#line 810 "MD03AD.f"
	i__1 = wrkopt, i__2 = (integer) dwork[jw1] + jw1 - 1;
#line 810 "MD03AD.f"
	wrkopt = max(i__1,i__2);
/* Computing MIN */
#line 811 "MD03AD.f"
	d__1 = gnorm, d__2 = sqrt(1e20);
#line 811 "MD03AD.f"
	par = min(d__1,d__2);
#line 812 "MD03AD.f"
    }
#line 813 "MD03AD.f"
    if (iflag < 0) {
#line 813 "MD03AD.f"
	goto L40;
#line 813 "MD03AD.f"
    }

/*        If requested, call FCN to enable printing of iterates. */

#line 818 "MD03AD.f"
    if (*nprint > 0) {
#line 819 "MD03AD.f"
	iflag = 0;
#line 820 "MD03AD.f"
	if ((iter - 1) % *nprint == 0) {
#line 821 "MD03AD.f"
	    i__1 = *ldwork - jw1 + 1;
#line 821 "MD03AD.f"
	    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[dpar1_offset], 
		    ldpar1, &dpar2[dpar2_offset], ldpar2, &x[1], nfev, &dwork[
		    e], &dwork[jac], &ldj, &dwork[jte], &dwork[jw1], &i__1, &
		    infol);

#line 825 "MD03AD.f"
	    if (iflag < 0) {
#line 825 "MD03AD.f"
		goto L40;
#line 825 "MD03AD.f"
	    }
#line 827 "MD03AD.f"
	}
#line 828 "MD03AD.f"
    }

/*        Beginning of the inner loop. */

#line 832 "MD03AD.f"
L20:

/*           Store the Levenberg factor in DWORK(E) (which is no longer */
/*           needed), to pass it to JPJ routine. */

#line 837 "MD03AD.f"
    dwork[e] = par;

/*           Solve (J'*J + PAR*I)*x = J'*e, and store x in DWORK(IW1). */
/*           Additional workspace: */
/*                      N*N + DW(JPJ),          if ALG = 'D', STOR = 'F'; */
/*                      N*( N + 1)/2 + DW(JPJ), if ALG = 'D', STOR = 'P'; */
/*                      3*N + DW(JPJ),          if ALG = 'I'. */

#line 845 "MD03AD.f"
    if (chol) {
#line 846 "MD03AD.f"
	dcopy_(n, &dwork[jte], &c__1, &dwork[iw1], &c__1);
#line 847 "MD03AD.f"
	i__1 = *ldwork - jwork + 1;
#line 847 "MD03AD.f"
	mb02xd_("Function", stor, uplo, (U_fp)jpj, m, n, &c__1, &ipar[1], 
		lipar, &dwork[e], &c__1, &dwork[jac], &ldj, &dwork[iw1], n, &
		dwork[dwjtj], n, &dwork[jwork], &i__1, &infol, (ftnlen)8, (
		ftnlen)1, (ftnlen)1);
#line 851 "MD03AD.f"
    } else {
#line 852 "MD03AD.f"
	i__1 = *n * 3;
#line 852 "MD03AD.f"
	d__1 = *cgtol * gnorm;
#line 852 "MD03AD.f"
	i__2 = *ldwork - jwork + 1;
#line 852 "MD03AD.f"
	mb02wd_("Function", (U_fp)jpj, n, &ipar[1], lipar, &dwork[e], &c__1, &
		i__1, &dwork[jac], &ldj, &dwork[jte], &c__1, &dwork[iw1], &
		c__1, &d__1, &dwork[jwork], &i__2, iwarn, &infol, (ftnlen)8);
#line 856 "MD03AD.f"
	itercg += (integer) dwork[jwork];
/* Computing MAX */
#line 857 "MD03AD.f"
	i__1 = *iwarn << 1;
#line 857 "MD03AD.f"
	iwarnl = max(i__1,iwarnl);
#line 858 "MD03AD.f"
    }

#line 860 "MD03AD.f"
    if (infol != 0) {
#line 861 "MD03AD.f"
	*info = 3;
#line 862 "MD03AD.f"
	return 0;
#line 863 "MD03AD.f"
    }

/*           Compute updated X. */

#line 867 "MD03AD.f"
    i__1 = *n - 1;
#line 867 "MD03AD.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 868 "MD03AD.f"
	dwork[iw2 + i__] = x[i__ + 1] - dwork[iw1 + i__];
#line 869 "MD03AD.f"
/* L30: */
#line 869 "MD03AD.f"
    }

/*           Evaluate the function at x - p and calculate its norm. */
/*           Workspace:  need:    SIZEJ + M + 3*N + LFCN1; */
/*                       prefer:  larger. */

#line 875 "MD03AD.f"
    iflag = 1;
#line 876 "MD03AD.f"
    i__1 = *ldwork - jw2 + 1;
#line 876 "MD03AD.f"
    (*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[dpar1_offset], ldpar1, &
	    dpar2[dpar2_offset], ldpar2, &dwork[iw2], &nfevl, &dwork[e], &
	    dwork[jac], &ldj, &dwork[jte], &dwork[jw2], &i__1, &infol);

#line 880 "MD03AD.f"
    if (infol != 0) {
#line 881 "MD03AD.f"
	*info = 1;
#line 882 "MD03AD.f"
	return 0;
#line 883 "MD03AD.f"
    }

#line 885 "MD03AD.f"
    ++(*nfev);
#line 886 "MD03AD.f"
    if (iflag < 0) {
#line 886 "MD03AD.f"
	goto L40;
#line 886 "MD03AD.f"
    }
#line 888 "MD03AD.f"
    fnorm1 = dnrm2_(m, &dwork[e], &c__1);

/*           Now, check whether this step was successful and update the */
/*           Levenberg factor. */

#line 893 "MD03AD.f"
    if (fnorm < fnorm1) {

/*              Unsuccessful step: increase PAR. */

#line 897 "MD03AD.f"
	actred = 1.;
#line 898 "MD03AD.f"
	if (par > 1e20) {
#line 899 "MD03AD.f"
	    if (par / 4. <= bignum) {
#line 899 "MD03AD.f"
		par *= 4.;
#line 899 "MD03AD.f"
	    }
#line 901 "MD03AD.f"
	} else {
#line 902 "MD03AD.f"
	    par *= 4.;
#line 903 "MD03AD.f"
	}

#line 905 "MD03AD.f"
    } else {

/*              Successful step: update PAR, X, and FNORM. */

/* Computing 2nd power */
#line 909 "MD03AD.f"
	d__1 = fnorm1 / fnorm;
#line 909 "MD03AD.f"
	actred = 1. - d__1 * d__1;
#line 910 "MD03AD.f"
	if ((fnorm - fnorm1) * (fnorm + fnorm1) < ddot_(n, &dwork[iw1], &c__1,
		 &dwork[jte], &c__1) * .125) {
#line 913 "MD03AD.f"
	    if (par > 1e20) {
#line 914 "MD03AD.f"
		if (par / 4. <= bignum) {
#line 914 "MD03AD.f"
		    par *= 4.;
#line 914 "MD03AD.f"
		}
#line 916 "MD03AD.f"
	    } else {
#line 917 "MD03AD.f"
		par *= 4.;
#line 918 "MD03AD.f"
	    }
#line 919 "MD03AD.f"
	} else {
/* Computing MAX */
#line 920 "MD03AD.f"
	    d__1 = par / 4.;
#line 920 "MD03AD.f"
	    par = max(d__1,smlnum);
#line 921 "MD03AD.f"
	}
#line 922 "MD03AD.f"
	dcopy_(n, &dwork[iw2], &c__1, &x[1], &c__1);
#line 923 "MD03AD.f"
	fnorm = fnorm1;
#line 924 "MD03AD.f"
    }

#line 926 "MD03AD.f"
    if (actred <= toldef || iter > *itmax || par > 1e20) {
#line 926 "MD03AD.f"
	goto L40;
#line 926 "MD03AD.f"
    }
#line 929 "MD03AD.f"
    if (actred <= epsmch) {
#line 930 "MD03AD.f"
	*iwarn = 4;
#line 931 "MD03AD.f"
	goto L40;
#line 932 "MD03AD.f"
    }

/*           End of the inner loop. Repeat if unsuccessful iteration. */

#line 936 "MD03AD.f"
    if (fnorm < fnorm1) {
#line 936 "MD03AD.f"
	goto L20;
#line 936 "MD03AD.f"
    }

/*        End of the outer loop. */

#line 941 "MD03AD.f"
    goto L10;

/*     END WHILE 10 */

#line 945 "MD03AD.f"
L40:

/*     Termination, either normal or user imposed. */

#line 949 "MD03AD.f"
    if (actred > toldef) {
#line 949 "MD03AD.f"
	*iwarn = 1;
#line 949 "MD03AD.f"
    }
#line 951 "MD03AD.f"
    if (iwarnl != 0) {
#line 951 "MD03AD.f"
	*iwarn = 2;
#line 951 "MD03AD.f"
    }

#line 954 "MD03AD.f"
    if (iflag < 0) {
#line 954 "MD03AD.f"
	*iwarn = iflag;
#line 954 "MD03AD.f"
    }
#line 956 "MD03AD.f"
    if (*nprint > 0) {
#line 957 "MD03AD.f"
	iflag = 0;
#line 958 "MD03AD.f"
	i__1 = *ldwork - jw1 + 1;
#line 958 "MD03AD.f"
	(*fcn)(&iflag, m, n, &ipar[1], lipar, &dpar1[dpar1_offset], ldpar1, &
		dpar2[dpar2_offset], ldpar2, &x[1], nfev, &dwork[e], &dwork[
		jac], &ldj, &dwork[jte], &dwork[jw1], &i__1, &infol);
#line 961 "MD03AD.f"
	if (iflag < 0) {
#line 961 "MD03AD.f"
	    *iwarn = iflag;
#line 961 "MD03AD.f"
	}
#line 963 "MD03AD.f"
    }

#line 965 "MD03AD.f"
    dwork[1] = (doublereal) wrkopt;
#line 966 "MD03AD.f"
    dwork[2] = fnorm;
#line 967 "MD03AD.f"
    dwork[3] = (doublereal) iter;
#line 968 "MD03AD.f"
    dwork[4] = (doublereal) itercg;
#line 969 "MD03AD.f"
    dwork[5] = par;

#line 971 "MD03AD.f"
    return 0;
/* *** Last line of MD03AD *** */
} /* md03ad_ */

