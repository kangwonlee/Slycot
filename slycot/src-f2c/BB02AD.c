#line 1 "BB02AD.f"
/* BB02AD.f -- translated by f2c (version 20100827).
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

#line 1 "BB02AD.f"
/* Table of constant values */

static doublereal c_b7 = 0.;
static integer c__1 = 1;
static doublereal c_b34 = 1.;
static doublereal c_b38 = -4.;
static doublereal c_b40 = 11.;
static integer c__5 = 5;
static doublereal c_b112 = -1.;
static doublereal c_b118 = 4.877;

/* Subroutine */ int bb02ad_(char *def, integer *nr, doublereal *dpar, 
	integer *ipar, logical *bpar, char *chpar, logical *vec, integer *n, 
	integer *m, integer *p, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *q, integer *
	ldq, doublereal *r__, integer *ldr, doublereal *s, integer *lds, 
	doublereal *x, integer *ldx, doublereal *dwork, integer *ldwork, 
	integer *info, ftnlen def_len, ftnlen chpar_len)
{
    /* Initialized data */

    static integer nex[4] = { 13,5,0,1 };
    static integer ndef[52]	/* was [4][13] */ = { 2,2,0,100,2,2,0,0,2,2,0,
	    0,3,3,0,0,4,4,0,0,4,0,0,0,4,0,0,0,5,0,0,0,6,0,0,0,9,0,0,0,11,0,0,
	    0,13,0,0,0,26 };
    static integer mdef[26]	/* was [2][13] */ = { 1,1,2,2,1,1,2,3,2,1,2,0,
	    4,0,2,0,2,0,3,0,2,0,2,0,6 };
    static integer pdef[26]	/* was [2][13] */ = { 1,2,2,2,2,2,3,3,4,1,4,0,
	    4,0,5,0,2,0,2,0,4,0,4,0,12 };
    static struct {
	char e_1[510];
	char fill_2[255];
	char e_3[765];
	char fill_4[510];
	char e_5[510];
	char fill_6[510];
	char e_7[510];
	char fill_8[510];
	char e_9[510];
	char fill_10[510];
	char e_11[255];
	char fill_12[765];
	char e_13[255];
	char fill_14[765];
	char e_15[255];
	char fill_16[765];
	char e_17[255];
	char fill_18[765];
	char e_19[255];
	char fill_20[765];
	char e_21[255];
	char fill_22[765];
	char e_23[255];
	char fill_24[765];
	char e_25[255];
	char fill_26[765];
	} equiv_26 = { "Van Dooren 1981, Ex. II: singular R matrix          "
		"                                                            "
		"                                                            "
		"                                                            "
		"                       Laub 1979, Ex. 2: uncontrollable-unob"
		"servable data                                               "
		"                                                            "
		"                                                            "
		"                                      ", {0}, "Pappas et al."
		" 1980, Ex. 3                                                "
		"                                                            "
		"                                                            "
		"                                                            "
		"  Ionescu/Weiss 1992 : singular R matrix, nonzero S matrix  "
		"                                                            "
		"                                                            "
		"                                                            "
		"                 Laub 1979, Ex. 3: increasingly ill-conditio"
		"ned R-matrix                                                "
		"                                                            "
		"                                                            "
		"                                ", {0}, "Jonckheere 1981: (A"
		",B) controllable, no solution X <= 0                        "
		"                                                            "
		"                                                            "
		"                                                        incr"
		"easingly bad scaled system as eps -> oo                     "
		"                                                            "
		"                                                            "
		"                                                            "
		"           ", {0}, "Sun 1998: R singular, Q non-definite    "
		"                                                            "
		"                                                            "
		"                                                            "
		"                                   Petkov et. al. 1989 : inc"
		"reasingly bad scaling as eps -> oo                          "
		"                                                            "
		"                                                            "
		"                                                  ", {0}, 
		"Ackerson/Fu 1970 : satellite control problem               "
		"                                                            "
		"                                                            "
		"                                                            "
		"                Pappas et al. 1980: process control of paper"
		" machine                                                    "
		"                                                            "
		"                                                            "
		"                               ", {0}, "Litkouhi 1983 : syst"
		"em with slow and fast modes                                 "
		"                                                            "
		"                                                            "
		"                                                       ", {0},
		 "Lu/Lin 1993, Ex. 4.3                                      "
		"                                                            "
		"                                                            "
		"                                                            "
		"                 ", {0}, "Gajic/Shen 1993, Section 2.7.4: ch"
		"emical plant                                                "
		"                                                            "
		"                                                            "
		"                                         ", {0}, "Davison/Wa"
		"ng 1974: nonzero S matrix                                   "
		"                                                            "
		"                                                            "
		"                                                            "
		"     ", {0}, "Patnaik et al. 1980: tubular ammonia reactor  "
		"                                                            "
		"                                                            "
		"                                                            "
		"                             ", {0}, "Sima 1996, Sec. 1.2.2:"
		" paper machine model error integrators                      "
		"                                                            "
		"                                                            "
		"                                                     ", {0}, 
		"Sima 1996, Ex. 2.6: paper machine model with with disturban"
		"ces                                                         "
		"                                                            "
		"                                                            "
		"                ", {0}, "Power plant model, Katayama et al.,"
		" 1985                                                       "
		"                                                            "
		"                                                            "
		"                                        " };

#define notes ((char *)&equiv_26)


    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, s_dim1, 
	    s_offset, x_dim1, x_offset, i__1, i__2;
    doublereal d__1, d__2;
    icilist ici__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    double sqrt(doublereal);
    integer s_wsfi(icilist *), do_fio(integer *, char *, ftnlen), e_wsfi(void)
	    , f_open(olist *), s_rsle(cilist *), do_lio(integer *, integer *, 
	    char *, ftnlen), e_rsle(void), f_clos(cllist *);

    /* Local variables */
    static integer i__, j, ios;
    static doublereal beta, temp;
    extern /* Subroutine */ int dspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen), ma02dd_(char *, 
	    char *, integer *, doublereal *, integer *, doublereal *, ftnlen, 
	    ftnlen), ma02ed_(char *, integer *, doublereal *, integer *, 
	    ftnlen);
    static doublereal alpha;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static char ident[4];
    static integer qdimm, rdimm;
    extern /* Subroutine */ int drscl_(integer *, doublereal *, doublereal *, 
	    integer *), dcopy_(integer *, doublereal *, integer *, doublereal 
	    *, integer *), dspmv_(char *, integer *, doublereal *, doublereal 
	    *, doublereal *, integer *, doublereal *, doublereal *, integer *,
	     ftnlen), dsymm_(char *, char *, integer *, integer *, doublereal 
	    *, doublereal *, integer *, doublereal *, integer *, doublereal *,
	     doublereal *, integer *, ftnlen, ftnlen), dsyrk_(char *, char *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer isymm, msymm, nsymm, psymm;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpptrf_(char *, integer *, 
	    doublereal *, integer *, ftnlen), dpptri_(char *, integer *, 
	    doublereal *, integer *, ftnlen);

    /* Fortran I/O blocks */
    static cilist io___12 = { 1, 1, 1, 0, 0 };
    static cilist io___14 = { 1, 1, 1, 0, 0 };
    static cilist io___15 = { 1, 1, 1, 0, 0 };
    static cilist io___16 = { 1, 1, 1, 0, 0 };
    static cilist io___17 = { 1, 1, 1, 0, 0 };
    static cilist io___18 = { 1, 1, 1, 0, 0 };
    static cilist io___19 = { 1, 1, 1, 0, 0 };



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

/*     To generate the benchmark examples for the numerical solution of */
/*     discrete-time algebraic Riccati equations (DAREs) of the form */

/*            T                T               T    -1  T       T */
/*     0  =  A X A  -  X  -  (A X B + S) (R + B X B)  (B X A + S )  +  Q */

/*     as presented in [1]. Here, A,Q,X are real N-by-N matrices, B,S are */
/*     N-by-M, and R is M-by-M. The matrices Q and R are symmetric and Q */
/*     may be given in factored form */

/*                   T */
/*     (I)    Q  =  C Q0 C . */

/*     Here, C is P-by-N and Q0 is P-by-P. If R is nonsingular and S = 0, */
/*     the DARE can be rewritten equivalently as */

/*                  T             -1 */
/*     0  =  X  -  A X (I_n + G X)  A  -  Q, */

/*     where I_n is the N-by-N identity matrix and */

/*                   -1  T */
/*     (II)   G = B R   B . */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DEF     CHARACTER */
/*             This parameter specifies if the default parameters are */
/*             to be used or not. */
/*             = 'N' or 'n' : The parameters given in the input vectors */
/*                            xPAR (x = 'D', 'I', 'B', 'CH') are used. */
/*             = 'D' or 'd' : The default parameters for the example */
/*                            are used. */
/*             This parameter is not meaningful if NR(1) = 1. */

/*     Input/Output Parameters */

/*     NR      (input) INTEGER array, dimension (2) */
/*             This array determines the example for which DAREX returns */
/*             data. NR(1) is the group of examples. */
/*             NR(1) = 1 : parameter-free problems of fixed size. */
/*             NR(1) = 2 : parameter-dependent problems of fixed size. */
/*             NR(1) = 3 : parameter-free problems of scalable size. */
/*             NR(1) = 4 : parameter-dependent problems of scalable size. */
/*             NR(2) is the number of the example in group NR(1). */
/*             Let NEXi be the number of examples in group i. Currently, */
/*             NEX1 = 13, NEX2 = 5, NEX3 = 0, NEX4 = 1. */
/*             1 <= NR(1) <= 4; */
/*             0 <= NR(2) <= NEXi, where i = NR(1). */

/*     DPAR    (input/output) DOUBLE PRECISION array, dimension (4) */
/*             Double precision parameter vector. For explanation of the */
/*             parameters see [1]. */
/*             DPAR(1) defines the parameter 'epsilon' for */
/*             examples NR = 2.2,2.3,2.4, the parameter 'tau' */
/*             for NR = 2.5, and the 1-by-1 matrix R for NR = 2.1,4.1. */
/*             For Example 2.5, DPAR(2) - DPAR(4) define in */
/*             consecutive order 'D', 'K', and 'r'. */
/*             NOTE that DPAR is overwritten with default values */
/*             if DEF = 'D' or 'd'. */

/*     IPAR    (input/output) INTEGER array, dimension (3) */
/*             On input, IPAR(1) determines the actual state dimension, */
/*             i.e., the order of the matrix A as follows: */
/*             NR(1) = 1, NR(1) = 2   : IPAR(1) is ignored. */
/*             NR = NR(1).NR(2) = 4.1 : IPAR(1) determines the order of */
/*                                      the output matrix A. */
/*             NOTE that IPAR(1) is overwritten for Examples 1.1-2.3. For */
/*             the other examples, IPAR(1) is overwritten if the default */
/*             parameters are to be used. */
/*             On output, IPAR(1) contains the order of the matrix A. */

/*             On input, IPAR(2) is the number of colums in the matrix B */
/*             and the order of the matrix R (in control problems, the */
/*             number of inputs of the system). Currently, IPAR(2) is */
/*             fixed for all examples and thus is not referenced on */
/*             input. */
/*             On output, IPAR(2) is the number of columns of the */
/*             matrix B from (I). */

/*             On input, IPAR(3) is the number of rows in the matrix C */
/*             (in control problems, the number of outputs of the */
/*             system). Currently, IPAR(3) is fixed for all examples */
/*             and thus is not referenced on input. */
/*             On output, IPAR(3) is the number of rows of the matrix C */
/*             from (I). */

/*             NOTE that IPAR(2) and IPAR(3) are overwritten and */
/*             IPAR(2) <= IPAR(1) and IPAR(3) <= IPAR(1) for all */
/*             examples. */

/*     BPAR    (input) LOGICAL array, dimension (7) */
/*             This array defines the form of the output of the examples */
/*             and the storage mode of the matrices Q, G or R. */
/*             BPAR(1) = .TRUE.  : Q is returned. */
/*             BPAR(1) = .FALSE. : Q is returned in factored form, i.e., */
/*                                 Q0 and C from (I) are returned. */
/*             BPAR(2) = .TRUE.  : The matrix returned in array Q (i.e., */
/*                                 Q if BPAR(1) = .TRUE. and Q0 if */
/*                                 BPAR(1) = .FALSE.) is stored as full */
/*                                 matrix. */
/*             BPAR(2) = .FALSE. : The matrix returned in array Q is */
/*                                 provided in packed storage mode. */
/*             BPAR(3) = .TRUE.  : If BPAR(2) = .FALSE., the matrix */
/*                                 returned in array Q is stored in upper */
/*                                 packed mode, i.e., the upper triangle */
/*                                 of a symmetric n-by-n matrix is stored */
/*                                 by columns, e.g., the matrix entry */
/*                                 Q(i,j) is stored in the array entry */
/*                                 Q(i+j*(j-1)/2) for i <= j. */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(3) = .FALSE. : If BPAR(2) = .FALSE., the matrix */
/*                                 returned in array Q is stored in lower */
/*                                 packed mode, i.e., the lower triangle */
/*                                 of a symmetric n-by-n matrix is stored */
/*                                 by columns, e.g., the matrix entry */
/*                                 Q(i,j) is stored in the array entry */
/*                                 Q(i+(2*n-j)*(j-1)/2) for j <= i. */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(4) = .TRUE.  : The product G in (II) is returned. */
/*             BPAR(4) = .FALSE. : G is returned in factored form, i.e., */
/*                                 B and R from (II) are returned. */
/*             BPAR(5) = .TRUE.  : The matrix returned in array R (i.e., */
/*                                 G if BPAR(4) = .TRUE. and R if */
/*                                 BPAR(4) = .FALSE.) is stored as full */
/*                                 matrix. */
/*             BPAR(5) = .FALSE. : The matrix returned in array R is */
/*                                 provided in packed storage mode. */
/*             BPAR(6) = .TRUE.  : If BPAR(5) = .FALSE., the matrix */
/*                                 returned in array R is stored in upper */
/*                                 packed mode (see above). */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(6) = .FALSE. : If BPAR(5) = .FALSE., the matrix */
/*                                 returned in array R is stored in lower */
/*                                 packed mode (see above). */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(7) = .TRUE.  : The coefficient matrix S of the DARE */
/*                                 is returned in array S. */
/*             BPAR(7) = .FALSE. : The coefficient matrix S of the DARE */
/*                                 is not returned. */
/*             NOTE that there are no default values for BPAR.  If all */
/*             entries are declared to be .TRUE., then matrices Q, G or R */
/*             are returned in conventional storage mode, i.e., as */
/*             N-by-N or M-by-M arrays where the array element Z(I,J) */
/*             contains the matrix entry Z_{i,j}. */

/*     CHPAR   (output) CHARACTER*255 */
/*             On output, this string contains short information about */
/*             the chosen example. */

/*     VEC     (output) LOGICAL array, dimension (10) */
/*             Flag vector which displays the availability of the output */
/*             data: */
/*             VEC(j), j=1,2,3, refer to N, M, and P, respectively, and */
/*             are always .TRUE. */
/*             VEC(4) refers to A and is always .TRUE. */
/*             VEC(5) is .TRUE. if BPAR(4) = .FALSE., i.e., the factors B */
/*             and R from (II) are returned. */
/*             VEC(6) is .TRUE. if BPAR(1) = .FALSE., i.e., the factors C */
/*             and Q0 from (I) are returned. */
/*             VEC(7) refers to Q and is always .TRUE. */
/*             VEC(8) refers to R and is always .TRUE. */
/*             VEC(9) is .TRUE. if BPAR(7) = .TRUE., i.e., the matrix S */
/*             is returned. */
/*             VEC(10) refers to X and is .TRUE. if the exact solution */
/*             matrix is available. */
/*             NOTE that VEC(i) = .FALSE. for i = 1 to 10 if on exit */
/*             INFO .NE. 0. */

/*     N       (output) INTEGER */
/*             The order of the matrices A, X, G if BPAR(4) = .TRUE., and */
/*             Q if BPAR(1) = .TRUE. */

/*     M       (output) INTEGER */
/*             The number of columns in the matrix B (or the dimension of */
/*             the control input space of the underlying dynamical */
/*             system). */

/*     P       (output) INTEGER */
/*             The number of rows in the matrix C (or the dimension of */
/*             the output space of the underlying dynamical system). */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             coefficient matrix A of the DARE. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             If (BPAR(4) = .FALSE.), then the leading N-by-M part */
/*             of this array contains the coefficient matrix B of */
/*             the DARE.  Otherwise, B is used as workspace. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= N. */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             If (BPAR(1) = .FALSE.), then the leading P-by-N part */
/*             of this array contains the matrix C of the factored */
/*             form (I) of Q.  Otherwise, C is used as workspace. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= P. */

/*     Q       (output) DOUBLE PRECISION array, dimension (NQ) */
/*             If (BPAR(1) = .TRUE.) and (BPAR(2) = .TRUE.), then */
/*             NQ = LDQ*N. */
/*             IF (BPAR(1) = .TRUE.) and (BPAR(2) = .FALSE.), then */
/*             NQ = N*(N+1)/2. */
/*             If (BPAR(1) = .FALSE.) and (BPAR(2) = .TRUE.), then */
/*             NQ = LDQ*P. */
/*             IF (BPAR(1) = .FALSE.) and (BPAR(2) = .FALSE.), then */
/*             NQ = P*(P+1)/2. */
/*             The symmetric matrix contained in array Q is stored */
/*             according to BPAR(2) and BPAR(3). */

/*     LDQ     INTEGER */
/*             If conventional storage mode is used for Q, i.e., */
/*             BPAR(2) = .TRUE., then Q is stored like a 2-dimensional */
/*             array with leading dimension LDQ. If packed symmetric */
/*             storage mode is used, then LDQ is irrelevant. */
/*             LDQ >= N if BPAR(1) = .TRUE.; */
/*             LDQ >= P if BPAR(1) = .FALSE.. */

/*     R       (output) DOUBLE PRECISION array, dimension (MR) */
/*             If (BPAR(4) = .TRUE.) and (BPAR(5) = .TRUE.), then */
/*             MR = LDR*N. */
/*             IF (BPAR(4) = .TRUE.) and (BPAR(5) = .FALSE.), then */
/*             MR = N*(N+1)/2. */
/*             If (BPAR(4) = .FALSE.) and (BPAR(5) = .TRUE.), then */
/*             MR = LDR*M. */
/*             IF (BPAR(4) = .FALSE.) and (BPAR(5) = .FALSE.), then */
/*             MR = M*(M+1)/2. */
/*             The symmetric matrix contained in array R is stored */
/*             according to BPAR(5) and BPAR(6). */

/*     LDR     INTEGER */
/*             If conventional storage mode is used for R, i.e., */
/*             BPAR(5) = .TRUE., then R is stored like a 2-dimensional */
/*             array with leading dimension LDR. If packed symmetric */
/*             storage mode is used, then LDR is irrelevant. */
/*             LDR >= N  if BPAR(4) =  .TRUE.; */
/*             LDR >= M  if BPAR(4) = .FALSE.. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,M) */
/*             If (BPAR(7) = .TRUE.), then the leading N-by-M part of */
/*             this array contains the coefficient matrix S of the DARE. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= 1, and */
/*             LDS >= N if BPAR(7) = .TRUE.. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,NX) */
/*             If an exact solution is available (NR = 1.1,1.3,1.4,2.1, */
/*             2.3,2.4,2.5,4.1), then NX = N and the leading N-by-N part */
/*             of this array contains the solution matrix X. */
/*             Otherwise, X is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= 1, and */
/*             LDX >= N if an exact solution is available. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= N*N. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0 : successful exit; */
/*             < 0 : if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1 : data file could not be opened or had wrong format; */
/*             = 2 : division by zero; */
/*             = 3 : G can not be computed as in (II) due to a singular R */
/*                   matrix. This error can only occur if */
/*                   BPAR(4) = .TRUE.. */

/*     REFERENCES */

/*     [1] Abels, J. and Benner, P. */
/*         DAREX - A Collection of Benchmark Examples for Discrete-Time */
/*         Algebraic Riccati Equations (Version 2.0). */
/*         SLICOT Working Note 1999-16, November 1999. Available from */
/*         http://www.win.tue.nl/niconet/NIC2/reports.html. */

/*     This is an updated and extended version of */

/*     [2] Benner, P., Laub, A.J., and Mehrmann, V. */
/*         A Collection of Benchmark Examples for the Numerical Solution */
/*         of Algebraic Riccati Equations II: Discrete-Time Case. */
/*         Technical Report SPC 95_23, Fak. f. Mathematik, */
/*         TU Chemnitz-Zwickau (Germany), December 1995. */

/*     FURTHER COMMENTS */

/*     Some benchmark examples read data from the data files provided */
/*     with the collection. */

/*     CONTRIBUTOR */

/*     Peter Benner (Universitaet Bremen), November 25, 1999. */

/*     For questions concerning the collection or for the submission of */
/*     test examples, please send e-mail to benner@math.uni-bremen.de. */

/*     REVISIONS */

/*     1999, December 23 (V. Sima). */

/*     KEYWORDS */

/*     Discrete-time algebraic Riccati equation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     . # of examples available , # of examples with fixed size. . */

/*     .. Scalar Arguments .. */

/*     .. Array Arguments .. */

/*     .. Local Scalars .. */

/*     ..Local Arrays .. */

/*     .. External Functions .. */
/*     . LAPACK . */

/*     .. External Subroutines .. */
/*     . BLAS . */
/*     . LAPACK . */
/*     . SLICOT . */

/*     .. Intrinsic Functions .. */

/*     .. Data Statements .. */
/*     . default values for dimensions . */
#line 396 "BB02AD.f"
    /* Parameter adjustments */
#line 396 "BB02AD.f"
    --nr;
#line 396 "BB02AD.f"
    --dpar;
#line 396 "BB02AD.f"
    --ipar;
#line 396 "BB02AD.f"
    --bpar;
#line 396 "BB02AD.f"
    --vec;
#line 396 "BB02AD.f"
    a_dim1 = *lda;
#line 396 "BB02AD.f"
    a_offset = 1 + a_dim1;
#line 396 "BB02AD.f"
    a -= a_offset;
#line 396 "BB02AD.f"
    b_dim1 = *ldb;
#line 396 "BB02AD.f"
    b_offset = 1 + b_dim1;
#line 396 "BB02AD.f"
    b -= b_offset;
#line 396 "BB02AD.f"
    c_dim1 = *ldc;
#line 396 "BB02AD.f"
    c_offset = 1 + c_dim1;
#line 396 "BB02AD.f"
    c__ -= c_offset;
#line 396 "BB02AD.f"
    --q;
#line 396 "BB02AD.f"
    --r__;
#line 396 "BB02AD.f"
    s_dim1 = *lds;
#line 396 "BB02AD.f"
    s_offset = 1 + s_dim1;
#line 396 "BB02AD.f"
    s -= s_offset;
#line 396 "BB02AD.f"
    x_dim1 = *ldx;
#line 396 "BB02AD.f"
    x_offset = 1 + x_dim1;
#line 396 "BB02AD.f"
    x -= x_offset;
#line 396 "BB02AD.f"
    --dwork;
#line 396 "BB02AD.f"

#line 396 "BB02AD.f"
    /* Function Body */
/*     . comments on examples . */

/*     .. Executable Statements .. */

#line 431 "BB02AD.f"
    *info = 0;
#line 432 "BB02AD.f"
    for (i__ = 1; i__ <= 10; ++i__) {
#line 433 "BB02AD.f"
	vec[i__] = FALSE_;
#line 434 "BB02AD.f"
/* L1: */
#line 434 "BB02AD.f"
    }

#line 436 "BB02AD.f"
    if (nr[1] >= 3) {
#line 437 "BB02AD.f"
	if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 437 "BB02AD.f"
	    ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
#line 437 "BB02AD.f"
	}
#line 438 "BB02AD.f"
	ipar[2] = 1;
#line 439 "BB02AD.f"
	ipar[3] = ipar[1];
#line 440 "BB02AD.f"
    } else {
#line 441 "BB02AD.f"
	ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
#line 442 "BB02AD.f"
	ipar[2] = mdef[nr[1] + (nr[2] << 1) - 3];
#line 443 "BB02AD.f"
	ipar[3] = pdef[nr[1] + (nr[2] << 1) - 3];
#line 444 "BB02AD.f"
    }

#line 446 "BB02AD.f"
    if (nr[1] >= 2 && ! (lsame_(def, "D", (ftnlen)1, (ftnlen)1) || lsame_(def,
	     "N", (ftnlen)1, (ftnlen)1))) {
#line 448 "BB02AD.f"
	*info = -1;
#line 449 "BB02AD.f"
    } else if (nr[1] < 1 || nr[1] > 4 || nr[2] < 0 || nr[2] > nex[nr[1] - 1]) 
	    {
#line 451 "BB02AD.f"
	*info = -2;
#line 452 "BB02AD.f"
    } else if (ipar[1] < 1) {
#line 453 "BB02AD.f"
	*info = -4;
#line 454 "BB02AD.f"
    } else if (ipar[1] > *lda) {
#line 455 "BB02AD.f"
	*info = -12;
#line 456 "BB02AD.f"
    } else if (ipar[1] > *ldb) {
#line 457 "BB02AD.f"
	*info = -14;
#line 458 "BB02AD.f"
    } else if (ipar[3] > *ldc) {
#line 459 "BB02AD.f"
	*info = -16;
#line 460 "BB02AD.f"
    } else if (bpar[2] && (! bpar[1] && ipar[3] > *ldq || bpar[1] && ipar[1] 
	    > *ldq)) {
#line 463 "BB02AD.f"
	*info = -18;
#line 464 "BB02AD.f"
    } else if (bpar[5] && (bpar[4] && ipar[1] > *ldr || ! bpar[4] && ipar[2] 
	    > *ldr)) {
#line 466 "BB02AD.f"
	*info = -20;
#line 467 "BB02AD.f"
    } else if (*lds < 1 || bpar[7] && ipar[1] > *lds) {
#line 468 "BB02AD.f"
	*info = -22;
#line 469 "BB02AD.f"
    } else if (*ldx < 1) {
#line 470 "BB02AD.f"
	*info = -24;
#line 471 "BB02AD.f"
    } else if (nr[1] == 1 && (nr[2] == 1 || nr[2] == 3 || nr[2] == 4) || nr[1]
	     == 2 && (nr[2] == 1 || nr[2] >= 3) || nr[1] == 4) {
/*    .. solution X available .. */
#line 476 "BB02AD.f"
	if (ipar[1] > *ldx) {
#line 477 "BB02AD.f"
	    *info = -24;
#line 478 "BB02AD.f"
	} else {
#line 479 "BB02AD.f"
	    dlaset_("A", &ipar[1], &ipar[1], &c_b7, &c_b7, &x[x_offset], ldx, 
		    (ftnlen)1);
#line 480 "BB02AD.f"
	}
#line 481 "BB02AD.f"
    } else if (*ldwork < *n * *n) {
#line 482 "BB02AD.f"
	*info = -26;
#line 483 "BB02AD.f"
    }
#line 484 "BB02AD.f"
    if (*info != 0) {
#line 485 "BB02AD.f"
	i__1 = -(*info);
#line 485 "BB02AD.f"
	xerbla_("BB02AD", &i__1, (ftnlen)6);
#line 486 "BB02AD.f"
	return 0;
#line 487 "BB02AD.f"
    }

#line 489 "BB02AD.f"
    nsymm = ipar[1] * (ipar[1] + 1) / 2;
#line 490 "BB02AD.f"
    msymm = ipar[2] * (ipar[2] + 1) / 2;
#line 491 "BB02AD.f"
    psymm = ipar[3] * (ipar[3] + 1) / 2;

#line 493 "BB02AD.f"
    dlaset_("A", &ipar[1], &ipar[1], &c_b7, &c_b7, &a[a_offset], lda, (ftnlen)
	    1);
#line 494 "BB02AD.f"
    dlaset_("A", &ipar[1], &ipar[2], &c_b7, &c_b7, &b[b_offset], ldb, (ftnlen)
	    1);
#line 495 "BB02AD.f"
    dlaset_("A", &ipar[3], &ipar[1], &c_b7, &c_b7, &c__[c_offset], ldc, (
	    ftnlen)1);
#line 496 "BB02AD.f"
    dlaset_("L", &psymm, &c__1, &c_b7, &c_b7, &q[1], &c__1, (ftnlen)1);
#line 497 "BB02AD.f"
    dlaset_("L", &msymm, &c__1, &c_b7, &c_b7, &r__[1], &c__1, (ftnlen)1);
#line 498 "BB02AD.f"
    if (bpar[7]) {
#line 498 "BB02AD.f"
	dlaset_("A", &ipar[1], &ipar[2], &c_b7, &c_b7, &s[s_offset], lds, (
		ftnlen)1);
#line 498 "BB02AD.f"
    }

#line 501 "BB02AD.f"
    if (nr[1] == 1) {

#line 503 "BB02AD.f"
	if (nr[2] == 1) {
#line 504 "BB02AD.f"
	    a[a_dim1 + 1] = 2.;
#line 505 "BB02AD.f"
	    a[a_dim1 + 2] = 1.;
#line 506 "BB02AD.f"
	    a[(a_dim1 << 1) + 1] = -1.;
#line 507 "BB02AD.f"
	    b[b_dim1 + 1] = 1.;
#line 508 "BB02AD.f"
	    q[1] = 1.;
#line 509 "BB02AD.f"
	    c__[(c_dim1 << 1) + 1] = 1.;
#line 510 "BB02AD.f"
	    r__[1] = 0.;
#line 511 "BB02AD.f"
	    dlaset_("A", &ipar[1], &ipar[1], &c_b7, &c_b34, &x[x_offset], ldx,
		     (ftnlen)1);
#line 512 "BB02AD.f"
	    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);

#line 514 "BB02AD.f"
	} else if (nr[2] == 2) {
#line 515 "BB02AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 516 "BB02AD.f"
	    a[(a_dim1 << 1) + 2] = -1.;
#line 517 "BB02AD.f"
	    b[b_dim1 + 1] = 1.;
#line 518 "BB02AD.f"
	    b[b_dim1 + 2] = 2.;
#line 519 "BB02AD.f"
	    b[(b_dim1 << 1) + 2] = 1.;
#line 520 "BB02AD.f"
	    r__[1] = 9.;
#line 521 "BB02AD.f"
	    r__[2] = 3.;
#line 522 "BB02AD.f"
	    r__[3] = 1.;
#line 523 "BB02AD.f"
	    dlaset_("A", &psymm, &c__1, &c_b38, &c_b38, &q[1], &psymm, (
		    ftnlen)1);
#line 524 "BB02AD.f"
	    q[3] = 7.;
#line 525 "BB02AD.f"
	    drscl_(&msymm, &c_b40, &q[1], &c__1);
#line 526 "BB02AD.f"
	    if (bpar[7]) {
#line 527 "BB02AD.f"
		s[s_dim1 + 1] = 3.;
#line 528 "BB02AD.f"
		s[s_dim1 + 2] = -1.;
#line 529 "BB02AD.f"
		s[(s_dim1 << 1) + 1] = 1.;
#line 530 "BB02AD.f"
		s[(s_dim1 << 1) + 2] = 7.;
#line 531 "BB02AD.f"
	    }
#line 532 "BB02AD.f"
	    s_copy(ident, "0100", (ftnlen)4, (ftnlen)4);

#line 534 "BB02AD.f"
	} else if (nr[2] == 3) {
#line 535 "BB02AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 536 "BB02AD.f"
	    b[b_dim1 + 2] = 1.;
#line 537 "BB02AD.f"
	    q[1] = 1.;
#line 538 "BB02AD.f"
	    q[2] = 2.;
#line 539 "BB02AD.f"
	    q[3] = 4.;
#line 540 "BB02AD.f"
	    x[x_dim1 + 1] = 1.;
#line 541 "BB02AD.f"
	    x[x_dim1 + 2] = 2.;
#line 542 "BB02AD.f"
	    x[(x_dim1 << 1) + 1] = 2.;
#line 543 "BB02AD.f"
	    x[(x_dim1 << 1) + 2] = sqrt(5.) + 2.;
#line 544 "BB02AD.f"
	    s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);

#line 546 "BB02AD.f"
	} else if (nr[2] == 4) {
#line 547 "BB02AD.f"
	    a[(a_dim1 << 1) + 1] = .1;
#line 548 "BB02AD.f"
	    a[a_dim1 * 3 + 2] = .01;
#line 549 "BB02AD.f"
	    b[b_dim1 + 1] = 1.;
#line 550 "BB02AD.f"
	    b[(b_dim1 << 1) + 3] = 1.;
#line 551 "BB02AD.f"
	    r__[3] = 1.;
#line 552 "BB02AD.f"
	    q[1] = 1e5;
#line 553 "BB02AD.f"
	    q[4] = 1e3;
#line 554 "BB02AD.f"
	    q[6] = -10.;
#line 555 "BB02AD.f"
	    x[x_dim1 + 1] = 1e5;
#line 556 "BB02AD.f"
	    x[(x_dim1 << 1) + 2] = 1e3;
#line 557 "BB02AD.f"
	    s_copy(ident, "0100", (ftnlen)4, (ftnlen)4);

#line 559 "BB02AD.f"
	} else if (nr[2] >= 5 && nr[2] <= 8 || nr[2] == 10 || nr[2] == 11 || 
		nr[2] == 13) {
#line 562 "BB02AD.f"
	    if (nr[2] < 10) {
#line 563 "BB02AD.f"
		ici__1.icierr = 0;
#line 563 "BB02AD.f"
		ici__1.icirnum = 1;
#line 563 "BB02AD.f"
		ici__1.icirlen = 11;
#line 563 "BB02AD.f"
		ici__1.iciunit = chpar;
#line 563 "BB02AD.f"
		ici__1.icifmt = "(A,I1,A,I1,A)";
#line 563 "BB02AD.f"
		s_wsfi(&ici__1);
#line 563 "BB02AD.f"
		do_fio(&c__1, "BB02", (ftnlen)4);
#line 563 "BB02AD.f"
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
#line 563 "BB02AD.f"
		do_fio(&c__1, "0", (ftnlen)1);
#line 563 "BB02AD.f"
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
#line 563 "BB02AD.f"
		do_fio(&c__1, ".dat", (ftnlen)4);
#line 563 "BB02AD.f"
		e_wsfi();
#line 565 "BB02AD.f"
		o__1.oerr = 1;
#line 565 "BB02AD.f"
		o__1.ounit = 1;
#line 565 "BB02AD.f"
		o__1.ofnmlen = 11;
#line 565 "BB02AD.f"
		o__1.ofnm = chpar;
#line 565 "BB02AD.f"
		o__1.orl = 0;
#line 565 "BB02AD.f"
		o__1.osta = "OLD";
#line 565 "BB02AD.f"
		o__1.oacc = 0;
#line 565 "BB02AD.f"
		o__1.ofm = 0;
#line 565 "BB02AD.f"
		o__1.oblnk = 0;
#line 565 "BB02AD.f"
		ios = f_open(&o__1);
#line 566 "BB02AD.f"
	    } else {
#line 567 "BB02AD.f"
		ici__1.icierr = 0;
#line 567 "BB02AD.f"
		ici__1.icirnum = 1;
#line 567 "BB02AD.f"
		ici__1.icirlen = 11;
#line 567 "BB02AD.f"
		ici__1.iciunit = chpar;
#line 567 "BB02AD.f"
		ici__1.icifmt = "(A,I1,I2,A)";
#line 567 "BB02AD.f"
		s_wsfi(&ici__1);
#line 567 "BB02AD.f"
		do_fio(&c__1, "BB02", (ftnlen)4);
#line 567 "BB02AD.f"
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
#line 567 "BB02AD.f"
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
#line 567 "BB02AD.f"
		do_fio(&c__1, ".dat", (ftnlen)4);
#line 567 "BB02AD.f"
		e_wsfi();
#line 569 "BB02AD.f"
		o__1.oerr = 1;
#line 569 "BB02AD.f"
		o__1.ounit = 1;
#line 569 "BB02AD.f"
		o__1.ofnmlen = 11;
#line 569 "BB02AD.f"
		o__1.ofnm = chpar;
#line 569 "BB02AD.f"
		o__1.orl = 0;
#line 569 "BB02AD.f"
		o__1.osta = "OLD";
#line 569 "BB02AD.f"
		o__1.oacc = 0;
#line 569 "BB02AD.f"
		o__1.ofm = 0;
#line 569 "BB02AD.f"
		o__1.oblnk = 0;
#line 569 "BB02AD.f"
		ios = f_open(&o__1);
#line 570 "BB02AD.f"
	    }
#line 571 "BB02AD.f"
	    if (ios != 0) {
#line 572 "BB02AD.f"
		*info = 1;
#line 573 "BB02AD.f"
	    } else {
#line 574 "BB02AD.f"
		if (! (nr[2] == 13)) {
#line 575 "BB02AD.f"
		    i__1 = ipar[1];
#line 575 "BB02AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 576 "BB02AD.f"
			ios = s_rsle(&io___12);
#line 576 "BB02AD.f"
			if (ios != 0) {
#line 576 "BB02AD.f"
			    goto L100001;
#line 576 "BB02AD.f"
			}
#line 576 "BB02AD.f"
			i__2 = ipar[1];
#line 576 "BB02AD.f"
			for (j = 1; j <= i__2; ++j) {
#line 576 "BB02AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				    a_dim1], (ftnlen)sizeof(doublereal));
#line 576 "BB02AD.f"
			    if (ios != 0) {
#line 576 "BB02AD.f"
				goto L100001;
#line 576 "BB02AD.f"
			    }
#line 576 "BB02AD.f"
			}
#line 576 "BB02AD.f"
			ios = e_rsle();
#line 576 "BB02AD.f"
L100001:
#line 577 "BB02AD.f"
			if (ios != 0) {
#line 577 "BB02AD.f"
			    *info = 1;
#line 577 "BB02AD.f"
			}
#line 578 "BB02AD.f"
/* L10: */
#line 578 "BB02AD.f"
		    }
#line 579 "BB02AD.f"
		    i__1 = ipar[1];
#line 579 "BB02AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 580 "BB02AD.f"
			ios = s_rsle(&io___14);
#line 580 "BB02AD.f"
			if (ios != 0) {
#line 580 "BB02AD.f"
			    goto L100002;
#line 580 "BB02AD.f"
			}
#line 580 "BB02AD.f"
			i__2 = ipar[2];
#line 580 "BB02AD.f"
			for (j = 1; j <= i__2; ++j) {
#line 580 "BB02AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				    b_dim1], (ftnlen)sizeof(doublereal));
#line 580 "BB02AD.f"
			    if (ios != 0) {
#line 580 "BB02AD.f"
				goto L100002;
#line 580 "BB02AD.f"
			    }
#line 580 "BB02AD.f"
			}
#line 580 "BB02AD.f"
			ios = e_rsle();
#line 580 "BB02AD.f"
L100002:
#line 581 "BB02AD.f"
			if (ios != 0) {
#line 581 "BB02AD.f"
			    *info = 1;
#line 581 "BB02AD.f"
			}
#line 582 "BB02AD.f"
/* L20: */
#line 582 "BB02AD.f"
		    }
#line 583 "BB02AD.f"
		}
#line 584 "BB02AD.f"
		if (nr[2] == 5) {
#line 585 "BB02AD.f"
		    q[1] = 1.87;
#line 586 "BB02AD.f"
		    q[4] = -.244;
#line 587 "BB02AD.f"
		    q[5] = .744;
#line 588 "BB02AD.f"
		    q[6] = .205;
#line 589 "BB02AD.f"
		    q[8] = .589;
#line 590 "BB02AD.f"
		    q[10] = 1.048;
#line 591 "BB02AD.f"
		} else if (nr[2] == 6) {
#line 592 "BB02AD.f"
		    q[1] = .01;
#line 593 "BB02AD.f"
		    q[5] = .01;
#line 594 "BB02AD.f"
		    q[8] = .01;
#line 595 "BB02AD.f"
		    q[10] = .01;
#line 596 "BB02AD.f"
		} else if (nr[2] == 7) {
#line 597 "BB02AD.f"
		    dlaset_("U", &ipar[3], &ipar[1], &c_b34, &c_b34, &c__[
			    c_offset], ldc, (ftnlen)1);
#line 598 "BB02AD.f"
		    c__[c_dim1 * 3 + 1] = 2.;
#line 599 "BB02AD.f"
		    c__[(c_dim1 << 2) + 1] = 4.;
#line 600 "BB02AD.f"
		    c__[(c_dim1 << 2) + 2] = 2.;
#line 601 "BB02AD.f"
		    q[1] = 2.;
#line 602 "BB02AD.f"
		    q[2] = -1.;
#line 603 "BB02AD.f"
		    q[5] = 2.;
#line 604 "BB02AD.f"
		    q[6] = -1.;
#line 605 "BB02AD.f"
		    q[8] = 2.;
#line 606 "BB02AD.f"
		} else if (nr[2] == 10) {
#line 607 "BB02AD.f"
		    c__[c_dim1 + 1] = 1.;
#line 608 "BB02AD.f"
		    c__[c_dim1 * 5 + 2] = 1.;
#line 609 "BB02AD.f"
		    q[1] = 50.;
#line 610 "BB02AD.f"
		    q[3] = 50.;
#line 611 "BB02AD.f"
		} else if (nr[2] == 11) {
#line 612 "BB02AD.f"
		    a[a_dim1 * 10 + 10] = 1.;
#line 613 "BB02AD.f"
		    a[a_dim1 * 11 + 11] = 1.;
#line 614 "BB02AD.f"
		    c__[c_dim1 * 6 + 1] = 15.;
#line 615 "BB02AD.f"
		    c__[c_dim1 * 7 + 2] = 7.;
#line 616 "BB02AD.f"
		    c__[(c_dim1 << 3) + 2] = -5.357;
#line 617 "BB02AD.f"
		    c__[c_dim1 * 9 + 2] = -3.943;
#line 618 "BB02AD.f"
		    c__[c_dim1 * 10 + 3] = 1.;
#line 619 "BB02AD.f"
		    c__[c_dim1 * 11 + 4] = 1.;
#line 620 "BB02AD.f"
		    q[1] = .5;
#line 621 "BB02AD.f"
		    q[5] = 5.;
#line 622 "BB02AD.f"
		    q[8] = .5;
#line 623 "BB02AD.f"
		    q[10] = 5.;
#line 624 "BB02AD.f"
		    r__[1] = 400.;
#line 625 "BB02AD.f"
		    r__[3] = 700.;
#line 626 "BB02AD.f"
		    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);

#line 628 "BB02AD.f"
		} else if (nr[2] == 13) {
#line 629 "BB02AD.f"
		    i__1 = ipar[1] - 6;
#line 629 "BB02AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 630 "BB02AD.f"
			ios = s_rsle(&io___15);
#line 630 "BB02AD.f"
			if (ios != 0) {
#line 630 "BB02AD.f"
			    goto L100003;
#line 630 "BB02AD.f"
			}
#line 630 "BB02AD.f"
			i__2 = ipar[1] - 6;
#line 630 "BB02AD.f"
			for (j = 1; j <= i__2; ++j) {
#line 630 "BB02AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				    a_dim1], (ftnlen)sizeof(doublereal));
#line 630 "BB02AD.f"
			    if (ios != 0) {
#line 630 "BB02AD.f"
				goto L100003;
#line 630 "BB02AD.f"
			    }
#line 630 "BB02AD.f"
			}
#line 630 "BB02AD.f"
			ios = e_rsle();
#line 630 "BB02AD.f"
L100003:
#line 632 "BB02AD.f"
			if (ios != 0) {
#line 632 "BB02AD.f"
			    *info = 1;
#line 632 "BB02AD.f"
			}
#line 633 "BB02AD.f"
/* L24: */
#line 633 "BB02AD.f"
		    }
#line 634 "BB02AD.f"
		    i__1 = ipar[1] - 6;
#line 634 "BB02AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 635 "BB02AD.f"
			ios = s_rsle(&io___16);
#line 635 "BB02AD.f"
			if (ios != 0) {
#line 635 "BB02AD.f"
			    goto L100004;
#line 635 "BB02AD.f"
			}
#line 635 "BB02AD.f"
			i__2 = ipar[2];
#line 635 "BB02AD.f"
			for (j = 1; j <= i__2; ++j) {
#line 635 "BB02AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				    b_dim1], (ftnlen)sizeof(doublereal));
#line 635 "BB02AD.f"
			    if (ios != 0) {
#line 635 "BB02AD.f"
				goto L100004;
#line 635 "BB02AD.f"
			    }
#line 635 "BB02AD.f"
			}
#line 635 "BB02AD.f"
			ios = e_rsle();
#line 635 "BB02AD.f"
L100004:
#line 637 "BB02AD.f"
			if (ios != 0) {
#line 637 "BB02AD.f"
			    *info = 1;
#line 637 "BB02AD.f"
			}
#line 638 "BB02AD.f"
/* L25: */
#line 638 "BB02AD.f"
		    }
#line 639 "BB02AD.f"
		    i__1 = ipar[2];
#line 639 "BB02AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 640 "BB02AD.f"
			ios = s_rsle(&io___17);
#line 640 "BB02AD.f"
			if (ios != 0) {
#line 640 "BB02AD.f"
			    goto L100005;
#line 640 "BB02AD.f"
			}
#line 640 "BB02AD.f"
			i__2 = ipar[1] - 6;
#line 640 "BB02AD.f"
			for (j = 1; j <= i__2; ++j) {
#line 640 "BB02AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&c__[i__ + j * 
				    c_dim1], (ftnlen)sizeof(doublereal));
#line 640 "BB02AD.f"
			    if (ios != 0) {
#line 640 "BB02AD.f"
				goto L100005;
#line 640 "BB02AD.f"
			    }
#line 640 "BB02AD.f"
			}
#line 640 "BB02AD.f"
			ios = e_rsle();
#line 640 "BB02AD.f"
L100005:
#line 642 "BB02AD.f"
			if (ios != 0) {
#line 642 "BB02AD.f"
			    *info = 1;
#line 642 "BB02AD.f"
			}
#line 643 "BB02AD.f"
/* L26: */
#line 643 "BB02AD.f"
		    }
#line 644 "BB02AD.f"
		    for (i__ = 1; i__ <= 6; ++i__) {
#line 645 "BB02AD.f"
			a[i__ + 20 + (i__ + 20) * a_dim1] = 1.;
#line 646 "BB02AD.f"
			c__[i__ + 6 + (i__ + 20) * c_dim1] = 1.;
#line 647 "BB02AD.f"
/* L27: */
#line 647 "BB02AD.f"
		    }
#line 648 "BB02AD.f"
		    j = 58;
#line 649 "BB02AD.f"
		    for (i__ = 7; i__ <= 12; ++i__) {
#line 650 "BB02AD.f"
			ios = s_rsle(&io___18);
#line 650 "BB02AD.f"
			if (ios != 0) {
#line 650 "BB02AD.f"
			    goto L100006;
#line 650 "BB02AD.f"
			}
#line 650 "BB02AD.f"
			ios = do_lio(&c__5, &c__1, (char *)&q[j], (ftnlen)
				sizeof(doublereal));
#line 650 "BB02AD.f"
			if (ios != 0) {
#line 650 "BB02AD.f"
			    goto L100006;
#line 650 "BB02AD.f"
			}
#line 650 "BB02AD.f"
			ios = e_rsle();
#line 650 "BB02AD.f"
L100006:
#line 651 "BB02AD.f"
			if (ios != 0) {
#line 651 "BB02AD.f"
			    *info = 1;
#line 651 "BB02AD.f"
			}
#line 652 "BB02AD.f"
			j += 13 - i__;
#line 653 "BB02AD.f"
/* L28: */
#line 653 "BB02AD.f"
		    }
#line 654 "BB02AD.f"
		    j = 1;
#line 655 "BB02AD.f"
		    for (i__ = 1; i__ <= 6; ++i__) {
#line 656 "BB02AD.f"
			ios = s_rsle(&io___19);
#line 656 "BB02AD.f"
			if (ios != 0) {
#line 656 "BB02AD.f"
			    goto L100007;
#line 656 "BB02AD.f"
			}
#line 656 "BB02AD.f"
			ios = do_lio(&c__5, &c__1, (char *)&r__[j], (ftnlen)
				sizeof(doublereal));
#line 656 "BB02AD.f"
			if (ios != 0) {
#line 656 "BB02AD.f"
			    goto L100007;
#line 656 "BB02AD.f"
			}
#line 656 "BB02AD.f"
			ios = e_rsle();
#line 656 "BB02AD.f"
L100007:
#line 657 "BB02AD.f"
			if (ios != 0) {
#line 657 "BB02AD.f"
			    *info = 1;
#line 657 "BB02AD.f"
			}
#line 658 "BB02AD.f"
			j += 7 - i__;
#line 659 "BB02AD.f"
/* L29: */
#line 659 "BB02AD.f"
		    }
#line 660 "BB02AD.f"
		    for (i__ = 1; i__ <= 6; ++i__) {
#line 661 "BB02AD.f"
			for (j = 1; j <= 20; ++j) {
#line 662 "BB02AD.f"
			    a[i__ + 20 + j * a_dim1] = -c__[i__ + j * c_dim1];
#line 663 "BB02AD.f"
/* L30: */
#line 663 "BB02AD.f"
			}
#line 664 "BB02AD.f"
/* L31: */
#line 664 "BB02AD.f"
		    }
#line 665 "BB02AD.f"
		    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);
#line 666 "BB02AD.f"
		}
#line 667 "BB02AD.f"
	    }
#line 668 "BB02AD.f"
	    cl__1.cerr = 0;
#line 668 "BB02AD.f"
	    cl__1.cunit = 1;
#line 668 "BB02AD.f"
	    cl__1.csta = 0;
#line 668 "BB02AD.f"
	    f_clos(&cl__1);
#line 669 "BB02AD.f"
	    if (nr[2] == 5 || nr[2] == 6) {
#line 670 "BB02AD.f"
		s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
#line 671 "BB02AD.f"
	    } else if (nr[2] == 7 || nr[2] == 10) {
#line 672 "BB02AD.f"
		s_copy(ident, "0001", (ftnlen)4, (ftnlen)4);
#line 673 "BB02AD.f"
	    } else if (nr[2] == 8) {
#line 674 "BB02AD.f"
		s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);
#line 675 "BB02AD.f"
	    }

#line 677 "BB02AD.f"
	} else if (nr[2] == 9) {
#line 678 "BB02AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 679 "BB02AD.f"
	    a[a_dim1 * 3 + 2] = 1.;
#line 680 "BB02AD.f"
	    a[a_dim1 * 5 + 4] = 1.;
#line 681 "BB02AD.f"
	    a[a_dim1 * 6 + 5] = 1.;
#line 682 "BB02AD.f"
	    b[b_dim1 + 3] = 1.;
#line 683 "BB02AD.f"
	    b[(b_dim1 << 1) + 6] = 1.;
#line 684 "BB02AD.f"
	    c__[c_dim1 + 1] = 1.;
#line 685 "BB02AD.f"
	    c__[(c_dim1 << 1) + 1] = 1.;
#line 686 "BB02AD.f"
	    c__[(c_dim1 << 2) + 2] = 1.;
#line 687 "BB02AD.f"
	    c__[c_dim1 * 5 + 2] = -1.;
#line 688 "BB02AD.f"
	    r__[1] = 3.;
#line 689 "BB02AD.f"
	    r__[3] = 1.;
#line 690 "BB02AD.f"
	    if (bpar[7]) {
#line 691 "BB02AD.f"
		s[s_dim1 + 1] = 1.;
#line 692 "BB02AD.f"
		s[s_dim1 + 2] = 1.;
#line 693 "BB02AD.f"
		s[s_dim1 + 4] = 1.;
#line 694 "BB02AD.f"
		s[s_dim1 + 5] = -1.;
#line 695 "BB02AD.f"
	    }
#line 696 "BB02AD.f"
	    s_copy(ident, "0010", (ftnlen)4, (ftnlen)4);
#line 697 "BB02AD.f"
	} else if (nr[2] == 12) {
#line 698 "BB02AD.f"
	    for (i__ = 1; i__ <= 10; ++i__) {
#line 699 "BB02AD.f"
		a[i__ + (i__ + 1) * a_dim1] = 1.;
#line 700 "BB02AD.f"
/* L32: */
#line 700 "BB02AD.f"
	    }
#line 701 "BB02AD.f"
	    a[a_dim1 * 7 + 6] = 0.;
#line 702 "BB02AD.f"
	    a[a_dim1 * 9 + 8] = 0.;
#line 703 "BB02AD.f"
	    a[a_dim1 * 12 + 12] = 1.;
#line 704 "BB02AD.f"
	    a[a_dim1 * 13 + 13] = 1.;
#line 705 "BB02AD.f"
	    a[a_dim1 + 12] = -3.318;
#line 706 "BB02AD.f"
	    a[a_dim1 + 13] = -1.5484;
#line 707 "BB02AD.f"
	    a[a_dim1 * 6 + 6] = .7788;
#line 708 "BB02AD.f"
	    a[a_dim1 * 7 + 8] = -.4724;
#line 709 "BB02AD.f"
	    a[a_dim1 * 7 + 13] = .3981;
#line 710 "BB02AD.f"
	    a[(a_dim1 << 3) + 8] = 1.3746;
#line 711 "BB02AD.f"
	    a[(a_dim1 << 3) + 13] = .5113;
#line 712 "BB02AD.f"
	    a[a_dim1 * 9 + 13] = 5.7865;
#line 713 "BB02AD.f"
	    a[a_dim1 * 11 + 11] = .8071;
#line 714 "BB02AD.f"
	    b[b_dim1 + 6] = 1.;
#line 715 "BB02AD.f"
	    b[(b_dim1 << 1) + 8] = 1.;
#line 716 "BB02AD.f"
	    c__[c_dim1 + 1] = 3.318;
#line 717 "BB02AD.f"
	    c__[c_dim1 + 2] = 1.5484;
#line 718 "BB02AD.f"
	    c__[c_dim1 * 7 + 2] = -.3981;
#line 719 "BB02AD.f"
	    c__[(c_dim1 << 3) + 2] = -.5113;
#line 720 "BB02AD.f"
	    c__[c_dim1 * 9 + 2] = -5.7865;
#line 721 "BB02AD.f"
	    c__[c_dim1 * 12 + 3] = 1.;
#line 722 "BB02AD.f"
	    c__[c_dim1 * 13 + 4] = 1.;
#line 723 "BB02AD.f"
	    q[1] = .5;
#line 724 "BB02AD.f"
	    q[5] = 5.;
#line 725 "BB02AD.f"
	    q[8] = .5;
#line 726 "BB02AD.f"
	    q[10] = 5.;
#line 727 "BB02AD.f"
	    r__[1] = 400.;
#line 728 "BB02AD.f"
	    r__[3] = 700.;
#line 729 "BB02AD.f"
	    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);
#line 730 "BB02AD.f"
	}

#line 732 "BB02AD.f"
    } else if (nr[1] == 2) {
#line 733 "BB02AD.f"
	if (nr[2] == 1) {
#line 734 "BB02AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 734 "BB02AD.f"
		dpar[1] = 1e6;
#line 734 "BB02AD.f"
	    }
#line 735 "BB02AD.f"
	    a[a_dim1 + 1] = 4.;
#line 736 "BB02AD.f"
	    a[a_dim1 + 2] = -4.5;
#line 737 "BB02AD.f"
	    a[(a_dim1 << 1) + 1] = 3.;
#line 738 "BB02AD.f"
	    a[(a_dim1 << 1) + 2] = -3.5;
#line 739 "BB02AD.f"
	    dlaset_("A", &ipar[1], &ipar[2], &c_b112, &c_b34, &b[b_offset], 
		    ldb, (ftnlen)1);
#line 740 "BB02AD.f"
	    r__[1] = dpar[1];
#line 741 "BB02AD.f"
	    q[1] = 9.;
#line 742 "BB02AD.f"
	    q[2] = 6.;
#line 743 "BB02AD.f"
	    q[3] = 4.;
#line 744 "BB02AD.f"
	    temp = (sqrt(dpar[1] * 4. + 1.) + 1.) / 2.;
#line 745 "BB02AD.f"
	    x[x_dim1 + 1] = temp * q[1];
#line 746 "BB02AD.f"
	    x[x_dim1 + 2] = temp * q[2];
#line 747 "BB02AD.f"
	    x[(x_dim1 << 1) + 1] = x[x_dim1 + 2];
#line 748 "BB02AD.f"
	    x[(x_dim1 << 1) + 2] = temp * q[3];
#line 749 "BB02AD.f"
	    s_copy(ident, "0100", (ftnlen)4, (ftnlen)4);

#line 751 "BB02AD.f"
	} else if (nr[2] == 2) {
#line 752 "BB02AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 752 "BB02AD.f"
		dpar[1] = 1e6;
#line 752 "BB02AD.f"
	    }
#line 753 "BB02AD.f"
	    if (dpar[1] == 0.) {
#line 754 "BB02AD.f"
		*info = 2;
#line 755 "BB02AD.f"
	    } else {
#line 756 "BB02AD.f"
		a[a_dim1 + 1] = .9512;
#line 757 "BB02AD.f"
		a[(a_dim1 << 1) + 2] = .9048;
#line 758 "BB02AD.f"
		dlaset_("A", &c__1, &ipar[2], &c_b118, &c_b118, &b[b_offset], 
			ldb, (ftnlen)1);
#line 759 "BB02AD.f"
		b[b_dim1 + 2] = -1.1895;
#line 760 "BB02AD.f"
		b[(b_dim1 << 1) + 2] = 3.569;
#line 761 "BB02AD.f"
		r__[1] = 1. / (dpar[1] * 3.);
#line 762 "BB02AD.f"
		r__[3] = dpar[1] * 3.;
#line 763 "BB02AD.f"
		q[1] = .005;
#line 764 "BB02AD.f"
		q[3] = .02;
#line 765 "BB02AD.f"
		s_copy(ident, "0100", (ftnlen)4, (ftnlen)4);
#line 766 "BB02AD.f"
	    }

#line 768 "BB02AD.f"
	} else if (nr[2] == 3) {
#line 769 "BB02AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 769 "BB02AD.f"
		dpar[1] = 1e6;
#line 769 "BB02AD.f"
	    }
#line 770 "BB02AD.f"
	    a[(a_dim1 << 1) + 1] = dpar[1];
#line 771 "BB02AD.f"
	    b[b_dim1 + 2] = 1.;
#line 772 "BB02AD.f"
	    x[x_dim1 + 1] = 1.;
#line 773 "BB02AD.f"
	    x[(x_dim1 << 1) + 2] = dpar[1] * dpar[1] + 1.;
#line 774 "BB02AD.f"
	    s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);

#line 776 "BB02AD.f"
	} else if (nr[2] == 4) {
#line 777 "BB02AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 777 "BB02AD.f"
		dpar[1] = 1e6;
#line 777 "BB02AD.f"
	    }
#line 778 "BB02AD.f"
	    a[(a_dim1 << 1) + 2] = 1.;
#line 779 "BB02AD.f"
	    a[a_dim1 * 3 + 3] = 3.;
#line 780 "BB02AD.f"
	    r__[1] = dpar[1];
#line 781 "BB02AD.f"
	    r__[4] = dpar[1];
#line 782 "BB02AD.f"
	    r__[6] = dpar[1];
/*     .. set C = V .. */
#line 784 "BB02AD.f"
	    temp = .66666666666666663;
#line 785 "BB02AD.f"
	    d__1 = -temp;
#line 785 "BB02AD.f"
	    d__2 = 1. - temp;
#line 785 "BB02AD.f"
	    dlaset_("A", &ipar[3], &ipar[1], &d__1, &d__2, &c__[c_offset], 
		    ldc, (ftnlen)1);
/*     .. and compute A <- C' A C */
#line 787 "BB02AD.f"
	    dsymm_("L", "L", &ipar[1], &ipar[1], &c_b34, &c__[c_offset], ldc, 
		    &a[a_offset], lda, &c_b7, &dwork[1], &ipar[1], (ftnlen)1, 
		    (ftnlen)1);
#line 789 "BB02AD.f"
	    dsymm_("R", "L", &ipar[1], &ipar[1], &c_b34, &c__[c_offset], ldc, 
		    &dwork[1], &ipar[1], &c_b7, &a[a_offset], lda, (ftnlen)1, 
		    (ftnlen)1);
#line 791 "BB02AD.f"
	    q[1] = dpar[1];
#line 792 "BB02AD.f"
	    q[4] = dpar[1];
#line 793 "BB02AD.f"
	    q[6] = dpar[1];
#line 794 "BB02AD.f"
	    x[x_dim1 + 1] = dpar[1];
#line 795 "BB02AD.f"
	    x[(x_dim1 << 1) + 2] = dpar[1] * (sqrt(5.) + 1.) / 2.;
#line 796 "BB02AD.f"
	    x[x_dim1 * 3 + 3] = dpar[1] * (sqrt(85.) + 9.) / 2.;
#line 797 "BB02AD.f"
	    dsymm_("L", "L", &ipar[1], &ipar[1], &c_b34, &c__[c_offset], ldc, 
		    &x[x_offset], ldx, &c_b7, &dwork[1], &ipar[1], (ftnlen)1, 
		    (ftnlen)1);
#line 799 "BB02AD.f"
	    dsymm_("R", "L", &ipar[1], &ipar[1], &c_b34, &c__[c_offset], ldc, 
		    &dwork[1], &ipar[1], &c_b7, &x[x_offset], ldx, (ftnlen)1, 
		    (ftnlen)1);
#line 801 "BB02AD.f"
	    s_copy(ident, "1000", (ftnlen)4, (ftnlen)4);

#line 803 "BB02AD.f"
	} else if (nr[2] == 5) {
#line 804 "BB02AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 805 "BB02AD.f"
		dpar[4] = .25;
#line 806 "BB02AD.f"
		dpar[3] = 1.;
#line 807 "BB02AD.f"
		dpar[2] = 1.;
#line 808 "BB02AD.f"
		dpar[1] = 1e8;
#line 809 "BB02AD.f"
	    }
#line 810 "BB02AD.f"
	    if (dpar[1] == 0.) {
#line 811 "BB02AD.f"
		*info = 2;
#line 812 "BB02AD.f"
	    } else {
#line 813 "BB02AD.f"
		temp = dpar[2] / dpar[1];
#line 814 "BB02AD.f"
		beta = dpar[3] * temp;
#line 815 "BB02AD.f"
		alpha = 1. - temp;
#line 816 "BB02AD.f"
		a[a_dim1 + 1] = alpha;
#line 817 "BB02AD.f"
		i__1 = ipar[1] - 1;
#line 817 "BB02AD.f"
		i__2 = ipar[1] - 1;
#line 817 "BB02AD.f"
		dlaset_("A", &i__1, &i__2, &c_b7, &c_b34, &a[a_dim1 + 2], lda,
			 (ftnlen)1);
#line 819 "BB02AD.f"
		b[b_dim1 + 1] = beta;
#line 820 "BB02AD.f"
		c__[(c_dim1 << 2) + 1] = 1.;
#line 821 "BB02AD.f"
		r__[1] = dpar[4];
#line 822 "BB02AD.f"
		if (beta == 0.) {
#line 823 "BB02AD.f"
		    *info = 2;
#line 824 "BB02AD.f"
		} else {
#line 825 "BB02AD.f"
		    dlaset_("A", &ipar[1], &ipar[1], &c_b7, &c_b34, &x[
			    x_offset], ldx, (ftnlen)1);
#line 826 "BB02AD.f"
		    beta *= beta;
#line 827 "BB02AD.f"
		    temp = dpar[4] * (alpha + 1.) * (alpha - 1.) + beta;
#line 828 "BB02AD.f"
		    x[x_dim1 + 1] = temp + sqrt(temp * temp + beta * 4. * 
			    dpar[4]);
#line 829 "BB02AD.f"
		    x[x_dim1 + 1] = x[x_dim1 + 1] / 2. / beta;
#line 830 "BB02AD.f"
		}
#line 831 "BB02AD.f"
		s_copy(ident, "0010", (ftnlen)4, (ftnlen)4);
#line 832 "BB02AD.f"
	    }
#line 833 "BB02AD.f"
	}

#line 835 "BB02AD.f"
    } else if (nr[1] == 4) {
#line 836 "BB02AD.f"
	if (nr[2] == 1) {
#line 837 "BB02AD.f"
	    if (lsame_(def, "D", (ftnlen)1, (ftnlen)1)) {
#line 837 "BB02AD.f"
		dpar[1] = 1.;
#line 837 "BB02AD.f"
	    }
#line 838 "BB02AD.f"
	    i__1 = ipar[1] - 1;
#line 838 "BB02AD.f"
	    i__2 = ipar[1] - 1;
#line 838 "BB02AD.f"
	    dlaset_("A", &i__1, &i__2, &c_b7, &c_b34, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
#line 839 "BB02AD.f"
	    b[ipar[1] + b_dim1] = 1.;
#line 840 "BB02AD.f"
	    r__[1] = dpar[1];
#line 841 "BB02AD.f"
	    i__1 = ipar[1];
#line 841 "BB02AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 842 "BB02AD.f"
		x[i__ + i__ * x_dim1] = (doublereal) i__;
#line 843 "BB02AD.f"
/* L40: */
#line 843 "BB02AD.f"
	    }
#line 844 "BB02AD.f"
	    s_copy(ident, "0110", (ftnlen)4, (ftnlen)4);
#line 845 "BB02AD.f"
	}
#line 846 "BB02AD.f"
    }

#line 848 "BB02AD.f"
    if (*info != 0) {
#line 848 "BB02AD.f"
	goto L2001;
#line 848 "BB02AD.f"
    }
/*     .. set up data in required format .. */

#line 851 "BB02AD.f"
    if (bpar[4]) {
/*     .. G is to be returned in product form .. */
#line 853 "BB02AD.f"
	rdimm = ipar[1];
#line 854 "BB02AD.f"
	if (*(unsigned char *)&ident[3] == '0') {
/*       .. invert R using Cholesky factorization, .. */
#line 856 "BB02AD.f"
	    dpptrf_("L", &ipar[2], &r__[1], info, (ftnlen)1);
#line 857 "BB02AD.f"
	    if (*info == 0) {
#line 858 "BB02AD.f"
		dpptri_("L", &ipar[2], &r__[1], info, (ftnlen)1);
#line 859 "BB02AD.f"
		if (*(unsigned char *)ident == '0') {
/*           .. B is not identity matrix .. */
#line 861 "BB02AD.f"
		    i__1 = ipar[1];
#line 861 "BB02AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 862 "BB02AD.f"
			dspmv_("L", &ipar[2], &c_b34, &r__[1], &b[i__ + 
				b_dim1], ldb, &c_b7, &dwork[(i__ - 1) * ipar[
				1] + 1], &c__1, (ftnlen)1);
#line 864 "BB02AD.f"
/* L100: */
#line 864 "BB02AD.f"
		    }
#line 865 "BB02AD.f"
		    dgemv_("T", &ipar[2], &ipar[1], &c_b34, &dwork[1], &ipar[
			    1], &b[b_dim1 + 1], ldb, &c_b7, &r__[1], &c__1, (
			    ftnlen)1);
#line 867 "BB02AD.f"
		    isymm = ipar[1] + 1;
#line 868 "BB02AD.f"
		    i__1 = ipar[1];
#line 868 "BB02AD.f"
		    for (i__ = 2; i__ <= i__1; ++i__) {
#line 869 "BB02AD.f"
			dgemv_("T", &ipar[2], &ipar[1], &c_b34, &dwork[1], &
				ipar[1], &b[i__ + b_dim1], ldb, &c_b7, &b[
				b_dim1 + 1], ldb, (ftnlen)1);
#line 871 "BB02AD.f"
			i__2 = ipar[1] - i__ + 1;
#line 871 "BB02AD.f"
			dcopy_(&i__2, &b[i__ * b_dim1 + 1], ldb, &r__[isymm], 
				&c__1);
#line 872 "BB02AD.f"
			isymm += ipar[1] - i__ + 1;
#line 873 "BB02AD.f"
/* L110: */
#line 873 "BB02AD.f"
		    }
#line 874 "BB02AD.f"
		}
#line 875 "BB02AD.f"
	    } else {
#line 876 "BB02AD.f"
		if (*info > 0) {
#line 877 "BB02AD.f"
		    *info = 3;
#line 878 "BB02AD.f"
		    goto L2001;
#line 879 "BB02AD.f"
		}
#line 880 "BB02AD.f"
	    }
#line 881 "BB02AD.f"
	} else {
/*       .. R = identity .. */
#line 883 "BB02AD.f"
	    if (*(unsigned char *)ident == '0') {
/*         .. B not identity matrix .. */
#line 885 "BB02AD.f"
		if (ipar[2] == 1) {
#line 886 "BB02AD.f"
		    dlaset_("L", &nsymm, &c__1, &c_b7, &c_b7, &r__[1], &c__1, 
			    (ftnlen)1);
#line 887 "BB02AD.f"
		    dspr_("L", &ipar[1], &c_b34, &b[b_offset], &c__1, &r__[1],
			     (ftnlen)1);
#line 888 "BB02AD.f"
		} else {
#line 889 "BB02AD.f"
		    dsyrk_("L", "N", &ipar[1], &ipar[2], &c_b34, &b[b_offset],
			     ldb, &c_b7, &dwork[1], &ipar[1], (ftnlen)1, (
			    ftnlen)1);
#line 891 "BB02AD.f"
		    ma02dd_("Pack", "Lower", &ipar[1], &dwork[1], &ipar[1], &
			    r__[1], (ftnlen)4, (ftnlen)5);
#line 892 "BB02AD.f"
		}
#line 893 "BB02AD.f"
	    } else {
/*         .. B = R = identity .. */
#line 895 "BB02AD.f"
		isymm = 1;
#line 896 "BB02AD.f"
		for (i__ = ipar[1]; i__ >= 1; --i__) {
#line 897 "BB02AD.f"
		    r__[isymm] = 1.;
#line 898 "BB02AD.f"
		    isymm += i__;
#line 899 "BB02AD.f"
/* L120: */
#line 899 "BB02AD.f"
		}
#line 900 "BB02AD.f"
	    }
#line 901 "BB02AD.f"
	}
#line 902 "BB02AD.f"
    } else {
#line 903 "BB02AD.f"
	rdimm = ipar[2];
#line 904 "BB02AD.f"
	if (*(unsigned char *)ident == '1') {
#line 904 "BB02AD.f"
	    dlaset_("A", &ipar[1], &ipar[2], &c_b7, &c_b34, &b[b_offset], ldb,
		     (ftnlen)1);
#line 904 "BB02AD.f"
	}
#line 906 "BB02AD.f"
	if (*(unsigned char *)&ident[3] == '1') {
#line 907 "BB02AD.f"
	    isymm = 1;
#line 908 "BB02AD.f"
	    for (i__ = ipar[2]; i__ >= 1; --i__) {
#line 909 "BB02AD.f"
		r__[isymm] = 1.;
#line 910 "BB02AD.f"
		isymm += i__;
#line 911 "BB02AD.f"
/* L130: */
#line 911 "BB02AD.f"
	    }
#line 912 "BB02AD.f"
	}
#line 913 "BB02AD.f"
    }

#line 915 "BB02AD.f"
    if (bpar[1]) {
/*     .. Q is to be returned in product form .. */
#line 917 "BB02AD.f"
	qdimm = ipar[1];
#line 918 "BB02AD.f"
	if (*(unsigned char *)&ident[2] == '0') {
#line 919 "BB02AD.f"
	    if (*(unsigned char *)&ident[1] == '0') {
/*         .. C is not identity matrix .. */
#line 921 "BB02AD.f"
		i__1 = ipar[1];
#line 921 "BB02AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 922 "BB02AD.f"
		    dspmv_("L", &ipar[3], &c_b34, &q[1], &c__[i__ * c_dim1 + 
			    1], &c__1, &c_b7, &dwork[(i__ - 1) * ipar[1] + 1],
			     &c__1, (ftnlen)1);
#line 924 "BB02AD.f"
/* L140: */
#line 924 "BB02AD.f"
		}
/*         .. use Q(1:IPAR(1)) as workspace and compute the first column */
/*            of Q at the end .. */
#line 927 "BB02AD.f"
		isymm = ipar[1] + 1;
#line 928 "BB02AD.f"
		i__1 = ipar[1];
#line 928 "BB02AD.f"
		for (i__ = 2; i__ <= i__1; ++i__) {
#line 929 "BB02AD.f"
		    dgemv_("T", &ipar[3], &ipar[1], &c_b34, &dwork[1], &ipar[
			    1], &c__[i__ * c_dim1 + 1], &c__1, &c_b7, &q[1], &
			    c__1, (ftnlen)1);
#line 931 "BB02AD.f"
		    i__2 = ipar[1] - i__ + 1;
#line 931 "BB02AD.f"
		    dcopy_(&i__2, &q[i__], &c__1, &q[isymm], &c__1);
#line 932 "BB02AD.f"
		    isymm += ipar[1] - i__ + 1;
#line 933 "BB02AD.f"
/* L150: */
#line 933 "BB02AD.f"
		}
#line 934 "BB02AD.f"
		dgemv_("T", &ipar[3], &ipar[1], &c_b34, &dwork[1], &ipar[1], &
			c__[c_dim1 + 1], &c__1, &c_b7, &q[1], &c__1, (ftnlen)
			1);
#line 936 "BB02AD.f"
	    }
#line 937 "BB02AD.f"
	} else {
/*       .. Q = identity .. */
#line 939 "BB02AD.f"
	    if (*(unsigned char *)&ident[1] == '0') {
/*         .. C is not identity matrix .. */
#line 941 "BB02AD.f"
		if (ipar[3] == 1) {
#line 942 "BB02AD.f"
		    dlaset_("L", &nsymm, &c__1, &c_b7, &c_b7, &q[1], &c__1, (
			    ftnlen)1);
#line 943 "BB02AD.f"
		    dspr_("L", &ipar[1], &c_b34, &c__[c_offset], ldc, &q[1], (
			    ftnlen)1);
#line 944 "BB02AD.f"
		} else {
#line 945 "BB02AD.f"
		    dsyrk_("L", "T", &ipar[1], &ipar[3], &c_b34, &c__[
			    c_offset], ldc, &c_b7, &dwork[1], &ipar[1], (
			    ftnlen)1, (ftnlen)1);
#line 947 "BB02AD.f"
		    ma02dd_("Pack", "Lower", &ipar[1], &dwork[1], &ipar[1], &
			    q[1], (ftnlen)4, (ftnlen)5);
#line 948 "BB02AD.f"
		}
#line 949 "BB02AD.f"
	    } else {
/*         .. C = Q = identity .. */
#line 951 "BB02AD.f"
		isymm = 1;
#line 952 "BB02AD.f"
		for (i__ = ipar[1]; i__ >= 1; --i__) {
#line 953 "BB02AD.f"
		    q[isymm] = 1.;
#line 954 "BB02AD.f"
		    isymm += i__;
#line 955 "BB02AD.f"
/* L160: */
#line 955 "BB02AD.f"
		}
#line 956 "BB02AD.f"
	    }
#line 957 "BB02AD.f"
	}
#line 958 "BB02AD.f"
    } else {
#line 959 "BB02AD.f"
	qdimm = ipar[3];
#line 960 "BB02AD.f"
	if (*(unsigned char *)&ident[1] == '1') {
#line 960 "BB02AD.f"
	    dlaset_("A", &ipar[3], &ipar[1], &c_b7, &c_b34, &c__[c_offset], 
		    ldc, (ftnlen)1);
#line 960 "BB02AD.f"
	}
#line 962 "BB02AD.f"
	if (*(unsigned char *)&ident[2] == '1') {
#line 963 "BB02AD.f"
	    isymm = 1;
#line 964 "BB02AD.f"
	    for (i__ = ipar[3]; i__ >= 1; --i__) {
#line 965 "BB02AD.f"
		q[isymm] = 1.;
#line 966 "BB02AD.f"
		isymm += i__;
#line 967 "BB02AD.f"
/* L170: */
#line 967 "BB02AD.f"
	    }
#line 968 "BB02AD.f"
	}
#line 969 "BB02AD.f"
    }

/*     .. unpack symmetric matrices if required .. */
#line 972 "BB02AD.f"
    if (bpar[2]) {
#line 973 "BB02AD.f"
	isymm = qdimm * (qdimm + 1) / 2;
#line 974 "BB02AD.f"
	dcopy_(&isymm, &q[1], &c__1, &dwork[1], &c__1);
#line 975 "BB02AD.f"
	ma02dd_("Unpack", "Lower", &qdimm, &q[1], ldq, &dwork[1], (ftnlen)6, (
		ftnlen)5);
#line 976 "BB02AD.f"
	ma02ed_("Lower", &qdimm, &q[1], ldq, (ftnlen)5);
#line 977 "BB02AD.f"
    } else if (bpar[3]) {
#line 978 "BB02AD.f"
	ma02dd_("Unpack", "Lower", &qdimm, &dwork[1], &qdimm, &q[1], (ftnlen)
		6, (ftnlen)5);
#line 979 "BB02AD.f"
	ma02ed_("Lower", &qdimm, &dwork[1], &qdimm, (ftnlen)5);
#line 980 "BB02AD.f"
	ma02dd_("Pack", "Upper", &qdimm, &dwork[1], &qdimm, &q[1], (ftnlen)4, 
		(ftnlen)5);
#line 981 "BB02AD.f"
    }
#line 982 "BB02AD.f"
    if (bpar[5]) {
#line 983 "BB02AD.f"
	isymm = rdimm * (rdimm + 1) / 2;
#line 984 "BB02AD.f"
	dcopy_(&isymm, &r__[1], &c__1, &dwork[1], &c__1);
#line 985 "BB02AD.f"
	ma02dd_("Unpack", "Lower", &rdimm, &r__[1], ldr, &dwork[1], (ftnlen)6,
		 (ftnlen)5);
#line 986 "BB02AD.f"
	ma02ed_("Lower", &rdimm, &r__[1], ldr, (ftnlen)5);
#line 987 "BB02AD.f"
    } else if (bpar[6]) {
#line 988 "BB02AD.f"
	ma02dd_("Unpack", "Lower", &rdimm, &dwork[1], &rdimm, &r__[1], (
		ftnlen)6, (ftnlen)5);
#line 989 "BB02AD.f"
	ma02ed_("Lower", &rdimm, &dwork[1], &rdimm, (ftnlen)5);
#line 990 "BB02AD.f"
	ma02dd_("Pack", "Upper", &rdimm, &dwork[1], &rdimm, &r__[1], (ftnlen)
		4, (ftnlen)5);
#line 991 "BB02AD.f"
    }

/*     ...set VEC... */
#line 994 "BB02AD.f"
    vec[1] = TRUE_;
#line 995 "BB02AD.f"
    vec[2] = TRUE_;
#line 996 "BB02AD.f"
    vec[3] = TRUE_;
#line 997 "BB02AD.f"
    vec[4] = TRUE_;
#line 998 "BB02AD.f"
    vec[5] = ! bpar[4];
#line 999 "BB02AD.f"
    vec[6] = ! bpar[1];
#line 1000 "BB02AD.f"
    vec[7] = TRUE_;
#line 1001 "BB02AD.f"
    vec[8] = TRUE_;
#line 1002 "BB02AD.f"
    vec[9] = bpar[7];
#line 1003 "BB02AD.f"
    if (nr[1] == 1 && (nr[2] == 1 || nr[2] == 3 || nr[2] == 4) || nr[1] == 2 
	    && (nr[2] == 1 || nr[2] >= 3) || nr[1] == 4) {
#line 1007 "BB02AD.f"
	vec[10] = TRUE_;
#line 1008 "BB02AD.f"
    }
#line 1009 "BB02AD.f"
    s_copy(chpar, notes + (nr[1] + (nr[2] << 2) - 5) * 255, (ftnlen)255, (
	    ftnlen)255);
#line 1010 "BB02AD.f"
    *n = ipar[1];
#line 1011 "BB02AD.f"
    *m = ipar[2];
#line 1012 "BB02AD.f"
    *p = ipar[3];

#line 1014 "BB02AD.f"
L2001:
#line 1015 "BB02AD.f"
    return 0;
/* *** Last line of BB02AD *** */
} /* bb02ad_ */

#undef notes


