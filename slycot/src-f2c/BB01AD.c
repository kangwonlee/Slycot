#line 1 "BB01AD.f"
/* BB01AD.f -- translated by f2c (version 20100827).
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

#line 1 "BB01AD.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static integer c__1 = 1;
static doublereal c_b30 = 1.;
static doublereal c_b31 = 2.;
static doublereal c_b33 = -1.;
static integer c__5 = 5;

/* Subroutine */ int bb01ad_(char *def, integer *nr, doublereal *dpar, 
	integer *ipar, logical *bpar, char *chpar, logical *vec, integer *n, 
	integer *m, integer *p, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, doublereal *g, integer *
	ldg, doublereal *q, integer *ldq, doublereal *x, integer *ldx, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen def_len, 
	ftnlen chpar_len)
{
    /* Initialized data */

    static integer nex[4] = { 6,9,2,4 };
    static integer ndef[36]	/* was [4][9] */ = { 2,2,20,21,2,2,64,100,4,2,
	    0,30,8,2,0,211,9,2,0,0,30,3,0,0,0,4,0,0,0,4,0,0,0,55 };
    static integer mdef[18]	/* was [2][9] */ = { 1,1,1,2,2,1,2,2,3,1,3,3,
	    0,1,0,1,0,2 };
    static integer pdef[18]	/* was [2][9] */ = { 2,1,2,1,4,2,8,2,9,2,5,3,
	    0,2,0,1,0,10 };
    static doublereal pardef[36]	/* was [4][9] */ = { 0.,1e-6,0.,1.,0.,
	    1e-8,0.,.01,0.,1e6,0.0,4.,0.,1e-7,0.0,0.,0.,0.,0.0,0.0,0.,1e6,0.0,
	    0.0,0.0,1e-6,0.0,0.0,0.0,1e-6,0.0,0.0,0.0,1. };
    static struct {
	char e_1[2550];
	char fill_2[255];
	char e_3[765];
	char fill_4[255];
	char e_5[765];
	char fill_6[510];
	char e_7[510];
	char fill_8[765];
	char e_9[255];
	char fill_10[765];
	char e_11[255];
	char fill_12[765];
	char e_13[255];
	char fill_14[510];
	} equiv_41 = { "Laub 1979, Ex.1                                     "
		"                                                            "
		"                                                            "
		"                                                            "
		"                       Arnold/Laub 1984, Ex.1: (A,B) unstabi"
		"lizable as EPS -> 0                                         "
		"                                                            "
		"                                                            "
		"                                      Laub 1979, Ex.4: strin"
		"g of high speed vehicles                                    "
		"                                                            "
		"                                                            "
		"                                                     Laub 19"
		"79, Ex.6: ill-conditioned Riccati equation                  "
		"                                                            "
		"                                                            "
		"                                                            "
		"        Laub 1979, Ex.2: uncontrollable-unobservable data   "
		"                                                            "
		"                                                            "
		"                                                            "
		"                       Arnold/Laub 1984, Ex.3: control weigh"
		"ting matrix singular as EPS -> 0                            "
		"                                                            "
		"                                                            "
		"                                      Laub 1979, Ex.5: circu"
		"lant matrices                                               "
		"                                                            "
		"                                                            "
		"                                                     Rosen/W"
		"ang 1992: lq control of 1-dimensional heat flow             "
		"                                                            "
		"                                                            "
		"                                                            "
		"        Beale/Shafai 1989: model of L-1011 aircraft         "
		"                                                            "
		"                                                            "
		"                                                            "
		"                       Kenney/Laub/Wette 1989, Ex.2: ARE ill"
		" conditioned for EPS -> oo                                  "
		"                                                            "
		"                                                            "
		"                                      ", {0}, "Hench et al. "
		"1995: coupled springs, dashpots and masses                  "
		"                                                            "
		"                                                            "
		"                                                            "
		"  Bhattacharyya et al. 1983: binary distillation column     "
		"                                                            "
		"                                                            "
		"                                                            "
		"                 Bai/Qian 1994: ill-conditioned Hamiltonian "
		"for EPS -> 0                                                "
		"                                                            "
		"                                                            "
		"                                ", {0}, "Lang/Penzl 1994: ro"
		"tating axle                                                 "
		"                                                            "
		"                                                            "
		"                                                        Patn"
		"aik et al. 1980: tubular ammonia reactor                    "
		"                                                            "
		"                                                            "
		"                                                            "
		"           Laub 1992: H-infinity problem, eigenvalues  +/- E"
		"PS +/- i                                                    "
		"                                                            "
		"                                                            "
		"                          ", {0}, "Davison/Gesing 1978: J-10"
		"0 jet engine                                                "
		"                                                            "
		"                                                            "
		"                                                  Petkov et "
		"al. 1987: increasingly badly scaled Hamiltonian as EPS -> oo"
		"                                                            "
		"                                                            "
		"                                                            "
		"     ", {0}, "Chow/Kokotovic 1976: magnetic tape control sys"
		"tem                                                         "
		"                                                            "
		"                                                            "
		"                             ", {0}, "Arnold/Laub 1984, Ex.2"
		": poor sep. of closed-loop spectrum as EPS -> 0             "
		"                                                            "
		"                                                            "
		"                                                     ", {0}, 
		"IFAC Benchmark Problem #90-06: LQG design for modified Boin"
		"g B-767 at flutter condition                                "
		"                                                            "
		"                                                            "
		"                " };

#define notes ((char *)&equiv_41)


    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, x_dim1, 
	    x_offset, i__1, i__2;
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
    double cos(doublereal);

    /* Local variables */
    static integer i__, j, k, l;
    static doublereal b1, b2, c1, c2;
    static integer ios, pos;
    static doublereal sum;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static doublereal temp;
    extern /* Subroutine */ int dspr_(char *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, ftnlen), ma02dd_(char *, 
	    char *, integer *, doublereal *, integer *, doublereal *, ftnlen, 
	    ftnlen), ma02ed_(char *, integer *, doublereal *, integer *, 
	    ftnlen), dscal_(integer *, doublereal *, doublereal *, integer *);
    static integer gdimm;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static char ident[4];
    static integer qdimm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal ttemp;
    extern /* Subroutine */ int dspmv_(char *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     integer *, ftnlen), dsymm_(char *, char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen, ftnlen), dsyrk_(
	    char *, char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static integer isymm, msymm, nsymm, psymm;
    extern doublereal dlapy2_(doublereal *, doublereal *);
    static doublereal appind;
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen), dpptrf_(char *, integer *, 
	    doublereal *, integer *, ftnlen), dpptri_(char *, integer *, 
	    doublereal *, integer *, ftnlen), dpttrf_(integer *, doublereal *,
	     doublereal *, integer *), dpttrs_(integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___15 = { 1, 1, 1, 0, 0 };
    static cilist io___17 = { 1, 1, 1, 0, 0 };
    static cilist io___19 = { 1, 1, 1, 0, 0 };
    static cilist io___20 = { 1, 1, 1, 0, 0 };
    static cilist io___22 = { 1, 1, 1, 0, 0 };
    static cilist io___24 = { 1, 1, 1, 0, 0 };
    static cilist io___25 = { 1, 1, 1, 0, 0 };
    static cilist io___26 = { 1, 1, 1, 0, 0 };
    static cilist io___27 = { 1, 1, 1, 0, 0 };
    static cilist io___28 = { 1, 1, 1, 0, 0 };
    static cilist io___29 = { 1, 1, 1, 0, 0 };
    static cilist io___30 = { 1, 1, 1, 0, 0 };
    static cilist io___37 = { 1, 1, 1, 0, 0 };



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
/*     continuous-time algebraic Riccati equations (CAREs) of the form */

/*       0 = Q + A'X + XA - XGX */

/*     corresponding to the Hamiltonian matrix */

/*            (  A   G  ) */
/*        H = (       T ). */
/*            (  Q  -A  ) */

/*     A,G,Q,X are real N-by-N matrices, Q and G are symmetric and may */
/*     be given in factored form */

/*                   -1 T                         T */
/*      (I)   G = B R  B  ,           (II)   Q = C W C . */

/*     Here, C is P-by-N, W P-by-P, B N-by-M, and R M-by-M, where W */
/*     and R are symmetric. In linear-quadratic optimal control problems, */
/*     usually W is positive semidefinite and R positive definite.  The */
/*     factorized form can be used if the CARE is solved using the */
/*     deflating subspaces of the extended Hamiltonian pencil */

/*                  (  A   0   B  )       (  I   0   0  ) */
/*                  (       T     )       (             ) */
/*        H - s K = (  Q   A   0  )  -  s (  0  -I   0  ) , */
/*                  (       T     )       (             ) */
/*                  (  0   B   R  )       (  0   0   0  ) */

/*     where I and 0 denote the identity and zero matrix, respectively, */
/*     of appropriate dimensions. */

/*     NOTE: the formulation of the CARE and the related matrix (pencils) */
/*           used here does not include CAREs as they arise in robust */
/*           control (H_infinity optimization). */

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
/*             This array determines the example for which CAREX returns */
/*             data. NR(1) is the group of examples. */
/*             NR(1) = 1 : parameter-free problems of fixed size. */
/*             NR(1) = 2 : parameter-dependent problems of fixed size. */
/*             NR(1) = 3 : parameter-free problems of scalable size. */
/*             NR(1) = 4 : parameter-dependent problems of scalable size. */
/*             NR(2) is the number of the example in group NR(1). */
/*             Let NEXi be the number of examples in group i. Currently, */
/*             NEX1 = 6, NEX2 = 9, NEX3 = 2, NEX4 = 4. */
/*             1 <= NR(1) <= 4; */
/*             1 <= NR(2) <= NEXi , where i = NR(1). */

/*     DPAR    (input/output) DOUBLE PRECISION array, dimension (7) */
/*             Double precision parameter vector. For explanation of the */
/*             parameters see [1]. */
/*             DPAR(1)           : defines the parameters */
/*                                 'delta' for NR(1) = 3, */
/*                                 'q' for NR(1).NR(2) = 4.1, */
/*                                 'a' for NR(1).NR(2) = 4.2, and */
/*                                 'mu' for NR(1).NR(2) = 4.3. */
/*             DPAR(2)           : defines parameters */
/*                                 'r' for NR(1).NR(2) = 4.1, */
/*                                 'b' for NR(1).NR(2) = 4.2, and */
/*                                 'delta' for NR(1).NR(2) = 4.3. */
/*             DPAR(3)           : defines parameters */
/*                                 'c' for NR(1).NR(2) = 4.2 and */
/*                                 'kappa' for NR(1).NR(2) = 4.3. */
/*             DPAR(j), j=4,5,6,7: These arguments are only used to */
/*                                 generate Example 4.2 and define in */
/*                                 consecutive order the intervals */
/*                                 ['beta_1', 'beta_2'], */
/*                                 ['gamma_1', 'gamma_2']. */
/*             NOTE that if DEF = 'D' or 'd', the values of DPAR entries */
/*             on input are ignored and, on output, they are overwritten */
/*             with the default parameters. */

/*     IPAR    (input/output) INTEGER array, dimension (3) */
/*             On input, IPAR(1) determines the actual state dimension, */
/*             i.e., the order of the matrix A as follows, where */
/*             NO = NR(1).NR(2). */
/*             NR(1) = 1 or 2.1-2.8: IPAR(1) is ignored. */
/*             NO = 2.9            : IPAR(1) = 1 generates the CARE for */
/*                                   optimal state feedback (default); */
/*                                   IPAR(1) = 2 generates the Kalman */
/*                                   filter CARE. */
/*             NO = 3.1            : IPAR(1) is the number of vehicles */
/*                                   (parameter 'l' in the description */
/*                                    in [1]). */
/*             NO = 3.2, 4.1 or 4.2: IPAR(1) is the order of the matrix */
/*                                   A. */
/*             NO = 4.3 or 4.4     : IPAR(1) determines the dimension of */
/*                                   the second-order system, i.e., the */
/*                                   order of the stiffness matrix for */
/*                                   Examples 4.3 and 4.4 (parameter 'l' */
/*                                   in the description in [1]). */

/*             The order of the output matrix A is N = 2*IPAR(1) for */
/*             Example 4.3 and N = 2*IPAR(1)-1 for Examples 3.1 and 4.4. */
/*             NOTE that IPAR(1) is overwritten for Examples 1.1-2.8. For */
/*             the other examples, IPAR(1) is overwritten if the default */
/*             parameters are to be used. */
/*             On output, IPAR(1) contains the order of the matrix A. */

/*             On input, IPAR(2) is the number of colums in the matrix B */
/*             in (I) (in control problems, the number of inputs of the */
/*             system). Currently, IPAR(2) is fixed or determined by */
/*             IPAR(1) for all examples and thus is not referenced on */
/*             input. */
/*             On output, IPAR(2) is the number of columns of the */
/*             matrix B from (I). */
/*             NOTE that currently IPAR(2) is overwritten and that */
/*             rank(G) <= IPAR(2). */

/*             On input, IPAR(3) is the number of rows in the matrix C */
/*             in (II) (in control problems, the number of outputs of the */
/*             system). Currently, IPAR(3) is fixed or determined by */
/*             IPAR(1) for all examples and thus is not referenced on */
/*             input. */
/*             On output, IPAR(3) contains the number of rows of the */
/*             matrix C in (II). */
/*             NOTE that currently IPAR(3) is overwritten and that */
/*             rank(Q) <= IPAR(3). */

/*     BPAR    (input) BOOLEAN array, dimension (6) */
/*             This array defines the form of the output of the examples */
/*             and the storage mode of the matrices G and Q. */
/*             BPAR(1) = .TRUE.  : G is returned. */
/*             BPAR(1) = .FALSE. : G is returned in factored form, i.e., */
/*                                 B and R from (I) are returned. */
/*             BPAR(2) = .TRUE.  : The matrix returned in array G (i.e., */
/*                                 G if BPAR(1) = .TRUE. and R if */
/*                                 BPAR(1) = .FALSE.) is stored as full */
/*                                 matrix. */
/*             BPAR(2) = .FALSE. : The matrix returned in array G is */
/*                                 provided in packed storage mode. */
/*             BPAR(3) = .TRUE.  : If BPAR(2) = .FALSE., the matrix */
/*                                 returned in array G is stored in upper */
/*                                 packed mode, i.e., the upper triangle */
/*                                 of a symmetric n-by-n matrix is stored */
/*                                 by columns, e.g., the matrix entry */
/*                                 G(i,j) is stored in the array entry */
/*                                 G(i+j*(j-1)/2) for i <= j. */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(3) = .FALSE. : If BPAR(2) = .FALSE., the matrix */
/*                                 returned in array G is stored in lower */
/*                                 packed mode, i.e., the lower triangle */
/*                                 of a symmetric n-by-n matrix is stored */
/*                                 by columns, e.g., the matrix entry */
/*                                 G(i,j) is stored in the array entry */
/*                                 G(i+(2*n-j)*(j-1)/2) for j <= i. */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(4) = .TRUE.  : Q is returned. */
/*             BPAR(4) = .FALSE. : Q is returned in factored form, i.e., */
/*                                 C and W from (II) are returned. */
/*             BPAR(5) = .TRUE.  : The matrix returned in array Q (i.e., */
/*                                 Q if BPAR(4) = .TRUE. and W if */
/*                                 BPAR(4) = .FALSE.) is stored as full */
/*                                 matrix. */
/*             BPAR(5) = .FALSE. : The matrix returned in array Q is */
/*                                 provided in packed storage mode. */
/*             BPAR(6) = .TRUE.  : If BPAR(5) = .FALSE., the matrix */
/*                                 returned in array Q is stored in upper */
/*                                 packed mode (see above). */
/*                                 Otherwise, this entry is ignored. */
/*             BPAR(6) = .FALSE. : If BPAR(5) = .FALSE., the matrix */
/*                                 returned in array Q is stored in lower */
/*                                 packed mode (see above). */
/*                                 Otherwise, this entry is ignored. */
/*             NOTE that there are no default values for BPAR.  If all */
/*             entries are declared to be .TRUE., then matrices G and Q */
/*             are returned in conventional storage mode, i.e., as */
/*             N-by-N arrays where the array element Z(I,J) contains the */
/*             matrix entry Z_{i,j}. */

/*     CHPAR   (input/output) CHARACTER*255 */
/*             On input, this is the name of a data file supplied by the */
/*             user. */
/*             In the current version, only Example 4.4 allows a */
/*             user-defined data file. This file must contain */
/*             consecutively DOUBLE PRECISION vectors mu, delta, gamma, */
/*             and kappa. The length of these vectors is determined by */
/*             the input value for IPAR(1). */
/*             If on entry, IPAR(1) = L, then mu and delta must each */
/*             contain L DOUBLE PRECISION values, and gamma and kappa */
/*             must each contain L-1 DOUBLE PRECISION values. */
/*             On output, this string contains short information about */
/*             the chosen example. */

/*     VEC     (output) LOGICAL array, dimension (9) */
/*             Flag vector which displays the availability of the output */
/*             data: */
/*             VEC(j), j=1,2,3, refer to N, M, and P, respectively, and */
/*             are always .TRUE. */
/*             VEC(4) refers to A and is always .TRUE. */
/*             VEC(5) is .TRUE. if BPAR(1) = .FALSE., i.e., the factors B */
/*             and R from (I) are returned. */
/*             VEC(6) is .TRUE. if BPAR(4) = .FALSE., i.e., the factors C */
/*             and W from (II) are returned. */
/*             VEC(7) refers to G and is always .TRUE. */
/*             VEC(8) refers to Q and is always .TRUE. */
/*             VEC(9) refers to X and is .TRUE. if the exact solution */
/*             matrix is available. */
/*             NOTE that VEC(i) = .FALSE. for i = 1 to 9 if on exit */
/*             INFO .NE. 0. */

/*     N       (output) INTEGER */
/*             The order of the matrices A, X, G if BPAR(1) = .TRUE., and */
/*             Q if BPAR(4) = .TRUE. */

/*     M       (output) INTEGER */
/*             The number of columns in the matrix B (or the dimension of */
/*             the control input space of the underlying dynamical */
/*             system). */

/*     P       (output) INTEGER */
/*             The number of rows in the matrix C (or the dimension of */
/*             the output space of the underlying dynamical system). */

/*     A       (output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array contains the */
/*             coefficient matrix A of the CARE. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= N. */

/*     B       (output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             If (BPAR(1) = .FALSE.), then the leading N-by-M part of */
/*             this array contains the matrix B of the factored form (I) */
/*             of G. Otherwise, B is used as workspace. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= N. */

/*     C       (output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             If (BPAR(4) = .FALSE.), then the leading P-by-N part of */
/*             this array contains the matrix C of the factored form (II) */
/*             of Q. Otherwise, C is used as workspace. */

/*     LDC     INTEGER */
/*             The leading dimension of array C. */
/*             LDC >= P, where P is the number of rows of the matrix C, */
/*             i.e., the output value of IPAR(3). (For all examples, */
/*             P <= N, where N equals the output value of the argument */
/*             IPAR(1), i.e., LDC >= LDA is always safe.) */

/*     G       (output) DOUBLE PRECISION array, dimension (NG) */
/*             If (BPAR(2) = .TRUE.)  then NG = LDG*N. */
/*             If (BPAR(2) = .FALSE.) then NG = N*(N+1)/2. */
/*             If (BPAR(1) = .TRUE.), then array G contains the */
/*             coefficient matrix G of the CARE. */
/*             If (BPAR(1) = .FALSE.), then array G contains the 'control */
/*             weighting matrix' R of G's factored form as in (I). (For */
/*             all examples, M <= N.) The symmetric matrix contained in */
/*             array G is stored according to BPAR(2) and BPAR(3). */

/*     LDG     INTEGER */
/*             If conventional storage mode is used for G, i.e., */
/*             BPAR(2) = .TRUE., then G is stored like a 2-dimensional */
/*             array with leading dimension LDG. If packed symmetric */
/*             storage mode is used, then LDG is not referenced. */
/*             LDG >= N if BPAR(2) = .TRUE.. */

/*     Q       (output) DOUBLE PRECISION array, dimension (NQ) */
/*             If (BPAR(5) = .TRUE.)  then NQ = LDQ*N. */
/*             If (BPAR(5) = .FALSE.) then NQ = N*(N+1)/2. */
/*             If (BPAR(4) = .TRUE.), then array Q contains the */
/*             coefficient matrix Q of the CARE. */
/*             If (BPAR(4) = .FALSE.), then array Q contains the 'output */
/*             weighting matrix' W of Q's factored form as in (II). */
/*             The symmetric matrix contained in array Q is stored */
/*             according to BPAR(5) and BPAR(6). */

/*     LDQ     INTEGER */
/*             If conventional storage mode is used for Q, i.e., */
/*             BPAR(5) = .TRUE., then Q is stored like a 2-dimensional */
/*             array with leading dimension LDQ. If packed symmetric */
/*             storage mode is used, then LDQ is not referenced. */
/*             LDQ >= N if BPAR(5) = .TRUE.. */

/*     X       (output) DOUBLE PRECISION array, dimension (LDX,IPAR(1)) */
/*             If an exact solution is available (NR = 1.1, 1.2, 2.1, */
/*             2.3-2.6, 3.2), then the leading N-by-N part of this array */
/*             contains the solution matrix X in conventional storage */
/*             mode. Otherwise, X is not referenced. */

/*     LDX     INTEGER */
/*             The leading dimension of array X.  LDX >= 1, and */
/*             LDX >= N if NR = 1.1, 1.2, 2.1, 2.3-2.6, 3.2. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= N*MAX(4,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0 : successful exit; */
/*             < 0 : if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1 : data file could not be opened or had wrong format; */
/*             = 2 : division by zero; */
/*             = 3 : G can not be computed as in (I) due to a singular R */
/*                   matrix. */

/*     REFERENCES */

/*     [1] Abels, J. and Benner, P. */
/*         CAREX - A Collection of Benchmark Examples for Continuous-Time */
/*         Algebraic Riccati Equations (Version 2.0). */
/*         SLICOT Working Note 1999-14, November 1999. Available from */
/*         http://www.win.tue.nl/niconet/NIC2/reports.html. */

/*     This is an updated and extended version of */

/*     [2] Benner, P., Laub, A.J., and Mehrmann, V. */
/*         A Collection of Benchmark Examples for the Numerical Solution */
/*         of Algebraic Riccati Equations I: Continuous-Time Case. */
/*         Technical Report SPC 95_22, Fak. f. Mathematik, */
/*         TU Chemnitz-Zwickau (Germany), October 1995. */

/*     NUMERICAL ASPECTS */

/*     If the original data as taken from the literature is given via */
/*     matrices G and Q, but factored forms are requested as output, then */
/*     these factors are obtained from Cholesky or LDL' decompositions of */
/*     G and Q, i.e., the output data will be corrupted by roundoff */
/*     errors. */

/*     FURTHER COMMENTS */

/*     Some benchmark examples read data from the data files provided */
/*     with the collection. */

/*     CONTRIBUTOR */

/*     Peter Benner (Universitaet Bremen), November 15, 1999. */

/*     For questions concerning the collection or for the submission of */
/*     test examples, please send e-mail to benner@math.uni-bremen.de. */

/*     REVISIONS */

/*     1999, December 23 (V. Sima). */

/*     KEYWORDS */

/*     Algebraic Riccati equation, Hamiltonian matrix. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     . # of examples available , # of examples with fixed size. . */

/*     .. Scalar Arguments .. */

/*     .. Array Arguments .. */

/*     .. Local Scalars .. */

/*     ..Local Arrays .. */

/*     .. External Functions .. */
/*     . BLAS . */
/*     . LAPACK . */

/*     .. External Subroutines .. */
/*     . BLAS . */
/*     . LAPACK . */
/*     . SLICOT . */

/*     .. Intrinsic Functions .. */

/*     .. Data Statements .. */
/*     . default values for dimensions . */
#line 446 "BB01AD.f"
    /* Parameter adjustments */
#line 446 "BB01AD.f"
    --nr;
#line 446 "BB01AD.f"
    --dpar;
#line 446 "BB01AD.f"
    --ipar;
#line 446 "BB01AD.f"
    --bpar;
#line 446 "BB01AD.f"
    --vec;
#line 446 "BB01AD.f"
    a_dim1 = *lda;
#line 446 "BB01AD.f"
    a_offset = 1 + a_dim1;
#line 446 "BB01AD.f"
    a -= a_offset;
#line 446 "BB01AD.f"
    b_dim1 = *ldb;
#line 446 "BB01AD.f"
    b_offset = 1 + b_dim1;
#line 446 "BB01AD.f"
    b -= b_offset;
#line 446 "BB01AD.f"
    c_dim1 = *ldc;
#line 446 "BB01AD.f"
    c_offset = 1 + c_dim1;
#line 446 "BB01AD.f"
    c__ -= c_offset;
#line 446 "BB01AD.f"
    --g;
#line 446 "BB01AD.f"
    --q;
#line 446 "BB01AD.f"
    x_dim1 = *ldx;
#line 446 "BB01AD.f"
    x_offset = 1 + x_dim1;
#line 446 "BB01AD.f"
    x -= x_offset;
#line 446 "BB01AD.f"
    --dwork;
#line 446 "BB01AD.f"

#line 446 "BB01AD.f"
    /* Function Body */
/*     . default values for parameters . */
/*     . comments on examples . */

/*     .. Executable Statements .. */

#line 490 "BB01AD.f"
    *info = 0;
#line 491 "BB01AD.f"
    for (i__ = 1; i__ <= 9; ++i__) {
#line 492 "BB01AD.f"
	vec[i__] = FALSE_;
#line 493 "BB01AD.f"
/* L5: */
#line 493 "BB01AD.f"
    }

#line 495 "BB01AD.f"
    if (nr[1] != 1 && ! (lsame_(def, "N", (ftnlen)1, (ftnlen)1) || lsame_(def,
	     "D", (ftnlen)1, (ftnlen)1))) {
#line 497 "BB01AD.f"
	*info = -1;
#line 498 "BB01AD.f"
    } else if (nr[1] < 1 || nr[2] < 1 || nr[1] > 4 || nr[2] > nex[nr[1] - 1]) 
	    {
#line 500 "BB01AD.f"
	*info = -2;
#line 501 "BB01AD.f"
    } else if (nr[1] > 2) {
#line 502 "BB01AD.f"
	if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
#line 502 "BB01AD.f"
	    ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
#line 502 "BB01AD.f"
	}
#line 503 "BB01AD.f"
	if (nr[1] == 3) {
#line 504 "BB01AD.f"
	    if (nr[2] == 1) {
#line 505 "BB01AD.f"
		ipar[2] = ipar[1];
#line 506 "BB01AD.f"
		ipar[3] = ipar[1] - 1;
#line 507 "BB01AD.f"
		ipar[1] = (ipar[1] << 1) - 1;
#line 508 "BB01AD.f"
	    } else if (nr[2] == 2) {
#line 509 "BB01AD.f"
		ipar[2] = ipar[1];
#line 510 "BB01AD.f"
		ipar[3] = ipar[1];
#line 511 "BB01AD.f"
	    } else {
#line 512 "BB01AD.f"
		ipar[2] = 1;
#line 513 "BB01AD.f"
		ipar[3] = 1;
#line 514 "BB01AD.f"
	    }
#line 515 "BB01AD.f"
	} else if (nr[1] == 4) {
#line 516 "BB01AD.f"
	    if (nr[2] == 3) {
#line 517 "BB01AD.f"
		l = ipar[1];
#line 518 "BB01AD.f"
		ipar[2] = 2;
#line 519 "BB01AD.f"
		ipar[3] = l << 1;
#line 520 "BB01AD.f"
		ipar[1] = l << 1;
#line 521 "BB01AD.f"
	    } else if (nr[2] == 4) {
#line 522 "BB01AD.f"
		l = ipar[1];
#line 523 "BB01AD.f"
		ipar[2] = l;
#line 524 "BB01AD.f"
		ipar[3] = l;
#line 525 "BB01AD.f"
		ipar[1] = (l << 1) - 1;
#line 526 "BB01AD.f"
	    } else {
#line 527 "BB01AD.f"
		ipar[2] = 1;
#line 528 "BB01AD.f"
		ipar[3] = 1;
#line 529 "BB01AD.f"
	    }
#line 530 "BB01AD.f"
	}
#line 531 "BB01AD.f"
    } else if (nr[1] == 2 && nr[2] == 9 && ipar[1] == 2) {
#line 533 "BB01AD.f"
	ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
#line 534 "BB01AD.f"
	ipar[2] = mdef[nr[1] + (nr[2] << 1) - 3];
#line 535 "BB01AD.f"
	ipar[3] = 3;
#line 536 "BB01AD.f"
    } else {
#line 537 "BB01AD.f"
	ipar[1] = ndef[nr[1] + (nr[2] << 2) - 5];
#line 538 "BB01AD.f"
	ipar[2] = mdef[nr[1] + (nr[2] << 1) - 3];
#line 539 "BB01AD.f"
	ipar[3] = pdef[nr[1] + (nr[2] << 1) - 3];
#line 540 "BB01AD.f"
    }
#line 541 "BB01AD.f"
    if (*info != 0) {
#line 541 "BB01AD.f"
	goto L7;
#line 541 "BB01AD.f"
    }

#line 543 "BB01AD.f"
    if (ipar[1] < 1) {
#line 544 "BB01AD.f"
	*info = -4;
#line 545 "BB01AD.f"
    } else if (ipar[1] > *lda) {
#line 546 "BB01AD.f"
	*info = -12;
#line 547 "BB01AD.f"
    } else if (ipar[1] > *ldb) {
#line 548 "BB01AD.f"
	*info = -14;
#line 549 "BB01AD.f"
    } else if (ipar[3] > *ldc) {
#line 550 "BB01AD.f"
	*info = -16;
#line 551 "BB01AD.f"
    } else if (bpar[2] && ipar[1] > *ldg) {
#line 552 "BB01AD.f"
	*info = -18;
#line 553 "BB01AD.f"
    } else if (bpar[5] && ipar[1] > *ldq) {
#line 554 "BB01AD.f"
	*info = -20;
#line 555 "BB01AD.f"
    } else if (*ldx < 1) {
#line 556 "BB01AD.f"
	*info = -22;
#line 557 "BB01AD.f"
    } else if (nr[1] == 1 && (nr[2] == 1 || nr[2] == 2)) {
#line 559 "BB01AD.f"
	if (ipar[1] > *ldx) {
#line 559 "BB01AD.f"
	    *info = -22;
#line 559 "BB01AD.f"
	}
#line 560 "BB01AD.f"
    } else if (nr[1] == 2 && nr[2] == 1) {
#line 561 "BB01AD.f"
	if (ipar[1] > *ldx) {
#line 561 "BB01AD.f"
	    *info = -22;
#line 561 "BB01AD.f"
	}
#line 562 "BB01AD.f"
    } else if (nr[1] == 2 && (nr[2] >= 3 && nr[2] <= 6)) {
#line 564 "BB01AD.f"
	if (ipar[1] > *ldx) {
#line 564 "BB01AD.f"
	    *info = -22;
#line 564 "BB01AD.f"
	}
#line 565 "BB01AD.f"
    } else if (nr[1] == 3 && nr[2] == 2) {
#line 566 "BB01AD.f"
	if (ipar[1] > *ldx) {
#line 566 "BB01AD.f"
	    *info = -22;
#line 566 "BB01AD.f"
	}
#line 567 "BB01AD.f"
    } else if (*ldwork < *n * max(4,*n)) {
#line 568 "BB01AD.f"
	*info = -24;
#line 569 "BB01AD.f"
    }

#line 571 "BB01AD.f"
L7:
#line 572 "BB01AD.f"
    if (*info != 0) {
#line 573 "BB01AD.f"
	i__1 = -(*info);
#line 573 "BB01AD.f"
	xerbla_("BB01AD", &i__1, (ftnlen)6);
#line 574 "BB01AD.f"
	return 0;
#line 575 "BB01AD.f"
    }

#line 577 "BB01AD.f"
    nsymm = ipar[1] * (ipar[1] + 1) / 2;
#line 578 "BB01AD.f"
    msymm = ipar[2] * (ipar[2] + 1) / 2;
#line 579 "BB01AD.f"
    psymm = ipar[3] * (ipar[3] + 1) / 2;
#line 580 "BB01AD.f"
    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
#line 580 "BB01AD.f"
	dpar[1] = pardef[nr[1] + (nr[2] << 2) - 5];
#line 580 "BB01AD.f"
    }

#line 582 "BB01AD.f"
    dlaset_("A", &ipar[1], &ipar[1], &c_b10, &c_b10, &a[a_offset], lda, (
	    ftnlen)1);
#line 583 "BB01AD.f"
    dlaset_("A", &ipar[1], &ipar[2], &c_b10, &c_b10, &b[b_offset], ldb, (
	    ftnlen)1);
#line 584 "BB01AD.f"
    dlaset_("A", &ipar[3], &ipar[1], &c_b10, &c_b10, &c__[c_offset], ldc, (
	    ftnlen)1);
#line 585 "BB01AD.f"
    dlaset_("L", &msymm, &c__1, &c_b10, &c_b10, &g[1], &c__1, (ftnlen)1);
#line 586 "BB01AD.f"
    dlaset_("L", &psymm, &c__1, &c_b10, &c_b10, &q[1], &c__1, (ftnlen)1);

#line 588 "BB01AD.f"
    if (nr[1] == 1) {
#line 589 "BB01AD.f"
	if (nr[2] == 1) {
#line 590 "BB01AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 591 "BB01AD.f"
	    b[b_dim1 + 2] = 1.;
#line 592 "BB01AD.f"
	    q[1] = 1.;
#line 593 "BB01AD.f"
	    q[3] = 2.;
#line 594 "BB01AD.f"
	    s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
#line 595 "BB01AD.f"
	    dlaset_("A", &ipar[1], &ipar[1], &c_b30, &c_b31, &x[x_offset], 
		    ldx, (ftnlen)1);

#line 597 "BB01AD.f"
	} else if (nr[2] == 2) {
#line 598 "BB01AD.f"
	    a[a_dim1 + 1] = 4.;
#line 599 "BB01AD.f"
	    a[a_dim1 + 2] = -4.5;
#line 600 "BB01AD.f"
	    a[(a_dim1 << 1) + 1] = 3.;
#line 601 "BB01AD.f"
	    a[(a_dim1 << 1) + 2] = -3.5;
#line 602 "BB01AD.f"
	    dlaset_("A", &ipar[1], &ipar[2], &c_b33, &c_b30, &b[b_offset], 
		    ldb, (ftnlen)1);
#line 603 "BB01AD.f"
	    q[1] = 9.;
#line 604 "BB01AD.f"
	    q[2] = 6.;
#line 605 "BB01AD.f"
	    q[3] = 4.;
#line 606 "BB01AD.f"
	    s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
#line 607 "BB01AD.f"
	    temp = sqrt(2.) + 1.;
#line 608 "BB01AD.f"
	    d__1 = temp * 6.;
#line 608 "BB01AD.f"
	    d__2 = temp * 4.;
#line 608 "BB01AD.f"
	    dlaset_("A", &ipar[1], &ipar[1], &d__1, &d__2, &x[x_offset], ldx, 
		    (ftnlen)1);
#line 610 "BB01AD.f"
	    x[x_dim1 + 1] = temp * 9.;

#line 612 "BB01AD.f"
	} else if (nr[2] >= 3 && nr[2] <= 6) {
#line 613 "BB01AD.f"
	    ici__1.icierr = 0;
#line 613 "BB01AD.f"
	    ici__1.icirnum = 1;
#line 613 "BB01AD.f"
	    ici__1.icirlen = 11;
#line 613 "BB01AD.f"
	    ici__1.iciunit = chpar;
#line 613 "BB01AD.f"
	    ici__1.icifmt = "(A,I1,A,I1,A)";
#line 613 "BB01AD.f"
	    s_wsfi(&ici__1);
#line 613 "BB01AD.f"
	    do_fio(&c__1, "BB01", (ftnlen)4);
#line 613 "BB01AD.f"
	    do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
#line 613 "BB01AD.f"
	    do_fio(&c__1, "0", (ftnlen)1);
#line 613 "BB01AD.f"
	    do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
#line 613 "BB01AD.f"
	    do_fio(&c__1, ".dat", (ftnlen)4);
#line 613 "BB01AD.f"
	    e_wsfi();
#line 615 "BB01AD.f"
	    if (nr[2] == 3 || nr[2] == 4) {
#line 616 "BB01AD.f"
		s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
#line 617 "BB01AD.f"
	    } else if (nr[2] == 5) {
#line 618 "BB01AD.f"
		s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);
#line 619 "BB01AD.f"
	    } else if (nr[2] == 6) {
#line 620 "BB01AD.f"
		s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);
#line 621 "BB01AD.f"
	    }
#line 622 "BB01AD.f"
	    o__1.oerr = 1;
#line 622 "BB01AD.f"
	    o__1.ounit = 1;
#line 622 "BB01AD.f"
	    o__1.ofnmlen = 11;
#line 622 "BB01AD.f"
	    o__1.ofnm = chpar;
#line 622 "BB01AD.f"
	    o__1.orl = 0;
#line 622 "BB01AD.f"
	    o__1.osta = "OLD";
#line 622 "BB01AD.f"
	    o__1.oacc = 0;
#line 622 "BB01AD.f"
	    o__1.ofm = 0;
#line 622 "BB01AD.f"
	    o__1.oblnk = 0;
#line 622 "BB01AD.f"
	    ios = f_open(&o__1);
#line 623 "BB01AD.f"
	    if (ios != 0) {
#line 624 "BB01AD.f"
		*info = 1;
#line 625 "BB01AD.f"
	    } else if (nr[2] <= 6) {
#line 626 "BB01AD.f"
		i__1 = ipar[1];
#line 626 "BB01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 627 "BB01AD.f"
		    ios = s_rsle(&io___15);
#line 627 "BB01AD.f"
		    if (ios != 0) {
#line 627 "BB01AD.f"
			goto L100001;
#line 627 "BB01AD.f"
		    }
#line 627 "BB01AD.f"
		    i__2 = ipar[1];
#line 627 "BB01AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 627 "BB01AD.f"
			ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				a_dim1], (ftnlen)sizeof(doublereal));
#line 627 "BB01AD.f"
			if (ios != 0) {
#line 627 "BB01AD.f"
			    goto L100001;
#line 627 "BB01AD.f"
			}
#line 627 "BB01AD.f"
		    }
#line 627 "BB01AD.f"
		    ios = e_rsle();
#line 627 "BB01AD.f"
L100001:
#line 629 "BB01AD.f"
		    if (ios != 0) {
#line 629 "BB01AD.f"
			*info = 1;
#line 629 "BB01AD.f"
		    }
#line 630 "BB01AD.f"
/* L10: */
#line 630 "BB01AD.f"
		}
#line 631 "BB01AD.f"
		i__1 = ipar[1];
#line 631 "BB01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 632 "BB01AD.f"
		    ios = s_rsle(&io___17);
#line 632 "BB01AD.f"
		    if (ios != 0) {
#line 632 "BB01AD.f"
			goto L100002;
#line 632 "BB01AD.f"
		    }
#line 632 "BB01AD.f"
		    i__2 = ipar[2];
#line 632 "BB01AD.f"
		    for (j = 1; j <= i__2; ++j) {
#line 632 "BB01AD.f"
			ios = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				b_dim1], (ftnlen)sizeof(doublereal));
#line 632 "BB01AD.f"
			if (ios != 0) {
#line 632 "BB01AD.f"
			    goto L100002;
#line 632 "BB01AD.f"
			}
#line 632 "BB01AD.f"
		    }
#line 632 "BB01AD.f"
		    ios = e_rsle();
#line 632 "BB01AD.f"
L100002:
#line 634 "BB01AD.f"
		    if (ios != 0) {
#line 634 "BB01AD.f"
			*info = 1;
#line 634 "BB01AD.f"
		    }
#line 635 "BB01AD.f"
/* L20: */
#line 635 "BB01AD.f"
		}
#line 636 "BB01AD.f"
		if (nr[2] <= 4) {
#line 637 "BB01AD.f"
		    i__1 = ipar[1];
#line 637 "BB01AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 638 "BB01AD.f"
			pos = (i__ - 1) * ipar[1];
#line 639 "BB01AD.f"
			ios = s_rsle(&io___19);
#line 639 "BB01AD.f"
			if (ios != 0) {
#line 639 "BB01AD.f"
			    goto L100003;
#line 639 "BB01AD.f"
			}
#line 639 "BB01AD.f"
			i__2 = ipar[1];
#line 639 "BB01AD.f"
			for (j = 1; j <= i__2; ++j) {
#line 639 "BB01AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&dwork[pos + j]
				    , (ftnlen)sizeof(doublereal));
#line 639 "BB01AD.f"
			    if (ios != 0) {
#line 639 "BB01AD.f"
				goto L100003;
#line 639 "BB01AD.f"
			    }
#line 639 "BB01AD.f"
			}
#line 639 "BB01AD.f"
			ios = e_rsle();
#line 639 "BB01AD.f"
L100003:
#line 641 "BB01AD.f"
/* L30: */
#line 641 "BB01AD.f"
			;
#line 641 "BB01AD.f"
		    }
#line 642 "BB01AD.f"
		    if (ios != 0) {
#line 643 "BB01AD.f"
			*info = 1;
#line 644 "BB01AD.f"
		    } else {
#line 645 "BB01AD.f"
			ma02dd_("Pack", "Lower", &ipar[1], &dwork[1], &ipar[1]
				, &q[1], (ftnlen)4, (ftnlen)5);
#line 646 "BB01AD.f"
		    }
#line 647 "BB01AD.f"
		} else if (nr[2] == 6) {
#line 648 "BB01AD.f"
		    i__1 = ipar[3];
#line 648 "BB01AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 649 "BB01AD.f"
			ios = s_rsle(&io___20);
#line 649 "BB01AD.f"
			if (ios != 0) {
#line 649 "BB01AD.f"
			    goto L100004;
#line 649 "BB01AD.f"
			}
#line 649 "BB01AD.f"
			i__2 = ipar[1];
#line 649 "BB01AD.f"
			for (j = 1; j <= i__2; ++j) {
#line 649 "BB01AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&c__[i__ + j * 
				    c_dim1], (ftnlen)sizeof(doublereal));
#line 649 "BB01AD.f"
			    if (ios != 0) {
#line 649 "BB01AD.f"
				goto L100004;
#line 649 "BB01AD.f"
			    }
#line 649 "BB01AD.f"
			}
#line 649 "BB01AD.f"
			ios = e_rsle();
#line 649 "BB01AD.f"
L100004:
#line 651 "BB01AD.f"
			if (ios != 0) {
#line 651 "BB01AD.f"
			    *info = 1;
#line 651 "BB01AD.f"
			}
#line 652 "BB01AD.f"
/* L35: */
#line 652 "BB01AD.f"
		    }
#line 653 "BB01AD.f"
		}
#line 654 "BB01AD.f"
		cl__1.cerr = 0;
#line 654 "BB01AD.f"
		cl__1.cunit = 1;
#line 654 "BB01AD.f"
		cl__1.csta = 0;
#line 654 "BB01AD.f"
		f_clos(&cl__1);
#line 655 "BB01AD.f"
	    }
#line 656 "BB01AD.f"
	}

#line 658 "BB01AD.f"
    } else if (nr[1] == 2) {
#line 659 "BB01AD.f"
	if (nr[2] == 1) {
#line 660 "BB01AD.f"
	    a[a_dim1 + 1] = 1.;
#line 661 "BB01AD.f"
	    a[(a_dim1 << 1) + 2] = -2.;
#line 662 "BB01AD.f"
	    b[b_dim1 + 1] = dpar[1];
#line 663 "BB01AD.f"
	    dlaset_("U", &ipar[3], &ipar[1], &c_b30, &c_b30, &c__[c_offset], 
		    ldc, (ftnlen)1);
#line 664 "BB01AD.f"
	    s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);
#line 665 "BB01AD.f"
	    if (dpar[1] != 0.) {
#line 666 "BB01AD.f"
		temp = dlapy2_(&c_b30, &dpar[1]);
#line 667 "BB01AD.f"
		x[x_dim1 + 1] = (temp + 1.) / dpar[1] / dpar[1];
#line 668 "BB01AD.f"
		x[x_dim1 + 2] = 1. / (temp + 2.);
#line 669 "BB01AD.f"
		x[(x_dim1 << 1) + 1] = x[x_dim1 + 2];
#line 670 "BB01AD.f"
		ttemp = dpar[1] * x[(x_dim1 << 1) + 1];
#line 671 "BB01AD.f"
		temp = (1. - ttemp) * (ttemp + 1.);
#line 672 "BB01AD.f"
		x[(x_dim1 << 1) + 2] = temp / 4.;
#line 673 "BB01AD.f"
	    } else {
#line 674 "BB01AD.f"
		*info = 2;
#line 675 "BB01AD.f"
	    }

#line 677 "BB01AD.f"
	} else if (nr[2] == 2) {
#line 678 "BB01AD.f"
	    a[a_dim1 + 1] = -.1;
#line 679 "BB01AD.f"
	    a[(a_dim1 << 1) + 2] = -.02;
#line 680 "BB01AD.f"
	    b[b_dim1 + 1] = .1;
#line 681 "BB01AD.f"
	    b[b_dim1 + 2] = .001;
#line 682 "BB01AD.f"
	    b[(b_dim1 << 1) + 2] = .01;
#line 683 "BB01AD.f"
	    dlaset_("L", &msymm, &c__1, &c_b30, &c_b30, &g[1], &msymm, (
		    ftnlen)1);
#line 684 "BB01AD.f"
	    g[1] += dpar[1];
#line 685 "BB01AD.f"
	    c__[c_dim1 + 1] = 10.;
#line 686 "BB01AD.f"
	    c__[(c_dim1 << 1) + 1] = 100.;
#line 687 "BB01AD.f"
	    s_copy(ident, "0010", (ftnlen)4, (ftnlen)4);

#line 689 "BB01AD.f"
	} else if (nr[2] == 3) {
#line 690 "BB01AD.f"
	    a[(a_dim1 << 1) + 1] = dpar[1];
#line 691 "BB01AD.f"
	    b[b_dim1 + 2] = 1.;
#line 692 "BB01AD.f"
	    s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);
#line 693 "BB01AD.f"
	    if (dpar[1] != 0.) {
#line 694 "BB01AD.f"
		temp = sqrt(dpar[1] * 2. + 1.);
#line 695 "BB01AD.f"
		dlaset_("A", &ipar[1], &ipar[1], &c_b30, &temp, &x[x_offset], 
			ldx, (ftnlen)1);
#line 696 "BB01AD.f"
		x[x_dim1 + 1] /= dpar[1];
#line 697 "BB01AD.f"
	    } else {
#line 698 "BB01AD.f"
		*info = 2;
#line 699 "BB01AD.f"
	    }

#line 701 "BB01AD.f"
	} else if (nr[2] == 4) {
#line 702 "BB01AD.f"
	    temp = dpar[1] + 1.;
#line 703 "BB01AD.f"
	    dlaset_("A", &ipar[1], &ipar[1], &c_b30, &temp, &a[a_offset], lda,
		     (ftnlen)1);
/* Computing 2nd power */
#line 704 "BB01AD.f"
	    d__1 = dpar[1];
#line 704 "BB01AD.f"
	    q[1] = d__1 * d__1;
#line 705 "BB01AD.f"
	    q[3] = q[1];
#line 706 "BB01AD.f"
	    s_copy(ident, "1101", (ftnlen)4, (ftnlen)4);
/* Computing 2nd power */
#line 707 "BB01AD.f"
	    d__1 = temp;
#line 707 "BB01AD.f"
	    x[x_dim1 + 1] = temp * 2. + sqrt(2.) * (sqrt(d__1 * d__1 + 1.) + 
		    dpar[1]);
#line 708 "BB01AD.f"
	    x[x_dim1 + 1] /= 2.;
#line 709 "BB01AD.f"
	    x[(x_dim1 << 1) + 2] = x[x_dim1 + 1];
#line 710 "BB01AD.f"
	    ttemp = x[x_dim1 + 1] - temp;
#line 711 "BB01AD.f"
	    if (ttemp != 0.) {
#line 712 "BB01AD.f"
		x[x_dim1 + 2] = x[x_dim1 + 1] / ttemp;
#line 713 "BB01AD.f"
		x[(x_dim1 << 1) + 1] = x[x_dim1 + 2];
#line 714 "BB01AD.f"
	    } else {
#line 715 "BB01AD.f"
		*info = 2;
#line 716 "BB01AD.f"
	    }

#line 718 "BB01AD.f"
	} else if (nr[2] == 5) {
#line 719 "BB01AD.f"
	    a[a_dim1 + 1] = 3. - dpar[1];
#line 720 "BB01AD.f"
	    a[a_dim1 + 2] = 4.;
#line 721 "BB01AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 722 "BB01AD.f"
	    a[(a_dim1 << 1) + 2] = 2. - dpar[1];
#line 723 "BB01AD.f"
	    dlaset_("L", &ipar[1], &ipar[2], &c_b30, &c_b30, &b[b_offset], 
		    ldb, (ftnlen)1);
#line 724 "BB01AD.f"
	    q[1] = dpar[1] * 4. - 11.;
#line 725 "BB01AD.f"
	    q[2] = dpar[1] * 2. - 5.;
#line 726 "BB01AD.f"
	    q[3] = dpar[1] * 2. - 2.;
#line 727 "BB01AD.f"
	    s_copy(ident, "0101", (ftnlen)4, (ftnlen)4);
#line 728 "BB01AD.f"
	    dlaset_("A", &ipar[1], &ipar[1], &c_b30, &c_b30, &x[x_offset], 
		    ldx, (ftnlen)1);
#line 729 "BB01AD.f"
	    x[x_dim1 + 1] = 2.;

#line 731 "BB01AD.f"
	} else if (nr[2] == 6) {
#line 732 "BB01AD.f"
	    if (dpar[1] != 0.) {
#line 733 "BB01AD.f"
		a[a_dim1 + 1] = dpar[1];
#line 734 "BB01AD.f"
		a[(a_dim1 << 1) + 2] = dpar[1] * 2.;
#line 735 "BB01AD.f"
		a[a_dim1 * 3 + 3] = dpar[1] * 3.;
/*     .. set C = V .. */
#line 737 "BB01AD.f"
		temp = .66666666666666663;
#line 738 "BB01AD.f"
		d__1 = -temp;
#line 738 "BB01AD.f"
		d__2 = 1. - temp;
#line 738 "BB01AD.f"
		dlaset_("A", &ipar[3], &ipar[1], &d__1, &d__2, &c__[c_offset],
			 ldc, (ftnlen)1);
#line 740 "BB01AD.f"
		dsymm_("L", "L", &ipar[1], &ipar[1], &c_b30, &c__[c_offset], 
			ldc, &a[a_offset], lda, &c_b10, &dwork[1], &ipar[1], (
			ftnlen)1, (ftnlen)1);
#line 742 "BB01AD.f"
		dsymm_("R", "L", &ipar[1], &ipar[1], &c_b30, &c__[c_offset], 
			ldc, &dwork[1], &ipar[1], &c_b10, &a[a_offset], lda, (
			ftnlen)1, (ftnlen)1);
/*     .. G = R ! .. */
#line 745 "BB01AD.f"
		g[1] = dpar[1];
#line 746 "BB01AD.f"
		g[4] = dpar[1];
#line 747 "BB01AD.f"
		g[6] = dpar[1];
#line 748 "BB01AD.f"
		q[1] = 1. / dpar[1];
#line 749 "BB01AD.f"
		q[4] = 1.;
#line 750 "BB01AD.f"
		q[6] = dpar[1];
#line 751 "BB01AD.f"
		s_copy(ident, "1000", (ftnlen)4, (ftnlen)4);
#line 752 "BB01AD.f"
		dlaset_("A", &ipar[1], &ipar[1], &c_b10, &c_b10, &x[x_offset],
			 ldx, (ftnlen)1);
/* Computing 2nd power */
#line 753 "BB01AD.f"
		d__1 = dpar[1];
#line 753 "BB01AD.f"
		temp = d__1 * d__1;
/* Computing 2nd power */
#line 754 "BB01AD.f"
		d__1 = temp;
#line 754 "BB01AD.f"
		x[x_dim1 + 1] = temp + sqrt(d__1 * d__1 + 1.);
/* Computing 2nd power */
#line 755 "BB01AD.f"
		d__1 = temp;
#line 755 "BB01AD.f"
		x[(x_dim1 << 1) + 2] = temp * 2. + sqrt(d__1 * d__1 * 4. + 
			dpar[1]);
#line 756 "BB01AD.f"
		x[x_dim1 * 3 + 3] = temp * 3. + dpar[1] * sqrt(temp * 9. + 1.)
			;
#line 757 "BB01AD.f"
		dsymm_("L", "L", &ipar[1], &ipar[1], &c_b30, &c__[c_offset], 
			ldc, &x[x_offset], ldx, &c_b10, &dwork[1], &ipar[1], (
			ftnlen)1, (ftnlen)1);
#line 759 "BB01AD.f"
		dsymm_("R", "L", &ipar[1], &ipar[1], &c_b30, &c__[c_offset], 
			ldc, &dwork[1], &ipar[1], &c_b10, &x[x_offset], ldx, (
			ftnlen)1, (ftnlen)1);
#line 761 "BB01AD.f"
	    } else {
#line 762 "BB01AD.f"
		*info = 2;
#line 763 "BB01AD.f"
	    }

#line 765 "BB01AD.f"
	} else if (nr[2] == 7) {
#line 766 "BB01AD.f"
	    if (dpar[1] != 0.) {
#line 767 "BB01AD.f"
		a[(a_dim1 << 1) + 1] = .4;
#line 768 "BB01AD.f"
		a[a_dim1 * 3 + 2] = .345;
#line 769 "BB01AD.f"
		a[(a_dim1 << 1) + 3] = -.524 / dpar[1];
#line 770 "BB01AD.f"
		a[a_dim1 * 3 + 3] = -.465 / dpar[1];
#line 771 "BB01AD.f"
		a[(a_dim1 << 2) + 3] = .262 / dpar[1];
#line 772 "BB01AD.f"
		a[(a_dim1 << 2) + 4] = -1. / dpar[1];
#line 773 "BB01AD.f"
		b[b_dim1 + 4] = 1. / dpar[1];
#line 774 "BB01AD.f"
		c__[c_dim1 + 1] = 1.;
#line 775 "BB01AD.f"
		c__[c_dim1 * 3 + 2] = 1.;
#line 776 "BB01AD.f"
		s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);
#line 777 "BB01AD.f"
	    } else {
#line 778 "BB01AD.f"
		*info = 2;
#line 779 "BB01AD.f"
	    }

#line 781 "BB01AD.f"
	} else if (nr[2] == 8) {
#line 782 "BB01AD.f"
	    a[a_dim1 + 1] = -dpar[1];
#line 783 "BB01AD.f"
	    a[a_dim1 + 2] = -1.;
#line 784 "BB01AD.f"
	    a[(a_dim1 << 1) + 1] = 1.;
#line 785 "BB01AD.f"
	    a[(a_dim1 << 1) + 2] = -dpar[1];
#line 786 "BB01AD.f"
	    a[a_dim1 * 3 + 3] = dpar[1];
#line 787 "BB01AD.f"
	    a[a_dim1 * 3 + 4] = -1.;
#line 788 "BB01AD.f"
	    a[(a_dim1 << 2) + 3] = 1.;
#line 789 "BB01AD.f"
	    a[(a_dim1 << 2) + 4] = dpar[1];
#line 790 "BB01AD.f"
	    dlaset_("L", &ipar[1], &ipar[2], &c_b30, &c_b30, &b[b_offset], 
		    ldb, (ftnlen)1);
#line 791 "BB01AD.f"
	    dlaset_("U", &ipar[3], &ipar[1], &c_b30, &c_b30, &c__[c_offset], 
		    ldc, (ftnlen)1);
#line 792 "BB01AD.f"
	    s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);

#line 794 "BB01AD.f"
	} else if (nr[2] == 9) {
#line 795 "BB01AD.f"
	    if (ipar[3] == 10) {
/*     .. read LQR CARE ... */
#line 797 "BB01AD.f"
		ici__1.icierr = 0;
#line 797 "BB01AD.f"
		ici__1.icirnum = 1;
#line 797 "BB01AD.f"
		ici__1.icirlen = 12;
#line 797 "BB01AD.f"
		ici__1.iciunit = chpar;
#line 797 "BB01AD.f"
		ici__1.icifmt = "(A,I1,A,I1,A)";
#line 797 "BB01AD.f"
		s_wsfi(&ici__1);
#line 797 "BB01AD.f"
		do_fio(&c__1, "BB01", (ftnlen)4);
#line 797 "BB01AD.f"
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
#line 797 "BB01AD.f"
		do_fio(&c__1, "0", (ftnlen)1);
#line 797 "BB01AD.f"
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
#line 797 "BB01AD.f"
		do_fio(&c__1, "1.dat", (ftnlen)5);
#line 797 "BB01AD.f"
		e_wsfi();
#line 799 "BB01AD.f"
		o__1.oerr = 1;
#line 799 "BB01AD.f"
		o__1.ounit = 1;
#line 799 "BB01AD.f"
		o__1.ofnmlen = 12;
#line 799 "BB01AD.f"
		o__1.ofnm = chpar;
#line 799 "BB01AD.f"
		o__1.orl = 0;
#line 799 "BB01AD.f"
		o__1.osta = "OLD";
#line 799 "BB01AD.f"
		o__1.oacc = 0;
#line 799 "BB01AD.f"
		o__1.ofm = 0;
#line 799 "BB01AD.f"
		o__1.oblnk = 0;
#line 799 "BB01AD.f"
		ios = f_open(&o__1);
#line 800 "BB01AD.f"
		if (ios != 0) {
#line 801 "BB01AD.f"
		    *info = 1;
#line 802 "BB01AD.f"
		} else {
#line 803 "BB01AD.f"
		    for (i__ = 1; i__ <= 27; i__ += 2) {
#line 804 "BB01AD.f"
			ios = s_rsle(&io___22);
#line 804 "BB01AD.f"
			if (ios != 0) {
#line 804 "BB01AD.f"
			    goto L100005;
#line 804 "BB01AD.f"
			}
#line 804 "BB01AD.f"
			for (j = 0; j <= 1; ++j) {
#line 804 "BB01AD.f"
			    for (k = 0; k <= 1; ++k) {
#line 804 "BB01AD.f"
				ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j 
					+ (i__ + k) * a_dim1], (ftnlen)sizeof(
					doublereal));
#line 804 "BB01AD.f"
				if (ios != 0) {
#line 804 "BB01AD.f"
				    goto L100005;
#line 804 "BB01AD.f"
				}
#line 804 "BB01AD.f"
			    }
#line 804 "BB01AD.f"
			}
#line 804 "BB01AD.f"
			ios = e_rsle();
#line 804 "BB01AD.f"
L100005:
#line 806 "BB01AD.f"
			if (ios != 0) {
#line 806 "BB01AD.f"
			    *info = 1;
#line 806 "BB01AD.f"
			}
#line 807 "BB01AD.f"
/* L36: */
#line 807 "BB01AD.f"
		    }
#line 808 "BB01AD.f"
		    for (i__ = 30; i__ <= 44; i__ += 2) {
#line 809 "BB01AD.f"
			ios = s_rsle(&io___24);
#line 809 "BB01AD.f"
			if (ios != 0) {
#line 809 "BB01AD.f"
			    goto L100006;
#line 809 "BB01AD.f"
			}
#line 809 "BB01AD.f"
			for (j = 0; j <= 1; ++j) {
#line 809 "BB01AD.f"
			    for (k = 0; k <= 1; ++k) {
#line 809 "BB01AD.f"
				ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j 
					+ (i__ + k) * a_dim1], (ftnlen)sizeof(
					doublereal));
#line 809 "BB01AD.f"
				if (ios != 0) {
#line 809 "BB01AD.f"
				    goto L100006;
#line 809 "BB01AD.f"
				}
#line 809 "BB01AD.f"
			    }
#line 809 "BB01AD.f"
			}
#line 809 "BB01AD.f"
			ios = e_rsle();
#line 809 "BB01AD.f"
L100006:
#line 811 "BB01AD.f"
			if (ios != 0) {
#line 811 "BB01AD.f"
			    *info = 1;
#line 811 "BB01AD.f"
			}
#line 812 "BB01AD.f"
/* L37: */
#line 812 "BB01AD.f"
		    }
#line 813 "BB01AD.f"
		    i__1 = ipar[1];
#line 813 "BB01AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 814 "BB01AD.f"
			ios = s_rsle(&io___25);
#line 814 "BB01AD.f"
			if (ios != 0) {
#line 814 "BB01AD.f"
			    goto L100007;
#line 814 "BB01AD.f"
			}
#line 814 "BB01AD.f"
			i__2 = ipar[1];
#line 814 "BB01AD.f"
			for (j = 46; j <= i__2; ++j) {
#line 814 "BB01AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&a[i__ + j * 
				    a_dim1], (ftnlen)sizeof(doublereal));
#line 814 "BB01AD.f"
			    if (ios != 0) {
#line 814 "BB01AD.f"
				goto L100007;
#line 814 "BB01AD.f"
			    }
#line 814 "BB01AD.f"
			}
#line 814 "BB01AD.f"
			ios = e_rsle();
#line 814 "BB01AD.f"
L100007:
#line 816 "BB01AD.f"
			if (ios != 0) {
#line 816 "BB01AD.f"
			    *info = 1;
#line 816 "BB01AD.f"
			}
#line 817 "BB01AD.f"
/* L38: */
#line 817 "BB01AD.f"
		    }
#line 818 "BB01AD.f"
		    a[a_dim1 * 29 + 29] = -5.301;
#line 819 "BB01AD.f"
		    b[b_dim1 + 48] = 8e5;
#line 820 "BB01AD.f"
		    b[(b_dim1 << 1) + 51] = 8e5;
#line 821 "BB01AD.f"
		    g[1] = 364.7;
#line 822 "BB01AD.f"
		    g[3] = 14.59;
#line 823 "BB01AD.f"
		    for (i__ = 1; i__ <= 6; ++i__) {
#line 824 "BB01AD.f"
			ios = s_rsle(&io___26);
#line 824 "BB01AD.f"
			if (ios != 0) {
#line 824 "BB01AD.f"
			    goto L100008;
#line 824 "BB01AD.f"
			}
#line 824 "BB01AD.f"
			for (j = 1; j <= 45; ++j) {
#line 824 "BB01AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&c__[i__ + j * 
				    c_dim1], (ftnlen)sizeof(doublereal));
#line 824 "BB01AD.f"
			    if (ios != 0) {
#line 824 "BB01AD.f"
				goto L100008;
#line 824 "BB01AD.f"
			    }
#line 824 "BB01AD.f"
			}
#line 824 "BB01AD.f"
			ios = e_rsle();
#line 824 "BB01AD.f"
L100008:
#line 826 "BB01AD.f"
			if (ios != 0) {
#line 826 "BB01AD.f"
			    *info = 1;
#line 826 "BB01AD.f"
			}
#line 827 "BB01AD.f"
/* L39: */
#line 827 "BB01AD.f"
		    }
#line 828 "BB01AD.f"
		    c__[c_dim1 * 47 + 7] = 1.;
#line 829 "BB01AD.f"
		    c__[c_dim1 * 46 + 8] = 1.;
#line 830 "BB01AD.f"
		    c__[c_dim1 * 50 + 9] = 1.;
#line 831 "BB01AD.f"
		    c__[c_dim1 * 49 + 10] = 1.;
#line 832 "BB01AD.f"
		    q[11] = 3.76e-14;
#line 833 "BB01AD.f"
		    q[20] = 1.2e-13;
#line 834 "BB01AD.f"
		    q[41] = 2.45e-12;
#line 835 "BB01AD.f"
		}
#line 836 "BB01AD.f"
	    } else {
/*     .. read Kalman filter CARE .. */
#line 838 "BB01AD.f"
		ici__1.icierr = 0;
#line 838 "BB01AD.f"
		ici__1.icirnum = 1;
#line 838 "BB01AD.f"
		ici__1.icirlen = 12;
#line 838 "BB01AD.f"
		ici__1.iciunit = chpar;
#line 838 "BB01AD.f"
		ici__1.icifmt = "(A,I1,A,I1,A)";
#line 838 "BB01AD.f"
		s_wsfi(&ici__1);
#line 838 "BB01AD.f"
		do_fio(&c__1, "BB01", (ftnlen)4);
#line 838 "BB01AD.f"
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
#line 838 "BB01AD.f"
		do_fio(&c__1, "0", (ftnlen)1);
#line 838 "BB01AD.f"
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
#line 838 "BB01AD.f"
		do_fio(&c__1, "2.dat", (ftnlen)5);
#line 838 "BB01AD.f"
		e_wsfi();
#line 840 "BB01AD.f"
		o__1.oerr = 1;
#line 840 "BB01AD.f"
		o__1.ounit = 1;
#line 840 "BB01AD.f"
		o__1.ofnmlen = 12;
#line 840 "BB01AD.f"
		o__1.ofnm = chpar;
#line 840 "BB01AD.f"
		o__1.orl = 0;
#line 840 "BB01AD.f"
		o__1.osta = "OLD";
#line 840 "BB01AD.f"
		o__1.oacc = 0;
#line 840 "BB01AD.f"
		o__1.ofm = 0;
#line 840 "BB01AD.f"
		o__1.oblnk = 0;
#line 840 "BB01AD.f"
		ios = f_open(&o__1);
#line 841 "BB01AD.f"
		if (ios != 0) {
#line 842 "BB01AD.f"
		    *info = 1;
#line 843 "BB01AD.f"
		} else {
#line 844 "BB01AD.f"
		    for (i__ = 1; i__ <= 27; i__ += 2) {
#line 845 "BB01AD.f"
			ios = s_rsle(&io___27);
#line 845 "BB01AD.f"
			if (ios != 0) {
#line 845 "BB01AD.f"
			    goto L100009;
#line 845 "BB01AD.f"
			}
#line 845 "BB01AD.f"
			for (j = 0; j <= 1; ++j) {
#line 845 "BB01AD.f"
			    for (k = 0; k <= 1; ++k) {
#line 845 "BB01AD.f"
				ios = do_lio(&c__5, &c__1, (char *)&a[i__ + k 
					+ (i__ + j) * a_dim1], (ftnlen)sizeof(
					doublereal));
#line 845 "BB01AD.f"
				if (ios != 0) {
#line 845 "BB01AD.f"
				    goto L100009;
#line 845 "BB01AD.f"
				}
#line 845 "BB01AD.f"
			    }
#line 845 "BB01AD.f"
			}
#line 845 "BB01AD.f"
			ios = e_rsle();
#line 845 "BB01AD.f"
L100009:
#line 847 "BB01AD.f"
			if (ios != 0) {
#line 847 "BB01AD.f"
			    *info = 1;
#line 847 "BB01AD.f"
			}
#line 848 "BB01AD.f"
/* L40: */
#line 848 "BB01AD.f"
		    }
#line 849 "BB01AD.f"
		    for (i__ = 30; i__ <= 44; i__ += 2) {
#line 850 "BB01AD.f"
			ios = s_rsle(&io___28);
#line 850 "BB01AD.f"
			if (ios != 0) {
#line 850 "BB01AD.f"
			    goto L100010;
#line 850 "BB01AD.f"
			}
#line 850 "BB01AD.f"
			for (j = 0; j <= 1; ++j) {
#line 850 "BB01AD.f"
			    for (k = 0; k <= 1; ++k) {
#line 850 "BB01AD.f"
				ios = do_lio(&c__5, &c__1, (char *)&a[i__ + k 
					+ (i__ + j) * a_dim1], (ftnlen)sizeof(
					doublereal));
#line 850 "BB01AD.f"
				if (ios != 0) {
#line 850 "BB01AD.f"
				    goto L100010;
#line 850 "BB01AD.f"
				}
#line 850 "BB01AD.f"
			    }
#line 850 "BB01AD.f"
			}
#line 850 "BB01AD.f"
			ios = e_rsle();
#line 850 "BB01AD.f"
L100010:
#line 852 "BB01AD.f"
			if (ios != 0) {
#line 852 "BB01AD.f"
			    *info = 1;
#line 852 "BB01AD.f"
			}
#line 853 "BB01AD.f"
/* L41: */
#line 853 "BB01AD.f"
		    }
#line 854 "BB01AD.f"
		    i__1 = ipar[1];
#line 854 "BB01AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 855 "BB01AD.f"
			ios = s_rsle(&io___29);
#line 855 "BB01AD.f"
			if (ios != 0) {
#line 855 "BB01AD.f"
			    goto L100011;
#line 855 "BB01AD.f"
			}
#line 855 "BB01AD.f"
			i__2 = ipar[1];
#line 855 "BB01AD.f"
			for (j = 46; j <= i__2; ++j) {
#line 855 "BB01AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&a[j + i__ * 
				    a_dim1], (ftnlen)sizeof(doublereal));
#line 855 "BB01AD.f"
			    if (ios != 0) {
#line 855 "BB01AD.f"
				goto L100011;
#line 855 "BB01AD.f"
			    }
#line 855 "BB01AD.f"
			}
#line 855 "BB01AD.f"
			ios = e_rsle();
#line 855 "BB01AD.f"
L100011:
#line 857 "BB01AD.f"
			if (ios != 0) {
#line 857 "BB01AD.f"
			    *info = 1;
#line 857 "BB01AD.f"
			}
#line 858 "BB01AD.f"
/* L42: */
#line 858 "BB01AD.f"
		    }
#line 859 "BB01AD.f"
		    a[a_dim1 * 29 + 29] = -5.301;
#line 860 "BB01AD.f"
		    i__1 = ipar[2];
#line 860 "BB01AD.f"
		    for (j = 1; j <= i__1; ++j) {
#line 861 "BB01AD.f"
			ios = s_rsle(&io___30);
#line 861 "BB01AD.f"
			if (ios != 0) {
#line 861 "BB01AD.f"
			    goto L100012;
#line 861 "BB01AD.f"
			}
#line 861 "BB01AD.f"
			i__2 = ipar[1];
#line 861 "BB01AD.f"
			for (i__ = 1; i__ <= i__2; ++i__) {
#line 861 "BB01AD.f"
			    ios = do_lio(&c__5, &c__1, (char *)&b[i__ + j * 
				    b_dim1], (ftnlen)sizeof(doublereal));
#line 861 "BB01AD.f"
			    if (ios != 0) {
#line 861 "BB01AD.f"
				goto L100012;
#line 861 "BB01AD.f"
			    }
#line 861 "BB01AD.f"
			}
#line 861 "BB01AD.f"
			ios = e_rsle();
#line 861 "BB01AD.f"
L100012:
#line 863 "BB01AD.f"
			if (ios != 0) {
#line 863 "BB01AD.f"
			    *info = 1;
#line 863 "BB01AD.f"
			}
#line 864 "BB01AD.f"
/* L43: */
#line 864 "BB01AD.f"
		    }
#line 865 "BB01AD.f"
		    g[1] = 6.85e-6;
#line 866 "BB01AD.f"
		    g[3] = 373.;
#line 867 "BB01AD.f"
		    c__[c_dim1 * 52 + 1] = .3713;
#line 868 "BB01AD.f"
		    c__[c_dim1 * 53 + 1] = 1.245;
#line 869 "BB01AD.f"
		    c__[c_dim1 * 48 + 2] = 8e5;
#line 870 "BB01AD.f"
		    c__[c_dim1 * 54 + 2] = 1.;
#line 871 "BB01AD.f"
		    c__[c_dim1 * 51 + 3] = 8e5;
#line 872 "BB01AD.f"
		    c__[c_dim1 * 55 + 3] = 1.;
#line 873 "BB01AD.f"
		    q[1] = 28224.;
#line 874 "BB01AD.f"
		    q[4] = 2.742e-5;
#line 875 "BB01AD.f"
		    q[6] = 6.854e-4;
#line 876 "BB01AD.f"
		}
#line 877 "BB01AD.f"
	    }
#line 878 "BB01AD.f"
	    cl__1.cerr = 0;
#line 878 "BB01AD.f"
	    cl__1.cunit = 1;
#line 878 "BB01AD.f"
	    cl__1.csta = 0;
#line 878 "BB01AD.f"
	    f_clos(&cl__1);
#line 879 "BB01AD.f"
	    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);
#line 880 "BB01AD.f"
	}

#line 882 "BB01AD.f"
    } else if (nr[1] == 3) {
#line 883 "BB01AD.f"
	if (nr[2] == 1) {
#line 884 "BB01AD.f"
	    i__1 = ipar[1];
#line 884 "BB01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 885 "BB01AD.f"
		if (i__ % 2 == 1) {
#line 886 "BB01AD.f"
		    a[i__ + i__ * a_dim1] = -1.;
#line 887 "BB01AD.f"
		    b[i__ + (i__ + 1) / 2 * b_dim1] = 1.;
#line 888 "BB01AD.f"
		} else {
#line 889 "BB01AD.f"
		    a[i__ + (i__ - 1) * a_dim1] = 1.;
#line 890 "BB01AD.f"
		    a[i__ + (i__ + 1) * a_dim1] = -1.;
#line 891 "BB01AD.f"
		    c__[i__ / 2 + i__ * c_dim1] = 1.;
#line 892 "BB01AD.f"
		}
#line 893 "BB01AD.f"
/* L45: */
#line 893 "BB01AD.f"
	    }
#line 894 "BB01AD.f"
	    isymm = 1;
#line 895 "BB01AD.f"
	    for (i__ = ipar[3]; i__ >= 1; --i__) {
#line 896 "BB01AD.f"
		q[isymm] = 10.;
#line 897 "BB01AD.f"
		isymm += i__;
#line 898 "BB01AD.f"
/* L50: */
#line 898 "BB01AD.f"
	    }
#line 899 "BB01AD.f"
	    s_copy(ident, "0001", (ftnlen)4, (ftnlen)4);

#line 901 "BB01AD.f"
	} else if (nr[2] == 2) {
#line 902 "BB01AD.f"
	    i__1 = ipar[1];
#line 902 "BB01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 903 "BB01AD.f"
		a[i__ + i__ * a_dim1] = -2.;
#line 904 "BB01AD.f"
		if (i__ < ipar[1]) {
#line 905 "BB01AD.f"
		    a[i__ + (i__ + 1) * a_dim1] = 1.;
#line 906 "BB01AD.f"
		    a[i__ + 1 + i__ * a_dim1] = 1.;
#line 907 "BB01AD.f"
		}
#line 908 "BB01AD.f"
/* L60: */
#line 908 "BB01AD.f"
	    }
#line 909 "BB01AD.f"
	    a[ipar[1] * a_dim1 + 1] = 1.;
#line 910 "BB01AD.f"
	    a[ipar[1] + a_dim1] = 1.;
#line 911 "BB01AD.f"
	    s_copy(ident, "1111", (ftnlen)4, (ftnlen)4);
#line 912 "BB01AD.f"
	    temp = 6.2831853071795862 / (doublereal) ipar[1];
#line 913 "BB01AD.f"
	    i__1 = ipar[1];
#line 913 "BB01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 914 "BB01AD.f"
		dwork[i__] = cos(temp * (doublereal) (i__ - 1));
#line 915 "BB01AD.f"
		dwork[ipar[1] + i__] = dwork[i__] * 2. - 2. + sqrt(dwork[i__] 
			* 4. * (dwork[i__] - 2.) + 5.);
#line 917 "BB01AD.f"
/* L70: */
#line 917 "BB01AD.f"
	    }
#line 918 "BB01AD.f"
	    i__1 = ipar[1];
#line 918 "BB01AD.f"
	    for (j = 1; j <= i__1; ++j) {
#line 919 "BB01AD.f"
		i__2 = ipar[1];
#line 919 "BB01AD.f"
		for (i__ = 1; i__ <= i__2; ++i__) {
#line 920 "BB01AD.f"
		    dwork[(ipar[1] << 1) + i__] = cos(temp * (doublereal) (
			    i__ - 1) * (doublereal) (j - 1));
#line 921 "BB01AD.f"
/* L80: */
#line 921 "BB01AD.f"
		}
#line 922 "BB01AD.f"
		x[j + x_dim1] = ddot_(&ipar[1], &dwork[ipar[1] + 1], &c__1, &
			dwork[(ipar[1] << 1) + 1], &c__1) / (doublereal) ipar[
			1];
#line 924 "BB01AD.f"
/* L90: */
#line 924 "BB01AD.f"
	    }
/*         .. set up circulant solution matrix .. */
#line 926 "BB01AD.f"
	    i__1 = ipar[1];
#line 926 "BB01AD.f"
	    for (i__ = 2; i__ <= i__1; ++i__) {
#line 927 "BB01AD.f"
		i__2 = ipar[1] - i__ + 1;
#line 927 "BB01AD.f"
		dcopy_(&i__2, &x[x_dim1 + 1], &c__1, &x[i__ + i__ * x_dim1], &
			c__1);
#line 928 "BB01AD.f"
		i__2 = i__ - 1;
#line 928 "BB01AD.f"
		dcopy_(&i__2, &x[ipar[1] - i__ + 2 + x_dim1], &c__1, &x[i__ * 
			x_dim1 + 1], &c__1);
#line 929 "BB01AD.f"
/* L100: */
#line 929 "BB01AD.f"
	    }
#line 930 "BB01AD.f"
	}

#line 932 "BB01AD.f"
    } else if (nr[1] == 4) {
#line 933 "BB01AD.f"
	if (nr[2] == 1) {
/*       .. set up remaining parameter .. */
#line 935 "BB01AD.f"
	    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
#line 936 "BB01AD.f"
		dpar[1] = 1.;
#line 937 "BB01AD.f"
		dpar[2] = 1.;
#line 938 "BB01AD.f"
	    }
#line 939 "BB01AD.f"
	    i__1 = ipar[1] - 1;
#line 939 "BB01AD.f"
	    i__2 = ipar[1] - 1;
#line 939 "BB01AD.f"
	    dlaset_("A", &i__1, &i__2, &c_b10, &c_b30, &a[(a_dim1 << 1) + 1], 
		    lda, (ftnlen)1);
#line 940 "BB01AD.f"
	    b[ipar[1] + b_dim1] = 1.;
#line 941 "BB01AD.f"
	    c__[c_dim1 + 1] = 1.;
#line 942 "BB01AD.f"
	    q[1] = dpar[1];
#line 943 "BB01AD.f"
	    g[1] = dpar[2];
#line 944 "BB01AD.f"
	    s_copy(ident, "0000", (ftnlen)4, (ftnlen)4);

#line 946 "BB01AD.f"
	} else if (nr[2] == 2) {
/*         .. set up remaining parameters .. */
#line 948 "BB01AD.f"
	    appind = (doublereal) (ipar[1] + 1);
#line 949 "BB01AD.f"
	    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
#line 950 "BB01AD.f"
		dpar[1] = pardef[nr[1] + (nr[2] << 2) - 5];
#line 951 "BB01AD.f"
		dpar[2] = 1.;
#line 952 "BB01AD.f"
		dpar[3] = 1.;
#line 953 "BB01AD.f"
		dpar[4] = .2;
#line 954 "BB01AD.f"
		dpar[5] = .3;
#line 955 "BB01AD.f"
		dpar[6] = .2;
#line 956 "BB01AD.f"
		dpar[7] = .3;
#line 957 "BB01AD.f"
	    }
/*         .. set up stiffness matrix .. */
#line 959 "BB01AD.f"
	    temp = -dpar[1] * appind;
#line 960 "BB01AD.f"
	    d__1 = temp * 2.;
#line 960 "BB01AD.f"
	    dlaset_("A", &ipar[1], &ipar[1], &c_b10, &d__1, &a[a_offset], lda,
		     (ftnlen)1);
#line 961 "BB01AD.f"
	    i__1 = ipar[1] - 1;
#line 961 "BB01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
#line 962 "BB01AD.f"
		a[i__ + 1 + i__ * a_dim1] = -temp;
#line 963 "BB01AD.f"
		a[i__ + (i__ + 1) * a_dim1] = -temp;
#line 964 "BB01AD.f"
/* L110: */
#line 964 "BB01AD.f"
	    }
/*         .. set up Gramian, stored by diagonals .. */
#line 966 "BB01AD.f"
	    temp = 1. / (appind * 6.);
#line 967 "BB01AD.f"
	    d__1 = temp * 4.;
#line 967 "BB01AD.f"
	    d__2 = temp * 4.;
#line 967 "BB01AD.f"
	    dlaset_("L", &ipar[1], &c__1, &d__1, &d__2, &dwork[1], &ipar[1], (
		    ftnlen)1);
#line 969 "BB01AD.f"
	    i__1 = ipar[1] - 1;
#line 969 "BB01AD.f"
	    dlaset_("L", &i__1, &c__1, &temp, &temp, &dwork[ipar[1] + 1], &
		    ipar[1], (ftnlen)1);
#line 971 "BB01AD.f"
	    dpttrf_(&ipar[1], &dwork[1], &dwork[ipar[1] + 1], info);
/*         .. A = (inverse of Gramian) * (stiffness matrix) .. */
#line 973 "BB01AD.f"
	    dpttrs_(&ipar[1], &ipar[1], &dwork[1], &dwork[ipar[1] + 1], &a[
		    a_offset], lda, info);
/*         .. compute B, C .. */
#line 976 "BB01AD.f"
	    i__1 = ipar[1];
#line 976 "BB01AD.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
#line 977 "BB01AD.f"
		d__1 = (doublereal) (i__ - 1) / appind;
#line 977 "BB01AD.f"
		b1 = max(d__1,dpar[4]);
/* Computing MIN */
#line 978 "BB01AD.f"
		d__1 = (doublereal) (i__ + 1) / appind;
#line 978 "BB01AD.f"
		b2 = min(d__1,dpar[5]);
/* Computing MAX */
#line 979 "BB01AD.f"
		d__1 = (doublereal) (i__ - 1) / appind;
#line 979 "BB01AD.f"
		c1 = max(d__1,dpar[6]);
/* Computing MIN */
#line 980 "BB01AD.f"
		d__1 = (doublereal) (i__ + 1) / appind;
#line 980 "BB01AD.f"
		c2 = min(d__1,dpar[7]);
#line 981 "BB01AD.f"
		if (b1 >= b2) {
#line 982 "BB01AD.f"
		    b[i__ + b_dim1] = 0.;
#line 983 "BB01AD.f"
		} else {
#line 984 "BB01AD.f"
		    b[i__ + b_dim1] = b2 - b1;
/* Computing MIN */
#line 985 "BB01AD.f"
		    d__1 = b2, d__2 = (doublereal) i__ / appind;
#line 985 "BB01AD.f"
		    temp = min(d__1,d__2);
#line 986 "BB01AD.f"
		    if (b1 < temp) {
/* Computing 2nd power */
#line 987 "BB01AD.f"
			d__1 = temp;
/* Computing 2nd power */
#line 987 "BB01AD.f"
			d__2 = b1;
#line 987 "BB01AD.f"
			b[i__ + b_dim1] += appind * (d__1 * d__1 - d__2 * 
				d__2) / 2.;
#line 988 "BB01AD.f"
			b[i__ + b_dim1] += (doublereal) i__ * (b1 - temp);
#line 989 "BB01AD.f"
		    }
/* Computing MAX */
#line 990 "BB01AD.f"
		    d__1 = b1, d__2 = (doublereal) i__ / appind;
#line 990 "BB01AD.f"
		    temp = max(d__1,d__2);
#line 991 "BB01AD.f"
		    if (temp < b2) {
/* Computing 2nd power */
#line 992 "BB01AD.f"
			d__1 = b2;
/* Computing 2nd power */
#line 992 "BB01AD.f"
			d__2 = temp;
#line 992 "BB01AD.f"
			b[i__ + b_dim1] -= appind * (d__1 * d__1 - d__2 * 
				d__2) / 2.;
#line 993 "BB01AD.f"
			b[i__ + b_dim1] -= (doublereal) i__ * (temp - b2);
#line 994 "BB01AD.f"
		    }
#line 995 "BB01AD.f"
		}
#line 996 "BB01AD.f"
		if (c1 >= c2) {
#line 997 "BB01AD.f"
		    c__[i__ * c_dim1 + 1] = 0.;
#line 998 "BB01AD.f"
		} else {
#line 999 "BB01AD.f"
		    c__[i__ * c_dim1 + 1] = c2 - c1;
/* Computing MIN */
#line 1000 "BB01AD.f"
		    d__1 = c2, d__2 = (doublereal) i__ / appind;
#line 1000 "BB01AD.f"
		    temp = min(d__1,d__2);
#line 1001 "BB01AD.f"
		    if (c1 < temp) {
/* Computing 2nd power */
#line 1002 "BB01AD.f"
			d__1 = temp;
/* Computing 2nd power */
#line 1002 "BB01AD.f"
			d__2 = c1;
#line 1002 "BB01AD.f"
			c__[i__ * c_dim1 + 1] += appind * (d__1 * d__1 - d__2 
				* d__2) / 2.;
#line 1003 "BB01AD.f"
			c__[i__ * c_dim1 + 1] += (doublereal) i__ * (c1 - 
				temp);
#line 1004 "BB01AD.f"
		    }
/* Computing MAX */
#line 1005 "BB01AD.f"
		    d__1 = c1, d__2 = (doublereal) i__ / appind;
#line 1005 "BB01AD.f"
		    temp = max(d__1,d__2);
#line 1006 "BB01AD.f"
		    if (temp < c2) {
/* Computing 2nd power */
#line 1007 "BB01AD.f"
			d__1 = c2;
/* Computing 2nd power */
#line 1007 "BB01AD.f"
			d__2 = temp;
#line 1007 "BB01AD.f"
			c__[i__ * c_dim1 + 1] -= appind * (d__1 * d__1 - d__2 
				* d__2) / 2.;
#line 1008 "BB01AD.f"
			c__[i__ * c_dim1 + 1] -= (doublereal) i__ * (temp - 
				c2);
#line 1009 "BB01AD.f"
		    }
#line 1010 "BB01AD.f"
		}
#line 1011 "BB01AD.f"
/* L120: */
#line 1011 "BB01AD.f"
	    }
#line 1012 "BB01AD.f"
	    dscal_(&ipar[1], &dpar[2], &b[b_dim1 + 1], &c__1);
#line 1013 "BB01AD.f"
	    dscal_(&ipar[1], &dpar[3], &c__[c_dim1 + 1], ldc);
#line 1014 "BB01AD.f"
	    dpttrs_(&ipar[1], &c__1, &dwork[1], &dwork[ipar[1] + 1], &b[
		    b_offset], ldb, info);
#line 1016 "BB01AD.f"
	    s_copy(ident, "0011", (ftnlen)4, (ftnlen)4);

#line 1018 "BB01AD.f"
	} else if (nr[2] == 3) {
/*         .. set up remaining parameters .. */
#line 1020 "BB01AD.f"
	    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
#line 1021 "BB01AD.f"
		dpar[1] = pardef[nr[1] + (nr[2] << 2) - 5];
#line 1022 "BB01AD.f"
		dpar[2] = 4.;
#line 1023 "BB01AD.f"
		dpar[3] = 1.;
#line 1024 "BB01AD.f"
	    }
#line 1025 "BB01AD.f"
	    if (dpar[1] != 0.) {
#line 1026 "BB01AD.f"
		dlaset_("A", &l, &l, &c_b10, &c_b30, &a[(l + 1) * a_dim1 + 1],
			 lda, (ftnlen)1);
#line 1027 "BB01AD.f"
		temp = dpar[3] / dpar[1];
#line 1028 "BB01AD.f"
		a[l + 1 + a_dim1] = -temp;
#line 1029 "BB01AD.f"
		a[l + 1 + (a_dim1 << 1)] = temp;
#line 1030 "BB01AD.f"
		a[ipar[1] + (l - 1) * a_dim1] = temp;
#line 1031 "BB01AD.f"
		a[ipar[1] + l * a_dim1] = -temp;
#line 1032 "BB01AD.f"
		ttemp = temp * 2.;
#line 1033 "BB01AD.f"
		i__1 = l - 1;
#line 1033 "BB01AD.f"
		for (i__ = 2; i__ <= i__1; ++i__) {
#line 1034 "BB01AD.f"
		    a[l + i__ + i__ * a_dim1] = -ttemp;
#line 1035 "BB01AD.f"
		    a[l + i__ + (i__ + 1) * a_dim1] = temp;
#line 1036 "BB01AD.f"
		    a[l + i__ + (i__ - 1) * a_dim1] = temp;
#line 1037 "BB01AD.f"
/* L130: */
#line 1037 "BB01AD.f"
		}
#line 1038 "BB01AD.f"
		d__1 = -dpar[2] / dpar[1];
#line 1038 "BB01AD.f"
		dlaset_("A", &l, &l, &c_b10, &d__1, &a[l + 1 + (l + 1) * 
			a_dim1], lda, (ftnlen)1);
#line 1040 "BB01AD.f"
		b[l + 1 + b_dim1] = 1. / dpar[1];
#line 1041 "BB01AD.f"
		b[ipar[1] + ipar[2] * b_dim1] = -1. / dpar[1];
#line 1042 "BB01AD.f"
		s_copy(ident, "0111", (ftnlen)4, (ftnlen)4);
#line 1043 "BB01AD.f"
	    } else {
#line 1044 "BB01AD.f"
		*info = 2;
#line 1045 "BB01AD.f"
	    }

#line 1047 "BB01AD.f"
	} else if (nr[2] == 4) {
#line 1048 "BB01AD.f"
	    if (! lsame_(def, "N", (ftnlen)1, (ftnlen)1)) {
#line 1048 "BB01AD.f"
		ici__1.icierr = 0;
#line 1048 "BB01AD.f"
		ici__1.icirnum = 1;
#line 1048 "BB01AD.f"
		ici__1.icirlen = 11;
#line 1048 "BB01AD.f"
		ici__1.iciunit = chpar;
#line 1048 "BB01AD.f"
		ici__1.icifmt = "(A,I1,A,I1,A)";
#line 1048 "BB01AD.f"
		s_wsfi(&ici__1);
#line 1048 "BB01AD.f"
		do_fio(&c__1, "BB01", (ftnlen)4);
#line 1048 "BB01AD.f"
		do_fio(&c__1, (char *)&nr[1], (ftnlen)sizeof(integer));
#line 1048 "BB01AD.f"
		do_fio(&c__1, "0", (ftnlen)1);
#line 1048 "BB01AD.f"
		do_fio(&c__1, (char *)&nr[2], (ftnlen)sizeof(integer));
#line 1048 "BB01AD.f"
		do_fio(&c__1, ".dat", (ftnlen)4);
#line 1048 "BB01AD.f"
		e_wsfi();
#line 1048 "BB01AD.f"
	    }
#line 1050 "BB01AD.f"
	    o__1.oerr = 1;
#line 1050 "BB01AD.f"
	    o__1.ounit = 1;
#line 1050 "BB01AD.f"
	    o__1.ofnmlen = 11;
#line 1050 "BB01AD.f"
	    o__1.ofnm = chpar;
#line 1050 "BB01AD.f"
	    o__1.orl = 0;
#line 1050 "BB01AD.f"
	    o__1.osta = "OLD";
#line 1050 "BB01AD.f"
	    o__1.oacc = 0;
#line 1050 "BB01AD.f"
	    o__1.ofm = 0;
#line 1050 "BB01AD.f"
	    o__1.oblnk = 0;
#line 1050 "BB01AD.f"
	    ios = f_open(&o__1);
#line 1051 "BB01AD.f"
	    if (ios != 0) {
#line 1052 "BB01AD.f"
		*info = 1;
#line 1053 "BB01AD.f"
	    } else {
#line 1054 "BB01AD.f"
		ios = s_rsle(&io___37);
#line 1054 "BB01AD.f"
		if (ios != 0) {
#line 1054 "BB01AD.f"
		    goto L100013;
#line 1054 "BB01AD.f"
		}
#line 1054 "BB01AD.f"
		i__1 = (l << 2) - 2;
#line 1054 "BB01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 1054 "BB01AD.f"
		    ios = do_lio(&c__5, &c__1, (char *)&dwork[i__], (ftnlen)
			    sizeof(doublereal));
#line 1054 "BB01AD.f"
		    if (ios != 0) {
#line 1054 "BB01AD.f"
			goto L100013;
#line 1054 "BB01AD.f"
		    }
#line 1054 "BB01AD.f"
		}
#line 1054 "BB01AD.f"
		ios = e_rsle();
#line 1054 "BB01AD.f"
L100013:
#line 1055 "BB01AD.f"
		if (ios != 0) {
#line 1055 "BB01AD.f"
		    *info = 1;
#line 1055 "BB01AD.f"
		}
#line 1056 "BB01AD.f"
	    }
#line 1057 "BB01AD.f"
	    cl__1.cerr = 0;
#line 1057 "BB01AD.f"
	    cl__1.cunit = 1;
#line 1057 "BB01AD.f"
	    cl__1.csta = 0;
#line 1057 "BB01AD.f"
	    f_clos(&cl__1);
#line 1058 "BB01AD.f"
	    if (*info == 0) {
#line 1059 "BB01AD.f"
		i__1 = l - 1;
#line 1059 "BB01AD.f"
		i__2 = l - 1;
#line 1059 "BB01AD.f"
		dlaset_("A", &i__1, &i__2, &c_b10, &c_b30, &a[l + 1 + (a_dim1 
			<< 1)], lda, (ftnlen)1);
#line 1060 "BB01AD.f"
		pos = (l << 1) + 1;
#line 1061 "BB01AD.f"
		a[(a_dim1 << 1) + 1] = -dwork[pos] / dwork[1];
#line 1062 "BB01AD.f"
		i__1 = l;
#line 1062 "BB01AD.f"
		for (i__ = 2; i__ <= i__1; ++i__) {
#line 1063 "BB01AD.f"
		    temp = dwork[pos] / dwork[i__ - 1];
#line 1064 "BB01AD.f"
		    ttemp = dwork[pos] / dwork[i__];
#line 1065 "BB01AD.f"
		    if (i__ > 2) {
#line 1065 "BB01AD.f"
			a[i__ - 1 + i__ * a_dim1] = temp;
#line 1065 "BB01AD.f"
		    }
#line 1066 "BB01AD.f"
		    a[i__ + i__ * a_dim1] = -(temp + ttemp);
#line 1067 "BB01AD.f"
		    if (i__ < l) {
#line 1067 "BB01AD.f"
			a[i__ + 1 + i__ * a_dim1] = ttemp;
#line 1067 "BB01AD.f"
		    }
#line 1068 "BB01AD.f"
		    ++pos;
#line 1069 "BB01AD.f"
/* L140: */
#line 1069 "BB01AD.f"
		}
#line 1070 "BB01AD.f"
		pos = l;
#line 1071 "BB01AD.f"
		temp = dwork[pos + 1] / dwork[1];
#line 1072 "BB01AD.f"
		a[a_dim1 + 1] = -temp;
#line 1073 "BB01AD.f"
		i__1 = l;
#line 1073 "BB01AD.f"
		for (i__ = 2; i__ <= i__1; ++i__) {
#line 1074 "BB01AD.f"
		    ttemp = temp;
#line 1075 "BB01AD.f"
		    temp = dwork[pos + i__] / dwork[i__];
#line 1076 "BB01AD.f"
		    sum = ttemp - temp;
#line 1077 "BB01AD.f"
		    a[i__ + a_dim1] = -sum;
#line 1078 "BB01AD.f"
		    a[i__ + i__ * a_dim1] -= temp;
#line 1079 "BB01AD.f"
		    i__2 = i__ - 2;
#line 1079 "BB01AD.f"
		    for (j = 2; j <= i__2; ++j) {
#line 1080 "BB01AD.f"
			a[i__ + j * a_dim1] = sum;
#line 1081 "BB01AD.f"
/* L150: */
#line 1081 "BB01AD.f"
		    }
#line 1082 "BB01AD.f"
		    if (i__ > 2) {
#line 1082 "BB01AD.f"
			a[i__ + (i__ - 1) * a_dim1] += sum;
#line 1082 "BB01AD.f"
		    }
#line 1083 "BB01AD.f"
/* L160: */
#line 1083 "BB01AD.f"
		}
#line 1084 "BB01AD.f"
		pos = l * 3;
#line 1085 "BB01AD.f"
		a[(l + 1) * a_dim1 + 1] = -dwork[l * 3] / dwork[1];
#line 1086 "BB01AD.f"
		i__1 = l;
#line 1086 "BB01AD.f"
		for (i__ = 2; i__ <= i__1; ++i__) {
#line 1087 "BB01AD.f"
		    temp = dwork[pos] / dwork[i__ - 1];
#line 1088 "BB01AD.f"
		    ttemp = dwork[pos] / dwork[i__];
#line 1089 "BB01AD.f"
		    if (i__ > 2) {
#line 1089 "BB01AD.f"
			a[i__ - 1 + (l + i__ - 1) * a_dim1] = temp;
#line 1089 "BB01AD.f"
		    }
#line 1090 "BB01AD.f"
		    a[i__ + (l + i__ - 1) * a_dim1] = -(temp + ttemp);
#line 1091 "BB01AD.f"
		    if (i__ < l) {
#line 1091 "BB01AD.f"
			a[i__ + 1 + (l + i__ - 1) * a_dim1] = ttemp;
#line 1091 "BB01AD.f"
		    }
#line 1092 "BB01AD.f"
		    ++pos;
#line 1093 "BB01AD.f"
/* L170: */
#line 1093 "BB01AD.f"
		}
#line 1094 "BB01AD.f"
		b[b_dim1 + 1] = 1. / dwork[1];
#line 1095 "BB01AD.f"
		i__1 = l;
#line 1095 "BB01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 1096 "BB01AD.f"
		    temp = 1. / dwork[i__];
#line 1097 "BB01AD.f"
		    if (i__ > 1) {
#line 1097 "BB01AD.f"
			b[i__ + i__ * b_dim1] = -temp;
#line 1097 "BB01AD.f"
		    }
#line 1098 "BB01AD.f"
		    if (i__ < l) {
#line 1098 "BB01AD.f"
			b[i__ + 1 + i__ * b_dim1] = temp;
#line 1098 "BB01AD.f"
		    }
#line 1099 "BB01AD.f"
/* L180: */
#line 1099 "BB01AD.f"
		}
#line 1100 "BB01AD.f"
		c__[c_dim1 + 1] = 1.;
#line 1101 "BB01AD.f"
		q[1] = 1.;
#line 1102 "BB01AD.f"
		pos = (l << 1) - 1;
#line 1103 "BB01AD.f"
		isymm = l + 1;
#line 1104 "BB01AD.f"
		i__1 = l;
#line 1104 "BB01AD.f"
		for (i__ = 2; i__ <= i__1; ++i__) {
#line 1105 "BB01AD.f"
		    temp = dwork[pos + i__];
#line 1106 "BB01AD.f"
		    ttemp = dwork[pos + l + i__ - 1];
#line 1107 "BB01AD.f"
		    c__[i__ + i__ * c_dim1] = temp;
#line 1108 "BB01AD.f"
		    c__[i__ + (l + i__ - 1) * c_dim1] = ttemp;
#line 1109 "BB01AD.f"
		    q[isymm] = 1. / (temp * temp + ttemp * ttemp);
#line 1110 "BB01AD.f"
		    isymm = isymm + l - i__ + 1;
#line 1111 "BB01AD.f"
/* L190: */
#line 1111 "BB01AD.f"
		}
#line 1112 "BB01AD.f"
		s_copy(ident, "0001", (ftnlen)4, (ftnlen)4);
#line 1113 "BB01AD.f"
	    }
#line 1114 "BB01AD.f"
	}
#line 1115 "BB01AD.f"
    }

#line 1117 "BB01AD.f"
    if (*info != 0) {
#line 1117 "BB01AD.f"
	goto L2001;
#line 1117 "BB01AD.f"
    }
/*     .. set up data in required format .. */

#line 1120 "BB01AD.f"
    if (bpar[1]) {
/*     .. G is to be returned in product form .. */
#line 1122 "BB01AD.f"
	gdimm = ipar[1];
#line 1123 "BB01AD.f"
	if (*(unsigned char *)&ident[3] == '0') {
/*       .. invert R using Cholesky factorization, store in G .. */
#line 1125 "BB01AD.f"
	    dpptrf_("L", &ipar[2], &g[1], info, (ftnlen)1);
#line 1126 "BB01AD.f"
	    if (*info == 0) {
#line 1127 "BB01AD.f"
		dpptri_("L", &ipar[2], &g[1], info, (ftnlen)1);
#line 1128 "BB01AD.f"
		if (*(unsigned char *)ident == '0') {
/*         .. B is not identity matrix .. */
#line 1130 "BB01AD.f"
		    i__1 = ipar[1];
#line 1130 "BB01AD.f"
		    for (i__ = 1; i__ <= i__1; ++i__) {
#line 1131 "BB01AD.f"
			dspmv_("L", &ipar[2], &c_b30, &g[1], &b[i__ + b_dim1],
				 ldb, &c_b10, &dwork[(i__ - 1) * ipar[1] + 1],
				 &c__1, (ftnlen)1);
#line 1133 "BB01AD.f"
/* L200: */
#line 1133 "BB01AD.f"
		    }
#line 1134 "BB01AD.f"
		    dgemv_("T", &ipar[2], &ipar[1], &c_b30, &dwork[1], &ipar[
			    1], &b[b_dim1 + 1], ldb, &c_b10, &g[1], &c__1, (
			    ftnlen)1);
#line 1136 "BB01AD.f"
		    isymm = ipar[1] + 1;
#line 1137 "BB01AD.f"
		    i__1 = ipar[1];
#line 1137 "BB01AD.f"
		    for (i__ = 2; i__ <= i__1; ++i__) {
#line 1138 "BB01AD.f"
			dgemv_("T", &ipar[2], &ipar[1], &c_b30, &dwork[1], &
				ipar[1], &b[i__ + b_dim1], ldb, &c_b10, &b[
				b_dim1 + 1], ldb, (ftnlen)1);
#line 1140 "BB01AD.f"
			i__2 = ipar[1] - i__ + 1;
#line 1140 "BB01AD.f"
			dcopy_(&i__2, &b[i__ * b_dim1 + 1], ldb, &g[isymm], &
				c__1);
#line 1141 "BB01AD.f"
			isymm += ipar[1] - i__ + 1;
#line 1142 "BB01AD.f"
/* L210: */
#line 1142 "BB01AD.f"
		    }
#line 1143 "BB01AD.f"
		}
#line 1144 "BB01AD.f"
	    } else {
#line 1145 "BB01AD.f"
		if (*info > 0) {
#line 1146 "BB01AD.f"
		    *info = 3;
#line 1147 "BB01AD.f"
		    goto L2001;
#line 1148 "BB01AD.f"
		}
#line 1149 "BB01AD.f"
	    }
#line 1150 "BB01AD.f"
	} else {
/*       .. R = identity .. */
#line 1152 "BB01AD.f"
	    if (*(unsigned char *)ident == '0') {
/*         .. B is not identity matrix .. */
#line 1154 "BB01AD.f"
		if (ipar[2] == 1) {
#line 1155 "BB01AD.f"
		    dlaset_("L", &nsymm, &c__1, &c_b10, &c_b10, &g[1], &c__1, 
			    (ftnlen)1);
#line 1156 "BB01AD.f"
		    dspr_("L", &ipar[1], &c_b30, &b[b_offset], &c__1, &g[1], (
			    ftnlen)1);
#line 1157 "BB01AD.f"
		} else {
#line 1158 "BB01AD.f"
		    dsyrk_("L", "N", &ipar[1], &ipar[2], &c_b30, &b[b_offset],
			     ldb, &c_b10, &dwork[1], &ipar[1], (ftnlen)1, (
			    ftnlen)1);
#line 1160 "BB01AD.f"
		    ma02dd_("Pack", "Lower", &ipar[1], &dwork[1], &ipar[1], &
			    g[1], (ftnlen)4, (ftnlen)5);
#line 1161 "BB01AD.f"
		}
#line 1162 "BB01AD.f"
	    } else {
/*         .. B = R = identity .. */
#line 1164 "BB01AD.f"
		isymm = 1;
#line 1165 "BB01AD.f"
		for (i__ = ipar[1]; i__ >= 1; --i__) {
#line 1166 "BB01AD.f"
		    g[isymm] = 1.;
#line 1167 "BB01AD.f"
		    isymm += i__;
#line 1168 "BB01AD.f"
/* L220: */
#line 1168 "BB01AD.f"
		}
#line 1169 "BB01AD.f"
	    }
#line 1170 "BB01AD.f"
	}
#line 1171 "BB01AD.f"
    } else {
#line 1172 "BB01AD.f"
	gdimm = ipar[2];
#line 1173 "BB01AD.f"
	if (*(unsigned char *)ident == '1') {
#line 1173 "BB01AD.f"
	    dlaset_("A", &ipar[1], &ipar[2], &c_b10, &c_b30, &b[b_offset], 
		    ldb, (ftnlen)1);
#line 1173 "BB01AD.f"
	}
#line 1175 "BB01AD.f"
	if (*(unsigned char *)&ident[3] == '1') {
#line 1176 "BB01AD.f"
	    isymm = 1;
#line 1177 "BB01AD.f"
	    for (i__ = ipar[2]; i__ >= 1; --i__) {
#line 1178 "BB01AD.f"
		g[isymm] = 1.;
#line 1179 "BB01AD.f"
		isymm += i__;
#line 1180 "BB01AD.f"
/* L230: */
#line 1180 "BB01AD.f"
	    }
#line 1181 "BB01AD.f"
	}
#line 1182 "BB01AD.f"
    }

#line 1184 "BB01AD.f"
    if (bpar[4]) {
/*     .. Q is to be returned in product form .. */
#line 1186 "BB01AD.f"
	qdimm = ipar[1];
#line 1187 "BB01AD.f"
	if (*(unsigned char *)&ident[2] == '0') {
#line 1188 "BB01AD.f"
	    if (*(unsigned char *)&ident[1] == '0') {
/*         .. C is not identity matrix .. */
#line 1190 "BB01AD.f"
		i__1 = ipar[1];
#line 1190 "BB01AD.f"
		for (i__ = 1; i__ <= i__1; ++i__) {
#line 1191 "BB01AD.f"
		    dspmv_("L", &ipar[3], &c_b30, &q[1], &c__[i__ * c_dim1 + 
			    1], &c__1, &c_b10, &dwork[(i__ - 1) * ipar[1] + 1]
			    , &c__1, (ftnlen)1);
#line 1193 "BB01AD.f"
/* L240: */
#line 1193 "BB01AD.f"
		}
/*         .. use Q(1:IPAR(1)) as workspace and compute the first column */
/*            of Q in the end .. */
#line 1196 "BB01AD.f"
		isymm = ipar[1] + 1;
#line 1197 "BB01AD.f"
		i__1 = ipar[1];
#line 1197 "BB01AD.f"
		for (i__ = 2; i__ <= i__1; ++i__) {
#line 1198 "BB01AD.f"
		    dgemv_("T", &ipar[3], &ipar[1], &c_b30, &dwork[1], &ipar[
			    1], &c__[i__ * c_dim1 + 1], &c__1, &c_b10, &q[1], 
			    &c__1, (ftnlen)1);
#line 1200 "BB01AD.f"
		    i__2 = ipar[1] - i__ + 1;
#line 1200 "BB01AD.f"
		    dcopy_(&i__2, &q[i__], &c__1, &q[isymm], &c__1);
#line 1201 "BB01AD.f"
		    isymm += ipar[1] - i__ + 1;
#line 1202 "BB01AD.f"
/* L250: */
#line 1202 "BB01AD.f"
		}
#line 1203 "BB01AD.f"
		dgemv_("T", &ipar[3], &ipar[1], &c_b30, &dwork[1], &ipar[1], &
			c__[c_dim1 + 1], &c__1, &c_b10, &q[1], &c__1, (ftnlen)
			1);
#line 1205 "BB01AD.f"
	    }
#line 1206 "BB01AD.f"
	} else {
/*       .. Q = identity .. */
#line 1208 "BB01AD.f"
	    if (*(unsigned char *)&ident[1] == '0') {
/*         .. C is not identity matrix .. */
#line 1210 "BB01AD.f"
		if (ipar[3] == 1) {
#line 1211 "BB01AD.f"
		    dlaset_("L", &nsymm, &c__1, &c_b10, &c_b10, &q[1], &c__1, 
			    (ftnlen)1);
#line 1212 "BB01AD.f"
		    dspr_("L", &ipar[1], &c_b30, &c__[c_offset], ldc, &q[1], (
			    ftnlen)1);
#line 1213 "BB01AD.f"
		} else {
#line 1214 "BB01AD.f"
		    dsyrk_("L", "T", &ipar[1], &ipar[3], &c_b30, &c__[
			    c_offset], ldc, &c_b10, &dwork[1], &ipar[1], (
			    ftnlen)1, (ftnlen)1);
#line 1216 "BB01AD.f"
		    ma02dd_("Pack", "Lower", &ipar[1], &dwork[1], &ipar[1], &
			    q[1], (ftnlen)4, (ftnlen)5);
#line 1217 "BB01AD.f"
		}
#line 1218 "BB01AD.f"
	    } else {
/*         .. C = Q = identity .. */
#line 1220 "BB01AD.f"
		isymm = 1;
#line 1221 "BB01AD.f"
		for (i__ = ipar[1]; i__ >= 1; --i__) {
#line 1222 "BB01AD.f"
		    q[isymm] = 1.;
#line 1223 "BB01AD.f"
		    isymm += i__;
#line 1224 "BB01AD.f"
/* L260: */
#line 1224 "BB01AD.f"
		}
#line 1225 "BB01AD.f"
	    }
#line 1226 "BB01AD.f"
	}
#line 1227 "BB01AD.f"
    } else {
#line 1228 "BB01AD.f"
	qdimm = ipar[3];
#line 1229 "BB01AD.f"
	if (*(unsigned char *)&ident[1] == '1') {
#line 1229 "BB01AD.f"
	    dlaset_("A", &ipar[3], &ipar[1], &c_b10, &c_b30, &c__[c_offset], 
		    ldc, (ftnlen)1);
#line 1229 "BB01AD.f"
	}
#line 1231 "BB01AD.f"
	if (*(unsigned char *)&ident[2] == '1') {
#line 1232 "BB01AD.f"
	    isymm = 1;
#line 1233 "BB01AD.f"
	    for (i__ = ipar[3]; i__ >= 1; --i__) {
#line 1234 "BB01AD.f"
		q[isymm] = 1.;
#line 1235 "BB01AD.f"
		isymm += i__;
#line 1236 "BB01AD.f"
/* L270: */
#line 1236 "BB01AD.f"
	    }
#line 1237 "BB01AD.f"
	}
#line 1238 "BB01AD.f"
    }

/*     .. unpack symmetric matrices if desired .. */
#line 1241 "BB01AD.f"
    if (bpar[2]) {
#line 1242 "BB01AD.f"
	isymm = gdimm * (gdimm + 1) / 2;
#line 1243 "BB01AD.f"
	dcopy_(&isymm, &g[1], &c__1, &dwork[1], &c__1);
#line 1244 "BB01AD.f"
	ma02dd_("Unpack", "Lower", &gdimm, &g[1], ldg, &dwork[1], (ftnlen)6, (
		ftnlen)5);
#line 1245 "BB01AD.f"
	ma02ed_("Lower", &gdimm, &g[1], ldg, (ftnlen)5);
#line 1246 "BB01AD.f"
    } else if (bpar[3]) {
#line 1247 "BB01AD.f"
	ma02dd_("Unpack", "Lower", &gdimm, &dwork[1], &gdimm, &g[1], (ftnlen)
		6, (ftnlen)5);
#line 1248 "BB01AD.f"
	ma02ed_("Lower", &gdimm, &dwork[1], &gdimm, (ftnlen)5);
#line 1249 "BB01AD.f"
	ma02dd_("Pack", "Upper", &gdimm, &dwork[1], &gdimm, &g[1], (ftnlen)4, 
		(ftnlen)5);
#line 1250 "BB01AD.f"
    }
#line 1251 "BB01AD.f"
    if (bpar[5]) {
#line 1252 "BB01AD.f"
	isymm = qdimm * (qdimm + 1) / 2;
#line 1253 "BB01AD.f"
	dcopy_(&isymm, &q[1], &c__1, &dwork[1], &c__1);
#line 1254 "BB01AD.f"
	ma02dd_("Unpack", "Lower", &qdimm, &q[1], ldq, &dwork[1], (ftnlen)6, (
		ftnlen)5);
#line 1255 "BB01AD.f"
	ma02ed_("Lower", &qdimm, &q[1], ldq, (ftnlen)5);
#line 1256 "BB01AD.f"
    } else if (bpar[6]) {
#line 1257 "BB01AD.f"
	ma02dd_("Unpack", "Lower", &qdimm, &dwork[1], &qdimm, &q[1], (ftnlen)
		6, (ftnlen)5);
#line 1258 "BB01AD.f"
	ma02ed_("Lower", &qdimm, &dwork[1], &qdimm, (ftnlen)5);
#line 1259 "BB01AD.f"
	ma02dd_("Pack", "Upper", &qdimm, &dwork[1], &qdimm, &q[1], (ftnlen)4, 
		(ftnlen)5);
#line 1260 "BB01AD.f"
    }

/*     ...set VEC... */
#line 1263 "BB01AD.f"
    vec[1] = TRUE_;
#line 1264 "BB01AD.f"
    vec[2] = TRUE_;
#line 1265 "BB01AD.f"
    vec[3] = TRUE_;
#line 1266 "BB01AD.f"
    vec[4] = TRUE_;
#line 1267 "BB01AD.f"
    vec[5] = ! bpar[1];
#line 1268 "BB01AD.f"
    vec[6] = ! bpar[4];
#line 1269 "BB01AD.f"
    vec[7] = TRUE_;
#line 1270 "BB01AD.f"
    vec[8] = TRUE_;
#line 1271 "BB01AD.f"
    if (nr[1] == 1) {
#line 1272 "BB01AD.f"
	if (nr[2] == 1 || nr[2] == 2) {
#line 1272 "BB01AD.f"
	    vec[9] = TRUE_;
#line 1272 "BB01AD.f"
	}
#line 1273 "BB01AD.f"
    } else if (nr[1] == 2) {
#line 1274 "BB01AD.f"
	if (nr[2] == 1 || nr[2] >= 3 && nr[2] <= 6) {
#line 1274 "BB01AD.f"
	    vec[9] = TRUE_;
#line 1274 "BB01AD.f"
	}
#line 1276 "BB01AD.f"
    } else if (nr[1] == 3) {
#line 1277 "BB01AD.f"
	if (nr[2] == 2) {
#line 1277 "BB01AD.f"
	    vec[9] = TRUE_;
#line 1277 "BB01AD.f"
	}
#line 1278 "BB01AD.f"
    }
#line 1279 "BB01AD.f"
    s_copy(chpar, notes + (nr[1] + (nr[2] << 2) - 5) * 255, (ftnlen)255, (
	    ftnlen)255);
#line 1280 "BB01AD.f"
    *n = ipar[1];
#line 1281 "BB01AD.f"
    *m = ipar[2];
#line 1282 "BB01AD.f"
    *p = ipar[3];
#line 1283 "BB01AD.f"
L2001:
#line 1284 "BB01AD.f"
    return 0;
/* *** Last line of BB01AD *** */
} /* bb01ad_ */

#undef notes


