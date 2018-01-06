#line 1 "MB03TS.f"
/* MB03TS.f -- translated by f2c (version 20100827).
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

#line 1 "MB03TS.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b8 = -1.;
static integer c__2 = 2;
static integer c__4 = 4;
static doublereal c_b21 = 0.;
static logical c_false = FALSE_;
static integer c_n1 = -1;
static integer c__3 = 3;
static doublereal c_b203 = 1.;

/* Subroutine */ int mb03ts_(logical *isham, logical *wantu, integer *n, 
	doublereal *a, integer *lda, doublereal *g, integer *ldg, doublereal *
	u1, integer *ldu1, doublereal *u2, integer *ldu2, integer *j1, 
	integer *n1, integer *n2, doublereal *dwork, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, u1_dim1, u1_offset, u2_dim1, 
	    u2_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static doublereal d__[16]	/* was [4][4] */;
    static integer k;
    static doublereal v[3], x[4]	/* was [2][2] */;
    static integer j2, j3, j4;
    static doublereal v1[3], v2[3], a11, a22, a33;
    static integer nd;
    static doublereal cs, sn, wi1, wi2, wr1, wr2, eps, tau, tau1, tau2;
    static logical lblk;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static integer ierr;
    static doublereal temp;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *), dsyr2_(char 
	    *, integer *, doublereal *, doublereal *, integer *, doublereal *,
	     integer *, doublereal *, integer *, ftnlen), mb01md_(char *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen), 
	    mb01nd_(char *, integer *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), dscal_(
	    integer *, doublereal *, doublereal *, integer *);
    static doublereal scale;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    static doublereal dnorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *), dsymv_(char *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen);
    static doublereal xnorm;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlasy2_(
	    logical *, logical *, integer *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    extern doublereal dlamch_(char *, ftnlen), dlange_(char *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *), dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dlaset_(char *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, ftnlen), dlarfx_(char *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    static doublereal thresh, smlnum;


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

/*     To swap diagonal blocks A11 and A22 of order 1 or 2 in the upper */
/*     quasi-triangular matrix A contained in a skew-Hamiltonian matrix */

/*                   [  A   G  ]          T */
/*             X  =  [       T ],   G = -G, */
/*                   [  0   A  ] */

/*     or in a Hamiltonian matrix */

/*                   [  A   G  ]          T */
/*             X  =  [       T ],   G =  G. */
/*                   [  0  -A  ] */

/*     This routine is a modified version of the LAPACK subroutine */
/*     DLAEX2. */

/*     The matrix A must be in Schur canonical form (as returned by the */
/*     LAPACK routine DHSEQR), that is, block upper triangular with */
/*     1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has */
/*     its diagonal elements equal and its off-diagonal elements of */
/*     opposite sign. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     ISHAM   LOGIGAL */
/*             Specifies the type of X: */
/*             = .TRUE.:   X is a Hamiltonian matrix; */
/*             = .FALSE.:  X is a skew-Hamiltonian matrix. */

/*     WANTU   LOGIGAL */
/*             = .TRUE.:   update the matrices U1 and U2 containing the */
/*                         Schur vectors; */
/*             = .FALSE.:  do not update U1 and U2. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi-triangular matrix A, in Schur */
/*             canonical form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the reordered matrix A, again in Schur canonical form. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper triangular part of the symmetric */
/*             matrix G, if ISHAM = .TRUE., or the strictly upper */
/*             triangular part of the skew-symmetric matrix G, otherwise. */
/*             The rest of this array is not referenced. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the upper or strictly upper triangular part of the */
/*             symmetric or skew-symmetric matrix G, respectively, */
/*             updated by the orthogonal transformation which reorders A. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On entry, if WANTU = .TRUE., the leading N-by-N part of */
/*             this array must contain the matrix U1. */
/*             On exit, if WANTU = .TRUE., the leading N-by-N part of */
/*             this array contains U1, postmultiplied by the orthogonal */
/*             transformation which reorders A. See the description in */
/*             the SLICOT subroutine MB03TD for further details. */
/*             If WANTU = .FALSE., this array is not referenced. */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1. */
/*             LDU1 >= MAX(1,N),  if WANTU = .TRUE.; */
/*             LDU1 >= 1,         otherwise. */

/*     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On entry, if WANTU = .TRUE., the leading N-by-N part of */
/*             this array must contain the matrix U2. */
/*             On exit, if WANTU = .TRUE., the leading N-by-N part of */
/*             this array contains U2, postmultiplied by the orthogonal */
/*             transformation which reorders A. */
/*             If WANTU = .FALSE., this array is not referenced. */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2. */
/*             LDU2 >= MAX(1,N),  if WANTU = .TRUE.; */
/*             LDU2 >= 1,         otherwise. */

/*     J1      (input) INTEGER */
/*             The index of the first row of the first block A11. */
/*             If J1+N1 < N, then A11 is swapped with the block starting */
/*             at (J1+N1+1)-th diagonal element. */
/*             If J1+N1 > N, then A11 is the last block in A and swapped */
/*             with -A11', if ISHAM = .TRUE., */
/*             or    A11', if ISHAM = .FALSE.. */

/*     N1      (input) INTEGER */
/*             The order of the first block A11. N1 = 0, 1 or 2. */

/*     N2      (input) INTEGER */
/*             The order of the second block A22. N2 = 0, 1 or 2. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (N) */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  the transformed matrix A would be too far from Schur */
/*                   form; the blocks are not swapped and A, G, U1 and */
/*                   U2 are unchanged. */

/*     REFERENCES */

/*     [1] Bai, Z., and Demmel, J.W. */
/*        On swapping diagonal blocks in real Schur form. */
/*        Linear Algebra Appl., 186, pp. 73-95, 1993. */

/*     [2] Benner, P., Kressner, D., and Mehrmann, V. */
/*         Skew-Hamiltonian and Hamiltonian Eigenvalue Problems: Theory, */
/*         Algorithms and Applications. Techn. Report, TU Berlin, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHAEX2). */

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

#line 196 "MB03TS.f"
    /* Parameter adjustments */
#line 196 "MB03TS.f"
    a_dim1 = *lda;
#line 196 "MB03TS.f"
    a_offset = 1 + a_dim1;
#line 196 "MB03TS.f"
    a -= a_offset;
#line 196 "MB03TS.f"
    g_dim1 = *ldg;
#line 196 "MB03TS.f"
    g_offset = 1 + g_dim1;
#line 196 "MB03TS.f"
    g -= g_offset;
#line 196 "MB03TS.f"
    u1_dim1 = *ldu1;
#line 196 "MB03TS.f"
    u1_offset = 1 + u1_dim1;
#line 196 "MB03TS.f"
    u1 -= u1_offset;
#line 196 "MB03TS.f"
    u2_dim1 = *ldu2;
#line 196 "MB03TS.f"
    u2_offset = 1 + u2_dim1;
#line 196 "MB03TS.f"
    u2 -= u2_offset;
#line 196 "MB03TS.f"
    --dwork;
#line 196 "MB03TS.f"

#line 196 "MB03TS.f"
    /* Function Body */
#line 196 "MB03TS.f"
    *info = 0;

/*     Quick return if possible. */

#line 200 "MB03TS.f"
    if (*n == 0 || *n1 == 0 || *n2 == 0) {
#line 200 "MB03TS.f"
	return 0;
#line 200 "MB03TS.f"
    }
#line 202 "MB03TS.f"
    lblk = *j1 + *n1 > *n;

#line 204 "MB03TS.f"
    j2 = *j1 + 1;
#line 205 "MB03TS.f"
    j3 = *j1 + 2;
#line 206 "MB03TS.f"
    j4 = *j1 + 3;

#line 208 "MB03TS.f"
    if (lblk && *n1 == 1) {

#line 210 "MB03TS.f"
	if (*isham) {
#line 211 "MB03TS.f"
	    a11 = a[*n + *n * a_dim1];
#line 212 "MB03TS.f"
	    d__1 = a11 * -2.;
#line 212 "MB03TS.f"
	    dlartg_(&g[*n + *n * g_dim1], &d__1, &cs, &sn, &temp);
#line 213 "MB03TS.f"
	    i__1 = *n - 1;
#line 213 "MB03TS.f"
	    drot_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1], &
		    c__1, &cs, &sn);
#line 214 "MB03TS.f"
	    a[*n + *n * a_dim1] = -a11;
#line 215 "MB03TS.f"
	    if (*wantu) {
#line 215 "MB03TS.f"
		drot_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 1], 
			&c__1, &cs, &sn);
#line 215 "MB03TS.f"
	    }
#line 217 "MB03TS.f"
	} else {
#line 218 "MB03TS.f"
	    i__1 = *n - 1;
#line 218 "MB03TS.f"
	    dswap_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1], &
		    c__1);
#line 219 "MB03TS.f"
	    i__1 = *n - 1;
#line 219 "MB03TS.f"
	    dscal_(&i__1, &c_b8, &a[*n * a_dim1 + 1], &c__1);
#line 220 "MB03TS.f"
	    if (*wantu) {
#line 221 "MB03TS.f"
		dswap_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 1],
			 &c__1);
#line 222 "MB03TS.f"
		dscal_(n, &c_b8, &u1[*n * u1_dim1 + 1], &c__1);
#line 223 "MB03TS.f"
	    }
#line 224 "MB03TS.f"
	}

#line 226 "MB03TS.f"
    } else if (lblk && *n1 == 2) {

#line 228 "MB03TS.f"
	if (*isham) {

/*           Reorder Hamiltonian matrix: */

/*                      [ A11  G11  ] */
/*                      [         T ]. */
/*                      [  0  -A11  ] */

#line 236 "MB03TS.f"
	    nd = 4;
#line 237 "MB03TS.f"
	    dlacpy_("Full", &c__2, &c__2, &a[*n - 1 + (*n - 1) * a_dim1], lda,
		     d__, &c__4, (ftnlen)4);
#line 238 "MB03TS.f"
	    dlaset_("All", &c__2, &c__2, &c_b21, &c_b21, &d__[2], &c__4, (
		    ftnlen)3);
#line 239 "MB03TS.f"
	    dlacpy_("Upper", &c__2, &c__2, &g[*n - 1 + (*n - 1) * g_dim1], 
		    ldg, &d__[8], &c__4, (ftnlen)5);
#line 240 "MB03TS.f"
	    d__[9] = d__[12];
#line 241 "MB03TS.f"
	    d__[10] = -d__[0];
#line 242 "MB03TS.f"
	    d__[11] = -d__[4];
#line 243 "MB03TS.f"
	    d__[14] = -d__[1];
#line 244 "MB03TS.f"
	    d__[15] = -d__[5];
#line 245 "MB03TS.f"
	    dnorm = dlange_("Max", &nd, &nd, d__, &c__4, &dwork[1], (ftnlen)3)
		    ;

/*           Compute machine-dependent threshold for test for accepting */
/*           swap. */

#line 250 "MB03TS.f"
	    eps = dlamch_("P", (ftnlen)1);
#line 251 "MB03TS.f"
	    smlnum = dlamch_("S", (ftnlen)1) / eps;
/* Computing MAX */
#line 252 "MB03TS.f"
	    d__1 = eps * 40. * dnorm;
#line 252 "MB03TS.f"
	    thresh = max(d__1,smlnum);

/*           Solve A11*X + X*A11' = scale*G11 for X. */

#line 256 "MB03TS.f"
	    dlasy2_(&c_false, &c_false, &c_n1, &c__2, &c__2, d__, &c__4, &d__[
		    10], &c__4, &d__[8], &c__4, &scale, x, &c__2, &xnorm, &
		    ierr);

/*           Compute symplectic QR decomposition of */

/*                  (  -X11  -X12 ) */
/*                  (  -X21  -X22 ). */
/*                  ( scale    0  ) */
/*                  (    0  scale ) */

#line 266 "MB03TS.f"
	    temp = -x[0];
#line 267 "MB03TS.f"
	    dlartg_(&temp, &scale, v1, v2, x);
#line 268 "MB03TS.f"
	    d__1 = -x[1];
#line 268 "MB03TS.f"
	    dlartg_(x, &d__1, &v1[1], &v2[1], &temp);
#line 269 "MB03TS.f"
	    x[2] = -x[2];
#line 270 "MB03TS.f"
	    x[3] = -x[3];
#line 271 "MB03TS.f"
	    x[0] = 0.;
#line 272 "MB03TS.f"
	    x[1] = scale;
#line 273 "MB03TS.f"
	    drot_(&c__1, &x[2], &c__1, x, &c__1, v1, v2);
#line 274 "MB03TS.f"
	    drot_(&c__1, &x[2], &c__1, &x[3], &c__1, &v1[1], &v2[1]);
#line 275 "MB03TS.f"
	    drot_(&c__1, x, &c__1, &x[1], &c__1, &v1[1], &v2[1]);
#line 276 "MB03TS.f"
	    dlartg_(&x[3], &x[1], &v1[2], &v2[2], &temp);

/*           Perform swap provisionally on D. */

#line 280 "MB03TS.f"
	    drot_(&c__4, d__, &c__4, &d__[2], &c__4, v1, v2);
#line 281 "MB03TS.f"
	    drot_(&c__4, d__, &c__4, &d__[1], &c__4, &v1[1], &v2[1]);
#line 282 "MB03TS.f"
	    drot_(&c__4, &d__[2], &c__4, &d__[3], &c__4, &v1[1], &v2[1]);
#line 283 "MB03TS.f"
	    drot_(&c__4, &d__[1], &c__4, &d__[3], &c__4, &v1[2], &v2[2]);
#line 284 "MB03TS.f"
	    drot_(&c__4, d__, &c__1, &d__[8], &c__1, v1, v2);
#line 285 "MB03TS.f"
	    drot_(&c__4, d__, &c__1, &d__[4], &c__1, &v1[1], &v2[1]);
#line 286 "MB03TS.f"
	    drot_(&c__4, &d__[8], &c__1, &d__[12], &c__1, &v1[1], &v2[1]);
#line 287 "MB03TS.f"
	    drot_(&c__4, &d__[4], &c__1, &d__[12], &c__1, &v1[2], &v2[2]);

/*           Test whether to reject swap. */

/* Computing MAX */
#line 291 "MB03TS.f"
	    d__1 = abs(d__[2]), d__2 = abs(d__[6]), d__1 = max(d__1,d__2), 
		    d__2 = abs(d__[3]), d__1 = max(d__1,d__2), d__2 = abs(d__[
		    7]);
#line 291 "MB03TS.f"
	    if (max(d__1,d__2) > thresh) {
#line 291 "MB03TS.f"
		goto L50;
#line 291 "MB03TS.f"
	    }

#line 294 "MB03TS.f"
	    dlacpy_("All", &c__2, &c__2, d__, &c__4, &a[*n - 1 + (*n - 1) * 
		    a_dim1], lda, (ftnlen)3);
#line 295 "MB03TS.f"
	    dlacpy_("Upper", &c__2, &c__2, &d__[8], &c__4, &g[*n - 1 + (*n - 
		    1) * g_dim1], ldg, (ftnlen)5);

#line 297 "MB03TS.f"
	    if (*n > 2) {
#line 298 "MB03TS.f"
		i__1 = *n - 2;
#line 298 "MB03TS.f"
		drot_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &g[(*n - 1) * 
			g_dim1 + 1], &c__1, v1, v2);
#line 299 "MB03TS.f"
		i__1 = *n - 2;
#line 299 "MB03TS.f"
		drot_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &a[*n * a_dim1 
			+ 1], &c__1, &v1[1], &v2[1]);
#line 300 "MB03TS.f"
		i__1 = *n - 2;
#line 300 "MB03TS.f"
		drot_(&i__1, &g[(*n - 1) * g_dim1 + 1], &c__1, &g[*n * g_dim1 
			+ 1], &c__1, &v1[1], &v2[1]);
#line 301 "MB03TS.f"
		i__1 = *n - 2;
#line 301 "MB03TS.f"
		drot_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1], 
			&c__1, &v1[2], &v2[2]);
#line 302 "MB03TS.f"
	    }

#line 304 "MB03TS.f"
	    if (*wantu) {
#line 305 "MB03TS.f"
		drot_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u2[(*n - 1) * 
			u2_dim1 + 1], &c__1, v1, v2);
#line 306 "MB03TS.f"
		drot_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u1[*n * u1_dim1 
			+ 1], &c__1, &v1[1], &v2[1]);
#line 307 "MB03TS.f"
		drot_(n, &u2[(*n - 1) * u2_dim1 + 1], &c__1, &u2[*n * u2_dim1 
			+ 1], &c__1, &v1[1], &v2[1]);
#line 308 "MB03TS.f"
		drot_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 1], 
			&c__1, &v1[2], &v2[2]);
#line 309 "MB03TS.f"
	    }

#line 311 "MB03TS.f"
	} else {

#line 313 "MB03TS.f"
	    if ((d__1 = a[*n - 1 + *n * a_dim1], abs(d__1)) > (d__2 = a[*n + (
		    *n - 1) * a_dim1], abs(d__2))) {
#line 314 "MB03TS.f"
		temp = g[*n - 1 + *n * g_dim1];
#line 315 "MB03TS.f"
		dlartg_(&temp, &a[*n - 1 + *n * a_dim1], &cs, &sn, &g[*n - 1 
			+ *n * g_dim1]);
#line 316 "MB03TS.f"
		sn = -sn;
#line 317 "MB03TS.f"
		i__1 = *n - 2;
#line 317 "MB03TS.f"
		drot_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1], 
			&c__1, &cs, &sn);

#line 319 "MB03TS.f"
		a[*n - 1 + *n * a_dim1] = -sn * a[*n + (*n - 1) * a_dim1];
#line 320 "MB03TS.f"
		temp = -cs * a[*n + (*n - 1) * a_dim1];
#line 321 "MB03TS.f"
		a[*n + (*n - 1) * a_dim1] = g[*n - 1 + *n * g_dim1];
#line 322 "MB03TS.f"
		g[*n - 1 + *n * g_dim1] = temp;
#line 323 "MB03TS.f"
		if (*wantu) {
#line 323 "MB03TS.f"
		    drot_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 
			    1], &c__1, &cs, &sn);
#line 323 "MB03TS.f"
		}
#line 325 "MB03TS.f"
		i__1 = *n - 2;
#line 325 "MB03TS.f"
		dswap_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &g[(*n - 1) * 
			g_dim1 + 1], &c__1);
#line 326 "MB03TS.f"
		i__1 = *n - 2;
#line 326 "MB03TS.f"
		dscal_(&i__1, &c_b8, &a[(*n - 1) * a_dim1 + 1], &c__1);
#line 327 "MB03TS.f"
		if (*wantu) {
#line 328 "MB03TS.f"
		    dswap_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u2[(*n - 1)
			     * u2_dim1 + 1], &c__1);
#line 329 "MB03TS.f"
		    dscal_(n, &c_b8, &u1[(*n - 1) * u1_dim1 + 1], &c__1);
#line 330 "MB03TS.f"
		}
#line 331 "MB03TS.f"
	    } else {
#line 332 "MB03TS.f"
		temp = g[*n - 1 + *n * g_dim1];
#line 333 "MB03TS.f"
		dlartg_(&temp, &a[*n + (*n - 1) * a_dim1], &cs, &sn, &g[*n - 
			1 + *n * g_dim1]);
#line 334 "MB03TS.f"
		i__1 = *n - 2;
#line 334 "MB03TS.f"
		drot_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &g[(*n - 1) * 
			g_dim1 + 1], &c__1, &cs, &sn);
#line 335 "MB03TS.f"
		a[*n + (*n - 1) * a_dim1] = -sn * a[*n - 1 + *n * a_dim1];
#line 336 "MB03TS.f"
		a[*n - 1 + *n * a_dim1] = cs * a[*n - 1 + *n * a_dim1];
#line 337 "MB03TS.f"
		if (*wantu) {
#line 337 "MB03TS.f"
		    drot_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u2[(*n - 1) 
			    * u2_dim1 + 1], &c__1, &cs, &sn);
#line 337 "MB03TS.f"
		}
#line 339 "MB03TS.f"
		i__1 = *n - 1;
#line 339 "MB03TS.f"
		dswap_(&i__1, &a[*n * a_dim1 + 1], &c__1, &g[*n * g_dim1 + 1],
			 &c__1);
#line 340 "MB03TS.f"
		i__1 = *n - 1;
#line 340 "MB03TS.f"
		dscal_(&i__1, &c_b8, &a[*n * a_dim1 + 1], &c__1);
#line 341 "MB03TS.f"
		if (*wantu) {
#line 342 "MB03TS.f"
		    dswap_(n, &u1[*n * u1_dim1 + 1], &c__1, &u2[*n * u2_dim1 
			    + 1], &c__1);
#line 343 "MB03TS.f"
		    dscal_(n, &c_b8, &u1[*n * u1_dim1 + 1], &c__1);
#line 344 "MB03TS.f"
		}
#line 345 "MB03TS.f"
	    }
#line 346 "MB03TS.f"
	}

/*        Standardize new 2-by-2 block. */

#line 350 "MB03TS.f"
	dlanv2_(&a[*n - 1 + (*n - 1) * a_dim1], &a[*n - 1 + *n * a_dim1], &a[*
		n + (*n - 1) * a_dim1], &a[*n + *n * a_dim1], &wr1, &wi1, &
		wr2, &wi2, &cs, &sn);
#line 352 "MB03TS.f"
	i__1 = *n - 2;
#line 352 "MB03TS.f"
	drot_(&i__1, &a[(*n - 1) * a_dim1 + 1], &c__1, &a[*n * a_dim1 + 1], &
		c__1, &cs, &sn);
#line 353 "MB03TS.f"
	if (*isham) {
#line 354 "MB03TS.f"
	    temp = g[*n - 1 + *n * g_dim1];
#line 355 "MB03TS.f"
	    i__1 = *n - 1;
#line 355 "MB03TS.f"
	    drot_(&i__1, &g[(*n - 1) * g_dim1 + 1], &c__1, &g[*n * g_dim1 + 1]
		    , &c__1, &cs, &sn);
#line 356 "MB03TS.f"
	    tau = cs * temp + sn * g[*n + *n * g_dim1];
#line 357 "MB03TS.f"
	    g[*n + *n * g_dim1] = cs * g[*n + *n * g_dim1] - sn * temp;
#line 358 "MB03TS.f"
	    g[*n - 1 + (*n - 1) * g_dim1] = cs * g[*n - 1 + (*n - 1) * g_dim1]
		     + sn * tau;
#line 359 "MB03TS.f"
	    drot_(&c__1, &g[*n - 1 + *n * g_dim1], ldg, &g[*n + *n * g_dim1], 
		    ldg, &cs, &sn);
#line 360 "MB03TS.f"
	} else {
#line 361 "MB03TS.f"
	    i__1 = *n - 2;
#line 361 "MB03TS.f"
	    drot_(&i__1, &g[(*n - 1) * g_dim1 + 1], &c__1, &g[*n * g_dim1 + 1]
		    , &c__1, &cs, &sn);
#line 362 "MB03TS.f"
	}
#line 363 "MB03TS.f"
	if (*wantu) {
#line 364 "MB03TS.f"
	    drot_(n, &u1[(*n - 1) * u1_dim1 + 1], &c__1, &u1[*n * u1_dim1 + 1]
		    , &c__1, &cs, &sn);
#line 365 "MB03TS.f"
	    drot_(n, &u2[(*n - 1) * u2_dim1 + 1], &c__1, &u2[*n * u2_dim1 + 1]
		    , &c__1, &cs, &sn);
#line 366 "MB03TS.f"
	}

#line 368 "MB03TS.f"
    } else if (*n1 == 1 && *n2 == 1) {

/*        Swap two 1-by-1 blocks. */

#line 372 "MB03TS.f"
	a11 = a[*j1 + *j1 * a_dim1];
#line 373 "MB03TS.f"
	a22 = a[j2 + j2 * a_dim1];

/*        Determine the transformation to perform the interchange. */

#line 377 "MB03TS.f"
	d__1 = a22 - a11;
#line 377 "MB03TS.f"
	dlartg_(&a[*j1 + j2 * a_dim1], &d__1, &cs, &sn, &temp);

/*        Apply transformation to the matrix A. */

#line 381 "MB03TS.f"
	if (j3 <= *n) {
#line 381 "MB03TS.f"
	    i__1 = *n - *j1 - 1;
#line 381 "MB03TS.f"
	    drot_(&i__1, &a[*j1 + j3 * a_dim1], lda, &a[j2 + j3 * a_dim1], 
		    lda, &cs, &sn);
#line 381 "MB03TS.f"
	}
#line 383 "MB03TS.f"
	i__1 = *j1 - 1;
#line 383 "MB03TS.f"
	drot_(&i__1, &a[*j1 * a_dim1 + 1], &c__1, &a[j2 * a_dim1 + 1], &c__1, 
		&cs, &sn);

#line 385 "MB03TS.f"
	a[*j1 + *j1 * a_dim1] = a22;
#line 386 "MB03TS.f"
	a[j2 + j2 * a_dim1] = a11;

/*        Apply transformation to the matrix G. */

#line 390 "MB03TS.f"
	if (*isham) {
#line 391 "MB03TS.f"
	    temp = g[*j1 + j2 * g_dim1];
#line 392 "MB03TS.f"
	    drot_(j1, &g[*j1 * g_dim1 + 1], &c__1, &g[j2 * g_dim1 + 1], &c__1,
		     &cs, &sn);
#line 393 "MB03TS.f"
	    tau = cs * temp + sn * g[j2 + j2 * g_dim1];
#line 394 "MB03TS.f"
	    g[j2 + j2 * g_dim1] = cs * g[j2 + j2 * g_dim1] - sn * temp;
#line 395 "MB03TS.f"
	    g[*j1 + *j1 * g_dim1] = cs * g[*j1 + *j1 * g_dim1] + sn * tau;
#line 396 "MB03TS.f"
	    i__1 = *n - *j1;
#line 396 "MB03TS.f"
	    drot_(&i__1, &g[*j1 + j2 * g_dim1], ldg, &g[j2 + j2 * g_dim1], 
		    ldg, &cs, &sn);
#line 397 "MB03TS.f"
	} else {
#line 398 "MB03TS.f"
	    if (*n > *j1 + 1) {
#line 398 "MB03TS.f"
		i__1 = *n - *j1 - 1;
#line 398 "MB03TS.f"
		drot_(&i__1, &g[*j1 + (*j1 + 2) * g_dim1], ldg, &g[j2 + (*j1 
			+ 2) * g_dim1], ldg, &cs, &sn);
#line 398 "MB03TS.f"
	    }
#line 401 "MB03TS.f"
	    i__1 = *j1 - 1;
#line 401 "MB03TS.f"
	    drot_(&i__1, &g[*j1 * g_dim1 + 1], &c__1, &g[j2 * g_dim1 + 1], &
		    c__1, &cs, &sn);
#line 402 "MB03TS.f"
	}
#line 403 "MB03TS.f"
	if (*wantu) {

/*           Accumulate transformation in the matrices U1 and U2. */

#line 407 "MB03TS.f"
	    drot_(n, &u1[*j1 * u1_dim1 + 1], &c__1, &u1[j2 * u1_dim1 + 1], &
		    c__1, &cs, &sn);
#line 408 "MB03TS.f"
	    drot_(n, &u2[*j1 * u2_dim1 + 1], &c__1, &u2[j2 * u2_dim1 + 1], &
		    c__1, &cs, &sn);
#line 409 "MB03TS.f"
	}

#line 411 "MB03TS.f"
    } else {

/*        Swapping involves at least one 2-by-2 block. */

/*        Copy the diagonal block of order N1+N2 to the local array D */
/*        and compute its norm. */

#line 418 "MB03TS.f"
	nd = *n1 + *n2;
#line 419 "MB03TS.f"
	dlacpy_("Full", &nd, &nd, &a[*j1 + *j1 * a_dim1], lda, d__, &c__4, (
		ftnlen)4);
#line 420 "MB03TS.f"
	dnorm = dlange_("Max", &nd, &nd, d__, &c__4, &dwork[1], (ftnlen)3);

/*        Compute machine-dependent threshold for test for accepting */
/*        swap. */

#line 425 "MB03TS.f"
	eps = dlamch_("P", (ftnlen)1);
#line 426 "MB03TS.f"
	smlnum = dlamch_("S", (ftnlen)1) / eps;
/* Computing MAX */
#line 427 "MB03TS.f"
	d__1 = eps * 30. * dnorm;
#line 427 "MB03TS.f"
	thresh = max(d__1,smlnum);

/*        Solve A11*X - X*A22 = scale*A12 for X. */

#line 431 "MB03TS.f"
	dlasy2_(&c_false, &c_false, &c_n1, n1, n2, d__, &c__4, &d__[*n1 + 1 + 
		(*n1 + 1 << 2) - 5], &c__4, &d__[(*n1 + 1 << 2) - 4], &c__4, &
		scale, x, &c__2, &xnorm, &ierr);

/*        Swap the adjacent diagonal blocks. */

#line 437 "MB03TS.f"
	k = *n1 + *n1 + *n2 - 3;
#line 438 "MB03TS.f"
	switch (k) {
#line 438 "MB03TS.f"
	    case 1:  goto L10;
#line 438 "MB03TS.f"
	    case 2:  goto L20;
#line 438 "MB03TS.f"
	    case 3:  goto L30;
#line 438 "MB03TS.f"
	}

#line 440 "MB03TS.f"
L10:

/*        N1 = 1, N2 = 2: generate elementary reflector H so that: */

/*        ( scale, X11, X12 ) H = ( 0, 0, * ). */

#line 446 "MB03TS.f"
	v[0] = scale;
#line 447 "MB03TS.f"
	v[1] = x[0];
#line 448 "MB03TS.f"
	v[2] = x[2];
#line 449 "MB03TS.f"
	dlarfg_(&c__3, &v[2], v, &c__1, &tau);
#line 450 "MB03TS.f"
	v[2] = 1.;
#line 451 "MB03TS.f"
	a11 = a[*j1 + *j1 * a_dim1];

/*        Perform swap provisionally on diagonal block in D. */

#line 455 "MB03TS.f"
	dlarfx_("Left", &c__3, &c__3, v, &tau, d__, &c__4, &dwork[1], (ftnlen)
		4);
#line 456 "MB03TS.f"
	dlarfx_("Right", &c__3, &c__3, v, &tau, d__, &c__4, &dwork[1], (
		ftnlen)5);

/*        Test whether to reject swap. */

/* Computing MAX */
#line 460 "MB03TS.f"
	d__2 = abs(d__[2]), d__3 = abs(d__[6]), d__2 = max(d__2,d__3), d__3 = 
		(d__1 = d__[10] - a11, abs(d__1));
#line 460 "MB03TS.f"
	if (max(d__2,d__3) > thresh) {
#line 460 "MB03TS.f"
	    goto L50;
#line 460 "MB03TS.f"
	}

/*        Accept swap: apply transformation to the entire matrix A. */

#line 465 "MB03TS.f"
	i__1 = *n - *j1 + 1;
#line 465 "MB03TS.f"
	dlarfx_("Left", &c__3, &i__1, v, &tau, &a[*j1 + *j1 * a_dim1], lda, &
		dwork[1], (ftnlen)4);
#line 466 "MB03TS.f"
	dlarfx_("Right", &j2, &c__3, v, &tau, &a[*j1 * a_dim1 + 1], lda, &
		dwork[1], (ftnlen)5);

#line 468 "MB03TS.f"
	a[j3 + *j1 * a_dim1] = 0.;
#line 469 "MB03TS.f"
	a[j3 + j2 * a_dim1] = 0.;
#line 470 "MB03TS.f"
	a[j3 + j3 * a_dim1] = a11;

/*        Apply transformation to G. */

#line 474 "MB03TS.f"
	if (*isham) {
#line 475 "MB03TS.f"
	    i__1 = *j1 - 1;
#line 475 "MB03TS.f"
	    dlarfx_("Right", &i__1, &c__3, v, &tau, &g[*j1 * g_dim1 + 1], ldg,
		     &dwork[1], (ftnlen)5);
#line 476 "MB03TS.f"
	    dsymv_("Upper", &c__3, &tau, &g[*j1 + *j1 * g_dim1], ldg, v, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
#line 478 "MB03TS.f"
	    temp = tau * -.5 * ddot_(&c__3, &dwork[1], &c__1, v, &c__1);
#line 479 "MB03TS.f"
	    daxpy_(&c__3, &temp, v, &c__1, &dwork[1], &c__1);
#line 480 "MB03TS.f"
	    dsyr2_("Upper", &c__3, &c_b8, v, &c__1, &dwork[1], &c__1, &g[*j1 
		    + *j1 * g_dim1], ldg, (ftnlen)5);
#line 482 "MB03TS.f"
	    if (*n > *j1 + 2) {
#line 482 "MB03TS.f"
		i__1 = *n - *j1 - 2;
#line 482 "MB03TS.f"
		dlarfx_("Left", &c__3, &i__1, v, &tau, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 482 "MB03TS.f"
	    }
#line 485 "MB03TS.f"
	} else {
#line 486 "MB03TS.f"
	    i__1 = *j1 - 1;
#line 486 "MB03TS.f"
	    dlarfx_("Right", &i__1, &c__3, v, &tau, &g[*j1 * g_dim1 + 1], ldg,
		     &dwork[1], (ftnlen)5);
#line 487 "MB03TS.f"
	    mb01md_("Upper", &c__3, &tau, &g[*j1 + *j1 * g_dim1], ldg, v, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
#line 489 "MB03TS.f"
	    mb01nd_("Upper", &c__3, &c_b203, v, &c__1, &dwork[1], &c__1, &g[*
		    j1 + *j1 * g_dim1], ldg, (ftnlen)5);
#line 491 "MB03TS.f"
	    if (*n > *j1 + 2) {
#line 491 "MB03TS.f"
		i__1 = *n - *j1 - 2;
#line 491 "MB03TS.f"
		dlarfx_("Left", &c__3, &i__1, v, &tau, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 491 "MB03TS.f"
	    }
#line 494 "MB03TS.f"
	}

#line 496 "MB03TS.f"
	if (*wantu) {

/*           Accumulate transformation in the matrices U1 and U2. */

#line 500 "MB03TS.f"
	    dlarfx_("R", n, &c__3, v, &tau, &u1[*j1 * u1_dim1 + 1], ldu1, &
		    dwork[1], (ftnlen)1);
#line 501 "MB03TS.f"
	    dlarfx_("R", n, &c__3, v, &tau, &u2[*j1 * u2_dim1 + 1], ldu2, &
		    dwork[1], (ftnlen)1);
#line 502 "MB03TS.f"
	}
#line 503 "MB03TS.f"
	goto L40;

#line 505 "MB03TS.f"
L20:

/*        N1 = 2, N2 = 1: generate elementary reflector H so that: */

/*        H (  -X11 ) = ( * ) */
/*          (  -X21 ) = ( 0 ). */
/*          ( scale ) = ( 0 ) */

#line 513 "MB03TS.f"
	v[0] = -x[0];
#line 514 "MB03TS.f"
	v[1] = -x[1];
#line 515 "MB03TS.f"
	v[2] = scale;
#line 516 "MB03TS.f"
	dlarfg_(&c__3, v, &v[1], &c__1, &tau);
#line 517 "MB03TS.f"
	v[0] = 1.;
#line 518 "MB03TS.f"
	a33 = a[j3 + j3 * a_dim1];

/*        Perform swap provisionally on diagonal block in D. */

#line 522 "MB03TS.f"
	dlarfx_("L", &c__3, &c__3, v, &tau, d__, &c__4, &dwork[1], (ftnlen)1);
#line 523 "MB03TS.f"
	dlarfx_("R", &c__3, &c__3, v, &tau, d__, &c__4, &dwork[1], (ftnlen)1);

/*        Test whether to reject swap. */

/* Computing MAX */
#line 527 "MB03TS.f"
	d__2 = abs(d__[1]), d__3 = abs(d__[2]), d__2 = max(d__2,d__3), d__3 = 
		(d__1 = d__[0] - a33, abs(d__1));
#line 527 "MB03TS.f"
	if (max(d__2,d__3) > thresh) {
#line 527 "MB03TS.f"
	    goto L50;
#line 527 "MB03TS.f"
	}

/*        Accept swap: apply transformation to the entire matrix A. */

#line 532 "MB03TS.f"
	dlarfx_("Right", &j3, &c__3, v, &tau, &a[*j1 * a_dim1 + 1], lda, &
		dwork[1], (ftnlen)5);
#line 533 "MB03TS.f"
	i__1 = *n - *j1;
#line 533 "MB03TS.f"
	dlarfx_("Left", &c__3, &i__1, v, &tau, &a[*j1 + j2 * a_dim1], lda, &
		dwork[1], (ftnlen)4);

#line 535 "MB03TS.f"
	a[*j1 + *j1 * a_dim1] = a33;
#line 536 "MB03TS.f"
	a[j2 + *j1 * a_dim1] = 0.;
#line 537 "MB03TS.f"
	a[j3 + *j1 * a_dim1] = 0.;

/*        Apply transformation to G. */

#line 541 "MB03TS.f"
	if (*isham) {
#line 542 "MB03TS.f"
	    i__1 = *j1 - 1;
#line 542 "MB03TS.f"
	    dlarfx_("Right", &i__1, &c__3, v, &tau, &g[*j1 * g_dim1 + 1], ldg,
		     &dwork[1], (ftnlen)5);
#line 543 "MB03TS.f"
	    dsymv_("Upper", &c__3, &tau, &g[*j1 + *j1 * g_dim1], ldg, v, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
#line 545 "MB03TS.f"
	    temp = tau * -.5 * ddot_(&c__3, &dwork[1], &c__1, v, &c__1);
#line 546 "MB03TS.f"
	    daxpy_(&c__3, &temp, v, &c__1, &dwork[1], &c__1);
#line 547 "MB03TS.f"
	    dsyr2_("Upper", &c__3, &c_b8, v, &c__1, &dwork[1], &c__1, &g[*j1 
		    + *j1 * g_dim1], ldg, (ftnlen)5);
#line 549 "MB03TS.f"
	    if (*n > *j1 + 2) {
#line 549 "MB03TS.f"
		i__1 = *n - *j1 - 2;
#line 549 "MB03TS.f"
		dlarfx_("Left", &c__3, &i__1, v, &tau, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 549 "MB03TS.f"
	    }
#line 552 "MB03TS.f"
	} else {
#line 553 "MB03TS.f"
	    i__1 = *j1 - 1;
#line 553 "MB03TS.f"
	    dlarfx_("Right", &i__1, &c__3, v, &tau, &g[*j1 * g_dim1 + 1], ldg,
		     &dwork[1], (ftnlen)5);
#line 554 "MB03TS.f"
	    mb01md_("Upper", &c__3, &tau, &g[*j1 + *j1 * g_dim1], ldg, v, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
#line 556 "MB03TS.f"
	    mb01nd_("Upper", &c__3, &c_b203, v, &c__1, &dwork[1], &c__1, &g[*
		    j1 + *j1 * g_dim1], ldg, (ftnlen)5);
#line 558 "MB03TS.f"
	    if (*n > *j1 + 2) {
#line 558 "MB03TS.f"
		i__1 = *n - *j1 - 2;
#line 558 "MB03TS.f"
		dlarfx_("Left", &c__3, &i__1, v, &tau, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 558 "MB03TS.f"
	    }
#line 561 "MB03TS.f"
	}

#line 563 "MB03TS.f"
	if (*wantu) {

/*           Accumulate transformation in the matrices U1 and U2. */

#line 567 "MB03TS.f"
	    dlarfx_("R", n, &c__3, v, &tau, &u1[*j1 * u1_dim1 + 1], ldu1, &
		    dwork[1], (ftnlen)1);
#line 568 "MB03TS.f"
	    dlarfx_("R", n, &c__3, v, &tau, &u2[*j1 * u2_dim1 + 1], ldu2, &
		    dwork[1], (ftnlen)1);
#line 569 "MB03TS.f"
	}
#line 570 "MB03TS.f"
	goto L40;

#line 572 "MB03TS.f"
L30:

/*        N1 = 2, N2 = 2: generate elementary reflectors H(1) and H(2) so */
/*        that: */

/*        H(2) H(1) (  -X11  -X12 ) = (  *  * ) */
/*                  (  -X21  -X22 )   (  0  * ). */
/*                  ( scale    0  )   (  0  0 ) */
/*                  (    0  scale )   (  0  0 ) */

#line 582 "MB03TS.f"
	v1[0] = -x[0];
#line 583 "MB03TS.f"
	v1[1] = -x[1];
#line 584 "MB03TS.f"
	v1[2] = scale;
#line 585 "MB03TS.f"
	dlarfg_(&c__3, v1, &v1[1], &c__1, &tau1);
#line 586 "MB03TS.f"
	v1[0] = 1.;

#line 588 "MB03TS.f"
	temp = -tau1 * (x[2] + v1[1] * x[3]);
#line 589 "MB03TS.f"
	v2[0] = -temp * v1[1] - x[3];
#line 590 "MB03TS.f"
	v2[1] = -temp * v1[2];
#line 591 "MB03TS.f"
	v2[2] = scale;
#line 592 "MB03TS.f"
	dlarfg_(&c__3, v2, &v2[1], &c__1, &tau2);
#line 593 "MB03TS.f"
	v2[0] = 1.;

/*        Perform swap provisionally on diagonal block in D. */

#line 597 "MB03TS.f"
	dlarfx_("L", &c__3, &c__4, v1, &tau1, d__, &c__4, &dwork[1], (ftnlen)
		1);
#line 598 "MB03TS.f"
	dlarfx_("R", &c__4, &c__3, v1, &tau1, d__, &c__4, &dwork[1], (ftnlen)
		1);
#line 599 "MB03TS.f"
	dlarfx_("L", &c__3, &c__4, v2, &tau2, &d__[1], &c__4, &dwork[1], (
		ftnlen)1);
#line 600 "MB03TS.f"
	dlarfx_("R", &c__4, &c__3, v2, &tau2, &d__[4], &c__4, &dwork[1], (
		ftnlen)1);

/*        Test whether to reject swap. */

/* Computing MAX */
#line 604 "MB03TS.f"
	d__1 = abs(d__[2]), d__2 = abs(d__[6]), d__1 = max(d__1,d__2), d__2 = 
		abs(d__[3]), d__1 = max(d__1,d__2), d__2 = abs(d__[7]);
#line 604 "MB03TS.f"
	if (max(d__1,d__2) > thresh) {
#line 604 "MB03TS.f"
	    goto L50;
#line 604 "MB03TS.f"
	}

/*        Accept swap: apply transformation to the entire matrix A. */

#line 609 "MB03TS.f"
	i__1 = *n - *j1 + 1;
#line 609 "MB03TS.f"
	dlarfx_("L", &c__3, &i__1, v1, &tau1, &a[*j1 + *j1 * a_dim1], lda, &
		dwork[1], (ftnlen)1);
#line 610 "MB03TS.f"
	dlarfx_("R", &j4, &c__3, v1, &tau1, &a[*j1 * a_dim1 + 1], lda, &dwork[
		1], (ftnlen)1);
#line 611 "MB03TS.f"
	i__1 = *n - *j1 + 1;
#line 611 "MB03TS.f"
	dlarfx_("L", &c__3, &i__1, v2, &tau2, &a[j2 + *j1 * a_dim1], lda, &
		dwork[1], (ftnlen)1);
#line 612 "MB03TS.f"
	dlarfx_("R", &j4, &c__3, v2, &tau2, &a[j2 * a_dim1 + 1], lda, &dwork[
		1], (ftnlen)1);

#line 614 "MB03TS.f"
	a[j3 + *j1 * a_dim1] = 0.;
#line 615 "MB03TS.f"
	a[j3 + j2 * a_dim1] = 0.;
#line 616 "MB03TS.f"
	a[j4 + *j1 * a_dim1] = 0.;
#line 617 "MB03TS.f"
	a[j4 + j2 * a_dim1] = 0.;

/*        Apply transformation to G. */

#line 621 "MB03TS.f"
	if (*isham) {
#line 622 "MB03TS.f"
	    i__1 = *j1 - 1;
#line 622 "MB03TS.f"
	    dlarfx_("Right", &i__1, &c__3, v1, &tau1, &g[*j1 * g_dim1 + 1], 
		    ldg, &dwork[1], (ftnlen)5);
#line 624 "MB03TS.f"
	    dsymv_("Upper", &c__3, &tau1, &g[*j1 + *j1 * g_dim1], ldg, v1, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
#line 626 "MB03TS.f"
	    temp = tau1 * -.5 * ddot_(&c__3, &dwork[1], &c__1, v1, &c__1);
#line 627 "MB03TS.f"
	    daxpy_(&c__3, &temp, v1, &c__1, &dwork[1], &c__1);
#line 628 "MB03TS.f"
	    dsyr2_("Upper", &c__3, &c_b8, v1, &c__1, &dwork[1], &c__1, &g[*j1 
		    + *j1 * g_dim1], ldg, (ftnlen)5);
#line 630 "MB03TS.f"
	    if (*n > *j1 + 2) {
#line 630 "MB03TS.f"
		i__1 = *n - *j1 - 2;
#line 630 "MB03TS.f"
		dlarfx_("Left", &c__3, &i__1, v1, &tau1, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 630 "MB03TS.f"
	    }

#line 634 "MB03TS.f"
	    i__1 = j2 - 1;
#line 634 "MB03TS.f"
	    dlarfx_("Right", &i__1, &c__3, v2, &tau2, &g[j2 * g_dim1 + 1], 
		    ldg, &dwork[1], (ftnlen)5);
#line 636 "MB03TS.f"
	    dsymv_("Upper", &c__3, &tau2, &g[j2 + j2 * g_dim1], ldg, v2, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
#line 638 "MB03TS.f"
	    temp = tau2 * -.5 * ddot_(&c__3, &dwork[1], &c__1, v2, &c__1);
#line 639 "MB03TS.f"
	    daxpy_(&c__3, &temp, v2, &c__1, &dwork[1], &c__1);
#line 640 "MB03TS.f"
	    dsyr2_("Upper", &c__3, &c_b8, v2, &c__1, &dwork[1], &c__1, &g[j2 
		    + j2 * g_dim1], ldg, (ftnlen)5);
#line 642 "MB03TS.f"
	    if (*n > j2 + 2) {
#line 642 "MB03TS.f"
		i__1 = *n - j2 - 2;
#line 642 "MB03TS.f"
		dlarfx_("Left", &c__3, &i__1, v2, &tau2, &g[j2 + (j2 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 642 "MB03TS.f"
	    }
#line 645 "MB03TS.f"
	} else {
#line 646 "MB03TS.f"
	    i__1 = *j1 - 1;
#line 646 "MB03TS.f"
	    dlarfx_("Right", &i__1, &c__3, v1, &tau1, &g[*j1 * g_dim1 + 1], 
		    ldg, &dwork[1], (ftnlen)5);
#line 648 "MB03TS.f"
	    mb01md_("Upper", &c__3, &tau1, &g[*j1 + *j1 * g_dim1], ldg, v1, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
#line 650 "MB03TS.f"
	    mb01nd_("Upper", &c__3, &c_b203, v1, &c__1, &dwork[1], &c__1, &g[*
		    j1 + *j1 * g_dim1], ldg, (ftnlen)5);
#line 652 "MB03TS.f"
	    if (*n > *j1 + 2) {
#line 652 "MB03TS.f"
		i__1 = *n - *j1 - 2;
#line 652 "MB03TS.f"
		dlarfx_("Left", &c__3, &i__1, v1, &tau1, &g[*j1 + (*j1 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 652 "MB03TS.f"
	    }
#line 655 "MB03TS.f"
	    i__1 = j2 - 1;
#line 655 "MB03TS.f"
	    dlarfx_("Right", &i__1, &c__3, v2, &tau2, &g[j2 * g_dim1 + 1], 
		    ldg, &dwork[1], (ftnlen)5);
#line 657 "MB03TS.f"
	    mb01md_("Upper", &c__3, &tau2, &g[j2 + j2 * g_dim1], ldg, v2, &
		    c__1, &c_b21, &dwork[1], &c__1, (ftnlen)5);
#line 659 "MB03TS.f"
	    mb01nd_("Upper", &c__3, &c_b203, v2, &c__1, &dwork[1], &c__1, &g[
		    j2 + j2 * g_dim1], ldg, (ftnlen)5);
#line 661 "MB03TS.f"
	    if (*n > j2 + 2) {
#line 661 "MB03TS.f"
		i__1 = *n - j2 - 2;
#line 661 "MB03TS.f"
		dlarfx_("Left", &c__3, &i__1, v2, &tau2, &g[j2 + (j2 + 3) * 
			g_dim1], ldg, &dwork[1], (ftnlen)4);
#line 661 "MB03TS.f"
	    }
#line 664 "MB03TS.f"
	}

#line 666 "MB03TS.f"
	if (*wantu) {

/*           Accumulate transformation in the matrices U1 and U2. */

#line 670 "MB03TS.f"
	    dlarfx_("R", n, &c__3, v1, &tau1, &u1[*j1 * u1_dim1 + 1], ldu1, &
		    dwork[1], (ftnlen)1);
#line 671 "MB03TS.f"
	    dlarfx_("R", n, &c__3, v2, &tau2, &u1[j2 * u1_dim1 + 1], ldu1, &
		    dwork[1], (ftnlen)1);
#line 672 "MB03TS.f"
	    dlarfx_("R", n, &c__3, v1, &tau1, &u2[*j1 * u2_dim1 + 1], ldu2, &
		    dwork[1], (ftnlen)1);
#line 673 "MB03TS.f"
	    dlarfx_("R", n, &c__3, v2, &tau2, &u2[j2 * u2_dim1 + 1], ldu2, &
		    dwork[1], (ftnlen)1);
#line 674 "MB03TS.f"
	}

#line 676 "MB03TS.f"
L40:

#line 678 "MB03TS.f"
	if (*n2 == 2) {

/*           Standardize new 2-by-2 block A11. */

#line 682 "MB03TS.f"
	    dlanv2_(&a[*j1 + *j1 * a_dim1], &a[*j1 + j2 * a_dim1], &a[j2 + *
		    j1 * a_dim1], &a[j2 + j2 * a_dim1], &wr1, &wi1, &wr2, &
		    wi2, &cs, &sn);
#line 684 "MB03TS.f"
	    i__1 = *n - *j1 - 1;
#line 684 "MB03TS.f"
	    drot_(&i__1, &a[*j1 + (*j1 + 2) * a_dim1], lda, &a[j2 + (*j1 + 2) 
		    * a_dim1], lda, &cs, &sn);
#line 686 "MB03TS.f"
	    i__1 = *j1 - 1;
#line 686 "MB03TS.f"
	    drot_(&i__1, &a[*j1 * a_dim1 + 1], &c__1, &a[j2 * a_dim1 + 1], &
		    c__1, &cs, &sn);
#line 687 "MB03TS.f"
	    if (*isham) {
#line 688 "MB03TS.f"
		temp = g[*j1 + j2 * g_dim1];
#line 689 "MB03TS.f"
		drot_(j1, &g[*j1 * g_dim1 + 1], &c__1, &g[j2 * g_dim1 + 1], &
			c__1, &cs, &sn);
#line 690 "MB03TS.f"
		tau = cs * temp + sn * g[j2 + j2 * g_dim1];
#line 691 "MB03TS.f"
		g[j2 + j2 * g_dim1] = cs * g[j2 + j2 * g_dim1] - sn * temp;
#line 692 "MB03TS.f"
		g[*j1 + *j1 * g_dim1] = cs * g[*j1 + *j1 * g_dim1] + sn * tau;
#line 693 "MB03TS.f"
		i__1 = *n - *j1;
#line 693 "MB03TS.f"
		drot_(&i__1, &g[*j1 + j2 * g_dim1], ldg, &g[j2 + j2 * g_dim1],
			 ldg, &cs, &sn);
#line 694 "MB03TS.f"
	    } else {
#line 695 "MB03TS.f"
		if (*n > *j1 + 1) {
#line 695 "MB03TS.f"
		    i__1 = *n - *j1 - 1;
#line 695 "MB03TS.f"
		    drot_(&i__1, &g[*j1 + (*j1 + 2) * g_dim1], ldg, &g[j2 + (*
			    j1 + 2) * g_dim1], ldg, &cs, &sn);
#line 695 "MB03TS.f"
		}
#line 698 "MB03TS.f"
		i__1 = *j1 - 1;
#line 698 "MB03TS.f"
		drot_(&i__1, &g[*j1 * g_dim1 + 1], &c__1, &g[j2 * g_dim1 + 1],
			 &c__1, &cs, &sn);
#line 699 "MB03TS.f"
	    }
#line 700 "MB03TS.f"
	    if (*wantu) {
#line 701 "MB03TS.f"
		drot_(n, &u1[*j1 * u1_dim1 + 1], &c__1, &u1[j2 * u1_dim1 + 1],
			 &c__1, &cs, &sn);
#line 702 "MB03TS.f"
		drot_(n, &u2[*j1 * u2_dim1 + 1], &c__1, &u2[j2 * u2_dim1 + 1],
			 &c__1, &cs, &sn);
#line 703 "MB03TS.f"
	    }
#line 704 "MB03TS.f"
	}

#line 706 "MB03TS.f"
	if (*n1 == 2) {

/*           Standardize new 2-by-2 block A22. */

#line 710 "MB03TS.f"
	    j3 = *j1 + *n2;
#line 711 "MB03TS.f"
	    j4 = j3 + 1;
#line 712 "MB03TS.f"
	    dlanv2_(&a[j3 + j3 * a_dim1], &a[j3 + j4 * a_dim1], &a[j4 + j3 * 
		    a_dim1], &a[j4 + j4 * a_dim1], &wr1, &wi1, &wr2, &wi2, &
		    cs, &sn);
#line 714 "MB03TS.f"
	    if (j3 + 2 <= *n) {
#line 714 "MB03TS.f"
		i__1 = *n - j3 - 1;
#line 714 "MB03TS.f"
		drot_(&i__1, &a[j3 + (j3 + 2) * a_dim1], lda, &a[j4 + (j3 + 2)
			 * a_dim1], lda, &cs, &sn);
#line 714 "MB03TS.f"
	    }
#line 717 "MB03TS.f"
	    i__1 = j3 - 1;
#line 717 "MB03TS.f"
	    drot_(&i__1, &a[j3 * a_dim1 + 1], &c__1, &a[j4 * a_dim1 + 1], &
		    c__1, &cs, &sn);
#line 718 "MB03TS.f"
	    if (*isham) {
#line 719 "MB03TS.f"
		temp = g[j3 + j4 * g_dim1];
#line 720 "MB03TS.f"
		drot_(&j3, &g[j3 * g_dim1 + 1], &c__1, &g[j4 * g_dim1 + 1], &
			c__1, &cs, &sn);
#line 721 "MB03TS.f"
		tau = cs * temp + sn * g[j4 + j4 * g_dim1];
#line 722 "MB03TS.f"
		g[j4 + j4 * g_dim1] = cs * g[j4 + j4 * g_dim1] - sn * temp;
#line 723 "MB03TS.f"
		g[j3 + j3 * g_dim1] = cs * g[j3 + j3 * g_dim1] + sn * tau;
#line 724 "MB03TS.f"
		i__1 = *n - j3;
#line 724 "MB03TS.f"
		drot_(&i__1, &g[j3 + j4 * g_dim1], ldg, &g[j4 + j4 * g_dim1], 
			ldg, &cs, &sn);
#line 725 "MB03TS.f"
	    } else {
#line 726 "MB03TS.f"
		if (*n > j3 + 1) {
#line 726 "MB03TS.f"
		    i__1 = *n - j3 - 1;
#line 726 "MB03TS.f"
		    drot_(&i__1, &g[j3 + (j3 + 2) * g_dim1], ldg, &g[j4 + (j3 
			    + 2) * g_dim1], ldg, &cs, &sn);
#line 726 "MB03TS.f"
		}
#line 729 "MB03TS.f"
		i__1 = j3 - 1;
#line 729 "MB03TS.f"
		drot_(&i__1, &g[j3 * g_dim1 + 1], &c__1, &g[j4 * g_dim1 + 1], 
			&c__1, &cs, &sn);
#line 730 "MB03TS.f"
	    }
#line 731 "MB03TS.f"
	    if (*wantu) {
#line 732 "MB03TS.f"
		drot_(n, &u1[j3 * u1_dim1 + 1], &c__1, &u1[j4 * u1_dim1 + 1], 
			&c__1, &cs, &sn);
#line 733 "MB03TS.f"
		drot_(n, &u2[j3 * u2_dim1 + 1], &c__1, &u2[j4 * u2_dim1 + 1], 
			&c__1, &cs, &sn);
#line 734 "MB03TS.f"
	    }
#line 735 "MB03TS.f"
	}

#line 737 "MB03TS.f"
    }
#line 738 "MB03TS.f"
    return 0;

/*     Exit with INFO = 1 if swap was rejected. */

#line 742 "MB03TS.f"
L50:
#line 743 "MB03TS.f"
    *info = 1;
#line 744 "MB03TS.f"
    return 0;
/* *** Last line of MB03TS *** */
} /* mb03ts_ */

