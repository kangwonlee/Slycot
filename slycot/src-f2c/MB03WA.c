#line 1 "MB03WA.f"
/* MB03WA.f -- translated by f2c (version 20100827).
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

#line 1 "MB03WA.f"
/* Table of constant values */

static integer c__4 = 4;
static doublereal c_b5 = 0.;
static integer c__1 = 1;
static integer c__2 = 2;
static doublereal c_b42 = 1.;
static doublereal c_b48 = -1.;

/* Subroutine */ int mb03wa_(logical *wantq, logical *wantz, integer *n1, 
	integer *n2, doublereal *a, integer *lda, doublereal *b, integer *ldb,
	 doublereal *q, integer *ldq, doublereal *z__, integer *ldz, integer *
	info)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, q_dim1, q_offset, z_dim1, 
	    z_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static doublereal f, g;
    static integer i__, m;
    static doublereal s[16]	/* was [4][4] */, t[16]	/* was [4][4] */, be[
	    2], ai[2], ar[2], sa, sb, li[16]	/* was [4][4] */, ir[16]	
	    /* was [4][4] */, ss, ws, eps;
    static logical weak;
    static doublereal ddum, taul[4], dsum;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static doublereal taur[4], scpy[16]	/* was [4][4] */, tcpy[16]	/* 
	    was [4][4] */;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    static doublereal scale, bqra21, brqa21;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    static doublereal licop[16]	/* was [4][4] */;
    static integer linfo;
    static doublereal ircop[16]	/* was [4][4] */;
    extern /* Subroutine */ int mb03yt_(doublereal *, integer *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *);
    static doublereal dnorm;
    extern /* Subroutine */ int sb04ow_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *);
    static doublereal dwork[32];
    static integer iwork[4];
    extern /* Subroutine */ int dgeqr2_(integer *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *), dgerq2_(
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *), dorg2r_(integer *, integer *, integer *,
	     doublereal *, integer *, doublereal *, doublereal *, integer *), 
	    dorgr2_(integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *), dorm2r_(char *, char *, 
	    integer *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    ftnlen, ftnlen), dormr2_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, ftnlen, ftnlen);
    extern doublereal dlamch_(char *, ftnlen);
    static doublereal dscale;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dlaset_(char *, integer *, integer *, doublereal *,
	     doublereal *, doublereal *, integer *, ftnlen), dlassq_(integer *
	    , doublereal *, integer *, doublereal *, doublereal *);
    static logical dtrong;
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

/*     To swap adjacent diagonal blocks A11*B11 and A22*B22 of size */
/*     1-by-1 or 2-by-2 in an upper (quasi) triangular matrix product */
/*     A*B by an orthogonal equivalence transformation. */

/*     (A, B) must be in periodic real Schur canonical form (as returned */
/*     by SLICOT Library routine MB03XP), i.e., A is block upper */
/*     triangular with 1-by-1 and 2-by-2 diagonal blocks, and B is upper */
/*     triangular. */

/*     Optionally, the matrices Q and Z of generalized Schur vectors are */
/*     updated. */

/*         Q(in) * A(in) * Z(in)' = Q(out) * A(out) * Z(out)', */
/*         Z(in) * B(in) * Q(in)' = Z(out) * B(out) * Q(out)'. */

/*     This routine is largely based on the LAPACK routine DTGEX2 */
/*     developed by Bo Kagstrom and Peter Poromaa. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     WANTQ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Q as follows: */
/*             = .TRUE. :  The matrix Q is updated; */
/*             = .FALSE.:  the matrix Q is not required. */

/*     WANTZ   LOGICAL */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrix Z as follows: */
/*             = .TRUE. :  The matrix Z is updated; */
/*             = .FALSE.:  the matrix Z is not required. */

/*     Input/Output Parameters */

/*     N1      (input) INTEGER */
/*             The order of the first block A11*B11. N1 = 0, 1 or 2. */

/*     N2      (input) INTEGER */
/*             The order of the second block A22*B22. N2 = 0, 1 or 2. */

/*     A       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDA,N1+N2) */
/*             On entry, the leading (N1+N2)-by-(N1+N2) part of this */
/*             array must contain the matrix A. */
/*             On exit, the leading (N1+N2)-by-(N1+N2) part of this array */
/*             contains the matrix A of the reordered pair. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A. LDA >= MAX(1,N1+N2). */

/*     B       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDB,N1+N2) */
/*             On entry, the leading (N1+N2)-by-(N1+N2) part of this */
/*             array must contain the matrix B. */
/*             On exit, the leading (N1+N2)-by-(N1+N2) part of this array */
/*             contains the matrix B of the reordered pair. */

/*     LDB     INTEGER */
/*             The leading dimension of the array B. LDB >= MAX(1,N1+N2). */

/*     Q       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDQ,N1+N2) */
/*             On entry, if WANTQ = .TRUE., the leading */
/*             (N1+N2)-by-(N1+N2) part of this array must contain the */
/*             orthogonal matrix Q. */
/*             On exit, the leading (N1+N2)-by-(N1+N2) part of this array */
/*             contains the updated matrix Q. Q will be a rotation */
/*             matrix for N1=N2=1. */
/*             This array is not referenced if WANTQ = .FALSE.. */

/*     LDQ     INTEGER */
/*             The leading dimension of the array Q. LDQ >= 1. */
/*             If WANTQ = .TRUE., LDQ >= N1+N2. */

/*     Z       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDZ,N1+N2) */
/*             On entry, if WANTZ = .TRUE., the leading */
/*             (N1+N2)-by-(N1+N2) part of this array must contain the */
/*             orthogonal matrix Z. */
/*             On exit, the leading (N1+N2)-by-(N1+N2) part of this array */
/*             contains the updated matrix Z. Z will be a rotation */
/*             matrix for N1=N2=1. */
/*             This array is not referenced if WANTZ = .FALSE.. */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z. LDZ >= 1. */
/*             If WANTZ = .TRUE., LDZ >= N1+N2. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             = 1:  the transformed matrix (A, B) would be */
/*                   too far from periodic Schur form; the blocks are */
/*                   not swapped and (A,B) and (Q,Z) are unchanged. */

/*     METHOD */

/*     In the current code both weak and strong stability tests are */
/*     performed. The user can omit the strong stability test by changing */
/*     the internal logical parameter WANDS to .FALSE.. See ref. [2] for */
/*     details. */

/*     REFERENCES */

/*     [1] Kagstrom, B. */
/*         A direct method for reordering eigenvalues in the generalized */
/*         real Schur form of a regular matrix pair (A,B), in M.S. Moonen */
/*         et al (eds.), Linear Algebra for Large Scale and Real-Time */
/*         Applications, Kluwer Academic Publ., 1993, pp. 195-218. */

/*     [2] Kagstrom, B., and Poromaa, P. */
/*         Computing eigenspaces with specified eigenvalues of a regular */
/*         matrix pair (A, B) and condition estimation: Theory, */
/*         algorithms and software, Numer. Algorithms, 1996, vol. 12, */
/*         pp. 369-407. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DTGPX2). */

/*     KEYWORDS */

/*     Eigenvalue, periodic Schur form, reordering */

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

#line 196 "MB03WA.f"
    /* Parameter adjustments */
#line 196 "MB03WA.f"
    a_dim1 = *lda;
#line 196 "MB03WA.f"
    a_offset = 1 + a_dim1;
#line 196 "MB03WA.f"
    a -= a_offset;
#line 196 "MB03WA.f"
    b_dim1 = *ldb;
#line 196 "MB03WA.f"
    b_offset = 1 + b_dim1;
#line 196 "MB03WA.f"
    b -= b_offset;
#line 196 "MB03WA.f"
    q_dim1 = *ldq;
#line 196 "MB03WA.f"
    q_offset = 1 + q_dim1;
#line 196 "MB03WA.f"
    q -= q_offset;
#line 196 "MB03WA.f"
    z_dim1 = *ldz;
#line 196 "MB03WA.f"
    z_offset = 1 + z_dim1;
#line 196 "MB03WA.f"
    z__ -= z_offset;
#line 196 "MB03WA.f"

#line 196 "MB03WA.f"
    /* Function Body */
#line 196 "MB03WA.f"
    *info = 0;

/*     Quick return if possible. */
/*     For efficiency, the arguments are not checked. */

#line 201 "MB03WA.f"
    if (*n1 <= 0 || *n2 <= 0) {
#line 201 "MB03WA.f"
	return 0;
#line 201 "MB03WA.f"
    }
#line 203 "MB03WA.f"
    m = *n1 + *n2;

#line 205 "MB03WA.f"
    weak = FALSE_;
#line 206 "MB03WA.f"
    dtrong = FALSE_;

/*     Make a local copy of selected block. */

#line 210 "MB03WA.f"
    dlaset_("All", &c__4, &c__4, &c_b5, &c_b5, li, &c__4, (ftnlen)3);
#line 211 "MB03WA.f"
    dlaset_("All", &c__4, &c__4, &c_b5, &c_b5, ir, &c__4, (ftnlen)3);
#line 212 "MB03WA.f"
    dlacpy_("Full", &m, &m, &a[a_offset], lda, s, &c__4, (ftnlen)4);
#line 213 "MB03WA.f"
    dlacpy_("Full", &m, &m, &b[b_offset], ldb, t, &c__4, (ftnlen)4);

/*     Compute threshold for testing acceptance of swapping. */

#line 217 "MB03WA.f"
    eps = dlamch_("P", (ftnlen)1);
#line 218 "MB03WA.f"
    smlnum = dlamch_("S", (ftnlen)1) / eps;
#line 219 "MB03WA.f"
    dscale = 0.;
#line 220 "MB03WA.f"
    dsum = 1.;
#line 221 "MB03WA.f"
    dlacpy_("Full", &m, &m, s, &c__4, dwork, &m, (ftnlen)4);
#line 222 "MB03WA.f"
    i__1 = m * m;
#line 222 "MB03WA.f"
    dlassq_(&i__1, dwork, &c__1, &dscale, &dsum);
#line 223 "MB03WA.f"
    dlacpy_("Full", &m, &m, t, &c__4, dwork, &m, (ftnlen)4);
#line 224 "MB03WA.f"
    i__1 = m * m;
#line 224 "MB03WA.f"
    dlassq_(&i__1, dwork, &c__1, &dscale, &dsum);
#line 225 "MB03WA.f"
    dnorm = dscale * sqrt(dsum);
/* Computing MAX */
#line 226 "MB03WA.f"
    d__1 = eps * 10. * dnorm;
#line 226 "MB03WA.f"
    thresh = max(d__1,smlnum);

#line 228 "MB03WA.f"
    if (m == 2) {

/*        CASE 1: Swap 1-by-1 and 1-by-1 blocks. */

/*        Compute orthogonal QL and RQ that swap 1-by-1 and 1-by-1 blocks */
/*        using Givens rotations and perform the swap tentatively. */

#line 235 "MB03WA.f"
	f = s[5] * t[5] - t[0] * s[0];
#line 236 "MB03WA.f"
	g = -s[5] * t[4] - t[0] * s[4];
#line 237 "MB03WA.f"
	sb = abs(t[0]);
#line 238 "MB03WA.f"
	sa = abs(s[5]);
#line 239 "MB03WA.f"
	dlartg_(&f, &g, &ir[4], ir, &ddum);
#line 240 "MB03WA.f"
	ir[1] = -ir[4];
#line 241 "MB03WA.f"
	ir[5] = ir[0];
#line 242 "MB03WA.f"
	drot_(&c__2, s, &c__1, &s[4], &c__1, ir, &ir[1]);
#line 243 "MB03WA.f"
	drot_(&c__2, t, &c__4, &t[1], &c__4, ir, &ir[1]);
#line 244 "MB03WA.f"
	if (sa >= sb) {
#line 245 "MB03WA.f"
	    dlartg_(s, &s[1], li, &li[1], &ddum);
#line 246 "MB03WA.f"
	} else {
#line 247 "MB03WA.f"
	    dlartg_(&t[5], &t[1], li, &li[1], &ddum);
#line 248 "MB03WA.f"
	    li[1] = -li[1];
#line 249 "MB03WA.f"
	}
#line 250 "MB03WA.f"
	drot_(&c__2, s, &c__4, &s[1], &c__4, li, &li[1]);
#line 251 "MB03WA.f"
	drot_(&c__2, t, &c__1, &t[4], &c__1, li, &li[1]);
#line 252 "MB03WA.f"
	li[5] = li[0];
#line 253 "MB03WA.f"
	li[4] = -li[1];

/*        Weak stability test: */
/*           |S21| + |T21| <= O(EPS * F-norm((S, T))). */

#line 258 "MB03WA.f"
	ws = abs(s[1]) + abs(t[1]);
#line 259 "MB03WA.f"
	weak = ws <= thresh;
#line 260 "MB03WA.f"
	if (! weak) {
#line 260 "MB03WA.f"
	    goto L50;
#line 260 "MB03WA.f"
	}

#line 263 "MB03WA.f"
	if (TRUE_) {

/*           Strong stability test: */
/*             F-norm((A-QL'*S*QR, B-QR'*T*QL)) <= O(EPS*F-norm((A,B))). */

#line 268 "MB03WA.f"
	    dlacpy_("Full", &m, &m, &a[a_offset], lda, &dwork[m * m], &m, (
		    ftnlen)4);
#line 269 "MB03WA.f"
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, li, &
		    c__4, s, &c__4, &c_b5, dwork, &m, (ftnlen)12, (ftnlen)12);
#line 271 "MB03WA.f"
	    dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b48, dwork, &m,
		     ir, &c__4, &c_b42, &dwork[m * m], &m, (ftnlen)12, (
		    ftnlen)9);
#line 273 "MB03WA.f"
	    dscale = 0.;
#line 274 "MB03WA.f"
	    dsum = 1.;
#line 275 "MB03WA.f"
	    i__1 = m * m;
#line 275 "MB03WA.f"
	    dlassq_(&i__1, &dwork[m * m], &c__1, &dscale, &dsum);

#line 277 "MB03WA.f"
	    dlacpy_("Full", &m, &m, &b[b_offset], ldb, &dwork[m * m], &m, (
		    ftnlen)4);
#line 278 "MB03WA.f"
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, ir, &
		    c__4, t, &c__4, &c_b5, dwork, &m, (ftnlen)12, (ftnlen)12);
#line 280 "MB03WA.f"
	    dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b48, dwork, &m,
		     li, &c__4, &c_b42, &dwork[m * m], &m, (ftnlen)12, (
		    ftnlen)9);
#line 282 "MB03WA.f"
	    i__1 = m * m;
#line 282 "MB03WA.f"
	    dlassq_(&i__1, &dwork[m * m], &c__1, &dscale, &dsum);
#line 283 "MB03WA.f"
	    ss = dscale * sqrt(dsum);
#line 284 "MB03WA.f"
	    dtrong = ss <= thresh;
#line 285 "MB03WA.f"
	    if (! dtrong) {
#line 285 "MB03WA.f"
		goto L50;
#line 285 "MB03WA.f"
	    }
#line 287 "MB03WA.f"
	}

/*        Update A and B. */

#line 291 "MB03WA.f"
	dlacpy_("All", &m, &m, s, &c__4, &a[a_offset], lda, (ftnlen)3);
#line 292 "MB03WA.f"
	dlacpy_("All", &m, &m, t, &c__4, &b[b_offset], ldb, (ftnlen)3);

/*        Set  N1-by-N2 (2,1) - blocks to ZERO. */

#line 296 "MB03WA.f"
	a[a_dim1 + 2] = 0.;
#line 297 "MB03WA.f"
	b[b_dim1 + 2] = 0.;

/*        Accumulate transformations into Q and Z if requested. */

#line 301 "MB03WA.f"
	if (*wantq) {
#line 301 "MB03WA.f"
	    drot_(&c__2, &q[q_dim1 + 1], &c__1, &q[(q_dim1 << 1) + 1], &c__1, 
		    li, &li[1]);
#line 301 "MB03WA.f"
	}
#line 303 "MB03WA.f"
	if (*wantz) {
#line 303 "MB03WA.f"
	    drot_(&c__2, &z__[z_dim1 + 1], &c__1, &z__[(z_dim1 << 1) + 1], &
		    c__1, ir, &ir[1]);
#line 303 "MB03WA.f"
	}

/*        Exit with INFO = 0 if swap was successfully performed. */

#line 308 "MB03WA.f"
	return 0;

#line 310 "MB03WA.f"
    } else {

/*        CASE 2: Swap 1-by-1 and 2-by-2 blocks, or 2-by-2 */
/*                and 2-by-2 blocks. */

/*        Solve the periodic Sylvester equation */
/*                 S11 * R - L * S22 = SCALE * S12 */
/*                 T11 * L - R * T22 = SCALE * T12 */
/*        for R and L. Solutions in IR and LI. */

#line 320 "MB03WA.f"
	dlacpy_("Full", n1, n2, &t[(*n1 + 1 << 2) - 4], &c__4, li, &c__4, (
		ftnlen)4);
#line 321 "MB03WA.f"
	dlacpy_("Full", n1, n2, &s[(*n1 + 1 << 2) - 4], &c__4, &ir[*n2 + 1 + (
		*n1 + 1 << 2) - 5], &c__4, (ftnlen)4);
#line 323 "MB03WA.f"
	sb04ow_(n1, n2, s, &c__4, &s[*n1 + 1 + (*n1 + 1 << 2) - 5], &c__4, &
		ir[*n2 + 1 + (*n1 + 1 << 2) - 5], &c__4, t, &c__4, &t[*n1 + 1 
		+ (*n1 + 1 << 2) - 5], &c__4, li, &c__4, &scale, iwork, &
		linfo);
#line 326 "MB03WA.f"
	if (linfo != 0) {
#line 326 "MB03WA.f"
	    goto L50;
#line 326 "MB03WA.f"
	}

/*        Compute orthogonal matrix QL: */

/*                    QL' * LI = [ TL ] */
/*                               [ 0  ] */
/*        where */
/*                    LI =  [      -L              ]. */
/*                          [ SCALE * identity(N2) ] */

#line 337 "MB03WA.f"
	i__1 = *n2;
#line 337 "MB03WA.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 338 "MB03WA.f"
	    dscal_(n1, &c_b48, &li[(i__ << 2) - 4], &c__1);
#line 339 "MB03WA.f"
	    li[*n1 + i__ + (i__ << 2) - 5] = scale;
#line 340 "MB03WA.f"
/* L10: */
#line 340 "MB03WA.f"
	}
#line 341 "MB03WA.f"
	dgeqr2_(&m, n2, li, &c__4, taul, dwork, &linfo);
#line 342 "MB03WA.f"
	dorg2r_(&m, &m, n2, li, &c__4, taul, dwork, &linfo);

/*        Compute orthogonal matrix RQ: */

/*                    IR * RQ' =   [ 0  TR], */

/*         where IR = [ SCALE * identity(N1), R ]. */

#line 350 "MB03WA.f"
	i__1 = *n1;
#line 350 "MB03WA.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 351 "MB03WA.f"
	    ir[*n2 + i__ + (i__ << 2) - 5] = scale;
#line 352 "MB03WA.f"
/* L20: */
#line 352 "MB03WA.f"
	}
#line 353 "MB03WA.f"
	dgerq2_(n1, &m, &ir[*n2], &c__4, taur, dwork, &linfo);
#line 354 "MB03WA.f"
	dorgr2_(&m, &m, n1, ir, &c__4, taur, dwork, &linfo);

/*        Perform the swapping tentatively: */

#line 358 "MB03WA.f"
	dgemm_("Transpose", "No Transpose", &m, &m, &m, &c_b42, li, &c__4, s, 
		&c__4, &c_b5, dwork, &m, (ftnlen)9, (ftnlen)12);
#line 360 "MB03WA.f"
	dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b42, dwork, &m, ir,
		 &c__4, &c_b5, s, &c__4, (ftnlen)12, (ftnlen)9);
#line 362 "MB03WA.f"
	dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, ir, &c__4, 
		t, &c__4, &c_b5, dwork, &m, (ftnlen)12, (ftnlen)12);
#line 364 "MB03WA.f"
	dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, dwork, &m, 
		li, &c__4, &c_b5, t, &c__4, (ftnlen)12, (ftnlen)12);
#line 366 "MB03WA.f"
	dlacpy_("All", &m, &m, s, &c__4, scpy, &c__4, (ftnlen)3);
#line 367 "MB03WA.f"
	dlacpy_("All", &m, &m, t, &c__4, tcpy, &c__4, (ftnlen)3);
#line 368 "MB03WA.f"
	dlacpy_("All", &m, &m, ir, &c__4, ircop, &c__4, (ftnlen)3);
#line 369 "MB03WA.f"
	dlacpy_("All", &m, &m, li, &c__4, licop, &c__4, (ftnlen)3);

/*        Triangularize the B-part by a QR factorization. */
/*        Apply transformation (from left) to A-part, giving S. */

#line 374 "MB03WA.f"
	dgeqr2_(&m, &m, t, &c__4, taur, dwork, &linfo);
#line 375 "MB03WA.f"
	dorm2r_("Right", "No Transpose", &m, &m, &m, t, &c__4, taur, s, &c__4,
		 dwork, &linfo, (ftnlen)5, (ftnlen)12);
#line 377 "MB03WA.f"
	dorm2r_("Left", "Transpose", &m, &m, &m, t, &c__4, taur, ir, &c__4, 
		dwork, &linfo, (ftnlen)4, (ftnlen)9);

/*        Compute F-norm(S21) in BRQA21. (T21 is 0.) */

#line 382 "MB03WA.f"
	dscale = 0.;
#line 383 "MB03WA.f"
	dsum = 1.;
#line 384 "MB03WA.f"
	i__1 = *n2;
#line 384 "MB03WA.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 385 "MB03WA.f"
	    dlassq_(n1, &s[*n2 + 1 + (i__ << 2) - 5], &c__1, &dscale, &dsum);
#line 386 "MB03WA.f"
/* L30: */
#line 386 "MB03WA.f"
	}
#line 387 "MB03WA.f"
	brqa21 = dscale * sqrt(dsum);

/*        Triangularize the B-part by an RQ factorization. */
/*        Apply transformation (from right) to A-part, giving S. */

#line 392 "MB03WA.f"
	dgerq2_(&m, &m, tcpy, &c__4, taul, dwork, &linfo);
#line 393 "MB03WA.f"
	dormr2_("Left", "No Transpose", &m, &m, &m, tcpy, &c__4, taul, scpy, &
		c__4, dwork, &linfo, (ftnlen)4, (ftnlen)12);
#line 395 "MB03WA.f"
	dormr2_("Right", "Transpose", &m, &m, &m, tcpy, &c__4, taul, licop, &
		c__4, dwork, &linfo, (ftnlen)5, (ftnlen)9);

/*        Compute F-norm(S21) in BQRA21. (T21 is 0.) */

#line 400 "MB03WA.f"
	dscale = 0.;
#line 401 "MB03WA.f"
	dsum = 1.;
#line 402 "MB03WA.f"
	i__1 = *n2;
#line 402 "MB03WA.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 403 "MB03WA.f"
	    dlassq_(n1, &scpy[*n2 + 1 + (i__ << 2) - 5], &c__1, &dscale, &
		    dsum);
#line 404 "MB03WA.f"
/* L40: */
#line 404 "MB03WA.f"
	}
#line 405 "MB03WA.f"
	bqra21 = dscale * sqrt(dsum);

/*        Decide which method to use. */
/*          Weak stability test: */
/*             F-norm(S21) <= O(EPS * F-norm((S, T))) */

#line 411 "MB03WA.f"
	if (bqra21 <= brqa21 && bqra21 <= thresh) {
#line 412 "MB03WA.f"
	    dlacpy_("All", &m, &m, scpy, &c__4, s, &c__4, (ftnlen)3);
#line 413 "MB03WA.f"
	    dlacpy_("All", &m, &m, tcpy, &c__4, t, &c__4, (ftnlen)3);
#line 414 "MB03WA.f"
	    dlacpy_("All", &m, &m, ircop, &c__4, ir, &c__4, (ftnlen)3);
#line 415 "MB03WA.f"
	    dlacpy_("All", &m, &m, licop, &c__4, li, &c__4, (ftnlen)3);
#line 416 "MB03WA.f"
	} else if (brqa21 >= thresh) {
#line 417 "MB03WA.f"
	    goto L50;
#line 418 "MB03WA.f"
	}

/*        Set lower triangle of B-part to zero */

#line 422 "MB03WA.f"
	i__1 = m - 1;
#line 422 "MB03WA.f"
	i__2 = m - 1;
#line 422 "MB03WA.f"
	dlaset_("Lower", &i__1, &i__2, &c_b5, &c_b5, &t[1], &c__4, (ftnlen)5);

#line 424 "MB03WA.f"
	if (TRUE_) {

/*           Strong stability test: */
/*              F-norm((A-QL*S*QR', B-QR*T*QL')) <= O(EPS*F-norm((A,B))) */

#line 429 "MB03WA.f"
	    dlacpy_("All", &m, &m, &a[a_offset], lda, &dwork[m * m], &m, (
		    ftnlen)3);
#line 430 "MB03WA.f"
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, li, &
		    c__4, s, &c__4, &c_b5, dwork, &m, (ftnlen)12, (ftnlen)12);
#line 432 "MB03WA.f"
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b48, dwork, 
		    &m, ir, &c__4, &c_b42, &dwork[m * m], &m, (ftnlen)12, (
		    ftnlen)12);
#line 434 "MB03WA.f"
	    dscale = 0.;
#line 435 "MB03WA.f"
	    dsum = 1.;
#line 436 "MB03WA.f"
	    i__1 = m * m;
#line 436 "MB03WA.f"
	    dlassq_(&i__1, &dwork[m * m], &c__1, &dscale, &dsum);

#line 438 "MB03WA.f"
	    dlacpy_("All", &m, &m, &b[b_offset], ldb, &dwork[m * m], &m, (
		    ftnlen)3);
#line 439 "MB03WA.f"
	    dgemm_("Transpose", "No Transpose", &m, &m, &m, &c_b42, ir, &c__4,
		     t, &c__4, &c_b5, dwork, &m, (ftnlen)9, (ftnlen)12);
#line 441 "MB03WA.f"
	    dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b48, dwork, &m,
		     li, &c__4, &c_b42, &dwork[m * m], &m, (ftnlen)12, (
		    ftnlen)9);
#line 443 "MB03WA.f"
	    i__1 = m * m;
#line 443 "MB03WA.f"
	    dlassq_(&i__1, &dwork[m * m], &c__1, &dscale, &dsum);
#line 444 "MB03WA.f"
	    ss = dscale * sqrt(dsum);
#line 445 "MB03WA.f"
	    dtrong = ss <= thresh;
#line 446 "MB03WA.f"
	    if (! dtrong) {
#line 446 "MB03WA.f"
		goto L50;
#line 446 "MB03WA.f"
	    }

#line 449 "MB03WA.f"
	}

/*        If the swap is accepted ("weakly" and "strongly"), apply the */
/*        transformations and set N1-by-N2 (2,1)-block to zero. */

#line 454 "MB03WA.f"
	dlaset_("All", n1, n2, &c_b5, &c_b5, &s[*n2], &c__4, (ftnlen)3);

/*        Copy (S,T) to (A,B). */

#line 458 "MB03WA.f"
	dlacpy_("All", &m, &m, s, &c__4, &a[a_offset], lda, (ftnlen)3);
#line 459 "MB03WA.f"
	dlacpy_("All", &m, &m, t, &c__4, &b[b_offset], ldb, (ftnlen)3);
#line 460 "MB03WA.f"
	dlaset_("All", &c__4, &c__4, &c_b5, &c_b5, t, &c__4, (ftnlen)3);

/*        Standardize existing 2-by-2 blocks. */

#line 464 "MB03WA.f"
	dlaset_("All", &m, &m, &c_b5, &c_b5, dwork, &m, (ftnlen)3);
#line 465 "MB03WA.f"
	dwork[0] = 1.;
#line 466 "MB03WA.f"
	t[0] = 1.;
#line 467 "MB03WA.f"
	if (*n2 > 1) {
#line 468 "MB03WA.f"
	    mb03yt_(&a[a_offset], lda, &b[b_offset], ldb, ar, ai, be, dwork, &
		    dwork[1], t, &t[1]);
#line 470 "MB03WA.f"
	    dwork[m] = -dwork[1];
#line 471 "MB03WA.f"
	    dwork[m + 1] = dwork[0];
#line 472 "MB03WA.f"
	    t[*n2 + (*n2 << 2) - 5] = t[0];
#line 473 "MB03WA.f"
	    t[4] = -t[1];
#line 474 "MB03WA.f"
	}
#line 475 "MB03WA.f"
	dwork[m * m - 1] = 1.;
#line 476 "MB03WA.f"
	t[m + (m << 2) - 5] = 1.;

#line 478 "MB03WA.f"
	if (*n1 > 1) {
#line 479 "MB03WA.f"
	    mb03yt_(&a[*n2 + 1 + (*n2 + 1) * a_dim1], lda, &b[*n2 + 1 + (*n2 
		    + 1) * b_dim1], ldb, taur, taul, &dwork[m * m], &dwork[*
		    n2 * m + *n2], &dwork[*n2 * m + *n2 + 1], &t[*n2 + 1 + (*
		    n2 + 1 << 2) - 5], &t[m + (m - 1 << 2) - 5]);
#line 482 "MB03WA.f"
	    dwork[m * m - 1] = dwork[*n2 * m + *n2];
#line 483 "MB03WA.f"
	    dwork[m * m - 2] = -dwork[*n2 * m + *n2 + 1];
#line 484 "MB03WA.f"
	    t[m + (m << 2) - 5] = t[*n2 + 1 + (*n2 + 1 << 2) - 5];
#line 485 "MB03WA.f"
	    t[m - 1 + (m << 2) - 5] = -t[m + (m - 1 << 2) - 5];
#line 486 "MB03WA.f"
	}

#line 488 "MB03WA.f"
	dgemm_("Transpose", "No Transpose", n2, n1, n2, &c_b42, dwork, &m, &a[
		(*n2 + 1) * a_dim1 + 1], lda, &c_b5, &dwork[m * m], n2, (
		ftnlen)9, (ftnlen)12);
#line 490 "MB03WA.f"
	dlacpy_("All", n2, n1, &dwork[m * m], n2, &a[(*n2 + 1) * a_dim1 + 1], 
		lda, (ftnlen)3);
#line 491 "MB03WA.f"
	dgemm_("Transpose", "No Transpose", n2, n1, n2, &c_b42, t, &c__4, &b[(
		*n2 + 1) * b_dim1 + 1], ldb, &c_b5, &dwork[m * m], n2, (
		ftnlen)9, (ftnlen)12);
#line 494 "MB03WA.f"
	dlacpy_("All", n2, n1, &dwork[m * m], n2, &b[(*n2 + 1) * b_dim1 + 1], 
		ldb, (ftnlen)3);
#line 495 "MB03WA.f"
	dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, li, &c__4, 
		dwork, &m, &c_b5, &dwork[m * m], &m, (ftnlen)12, (ftnlen)12);
#line 497 "MB03WA.f"
	dlacpy_("All", &m, &m, &dwork[m * m], &m, li, &c__4, (ftnlen)3);
#line 498 "MB03WA.f"
	dgemm_("No Transpose", "No Transpose", n2, n1, n1, &c_b42, &a[(*n2 + 
		1) * a_dim1 + 1], lda, &t[*n2 + 1 + (*n2 + 1 << 2) - 5], &
		c__4, &c_b5, &dwork[m * m], &m, (ftnlen)12, (ftnlen)12);
#line 501 "MB03WA.f"
	dlacpy_("All", n2, n1, &dwork[m * m], &m, &a[(*n2 + 1) * a_dim1 + 1], 
		lda, (ftnlen)3);
#line 502 "MB03WA.f"
	dgemm_("No Transpose", "No Transpose", n2, n1, n1, &c_b42, &b[(*n2 + 
		1) * b_dim1 + 1], ldb, &dwork[*n2 * m + *n2], &m, &c_b5, &
		dwork[m * m], &m, (ftnlen)12, (ftnlen)12);
#line 505 "MB03WA.f"
	dlacpy_("All", n2, n1, &dwork[m * m], &m, &b[(*n2 + 1) * b_dim1 + 1], 
		ldb, (ftnlen)3);
#line 506 "MB03WA.f"
	dgemm_("Transpose", "No Transpose", &m, &m, &m, &c_b42, t, &c__4, ir, 
		&c__4, &c_b5, dwork, &m, (ftnlen)9, (ftnlen)12);
#line 508 "MB03WA.f"
	dlacpy_("All", &m, &m, dwork, &m, ir, &c__4, (ftnlen)3);

/*        Accumulate transformations into Q and Z if requested. */

#line 512 "MB03WA.f"
	if (*wantq) {
#line 513 "MB03WA.f"
	    dgemm_("No Transpose", "No Transpose", &m, &m, &m, &c_b42, &q[
		    q_offset], ldq, li, &c__4, &c_b5, dwork, &m, (ftnlen)12, (
		    ftnlen)12);
#line 515 "MB03WA.f"
	    dlacpy_("All", &m, &m, dwork, &m, &q[q_offset], ldq, (ftnlen)3);
#line 516 "MB03WA.f"
	}

#line 518 "MB03WA.f"
	if (*wantz) {
#line 519 "MB03WA.f"
	    dgemm_("No Transpose", "Transpose", &m, &m, &m, &c_b42, &z__[
		    z_offset], ldz, ir, &c__4, &c_b5, dwork, &m, (ftnlen)12, (
		    ftnlen)9);
#line 521 "MB03WA.f"
	    dlacpy_("Full", &m, &m, dwork, &m, &z__[z_offset], ldz, (ftnlen)4)
		    ;

#line 523 "MB03WA.f"
	}

/*        Exit with INFO = 0 if swap was successfully performed. */

#line 527 "MB03WA.f"
	return 0;

#line 529 "MB03WA.f"
    }

/*     Exit with INFO = 1 if swap was rejected. */

#line 533 "MB03WA.f"
L50:

#line 535 "MB03WA.f"
    *info = 1;
#line 536 "MB03WA.f"
    return 0;
/* *** Last line of MB03WA *** */
} /* mb03wa_ */

