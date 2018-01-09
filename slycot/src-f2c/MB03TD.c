#line 1 "MB03TD.f"
/* MB03TD.f -- translated by f2c (version 20100827).
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

#line 1 "MB03TD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int mb03td_(char *typ, char *compu, logical *select, logical 
	*lower, integer *n, doublereal *a, integer *lda, doublereal *g, 
	integer *ldg, doublereal *u1, integer *ldu1, doublereal *u2, integer *
	ldu2, doublereal *wr, doublereal *wi, integer *m, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen typ_len, ftnlen compu_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, g_dim1, g_offset, u1_dim1, u1_offset, u2_dim1, 
	    u2_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer k, ks, nbf, nbl, here;
    static logical pair;
    static integer ierr, ifst;
    static logical flow, swap;
    static integer ilst;
    static logical isham;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb03ts_(logical *, logical *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, integer *, integer *, integer 
	    *, doublereal *, integer *);
    static logical wantu;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer nbnext, wrkmin;


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

/*     To reorder a matrix X in skew-Hamiltonian Schur form: */

/*                   [  A   G  ]          T */
/*             X  =  [       T ],   G = -G, */
/*                   [  0   A  ] */

/*     or in Hamiltonian Schur form: */

/*                   [  A   G  ]          T */
/*             X  =  [       T ],   G =  G, */
/*                   [  0  -A  ] */

/*     where A is in upper quasi-triangular form, so that a selected */
/*     cluster of eigenvalues appears in the leading diagonal blocks */
/*     of the matrix A (in X) and the leading columns of [ U1; -U2 ] form */
/*     an orthonormal basis for the corresponding right invariant */
/*     subspace. */

/*     If X is skew-Hamiltonian, then each eigenvalue appears twice; one */
/*     copy corresponds to the j-th diagonal element and the other to the */
/*     (n+j)-th diagonal element of X. The logical array LOWER controls */
/*     which copy is to be reordered to the leading part of A. */

/*     If X is Hamiltonian then the eigenvalues appear in pairs */
/*     (lambda,-lambda); lambda corresponds to the j-th diagonal */
/*     element and -lambda to the (n+j)-th diagonal element of X. */
/*     The logical array LOWER controls whether lambda or -lambda is to */
/*     be reordered to the leading part of A. */

/*     The matrix A must be in Schur canonical form (as returned by the */
/*     LAPACK routine DHSEQR), that is, block upper triangular with */
/*     1-by-1 and 2-by-2 diagonal blocks; each 2-by-2 diagonal block has */
/*     its diagonal elements equal and its off-diagonal elements of */
/*     opposite sign. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     TYP     CHARACTER*1 */
/*             Specifies the type of the input matrix X: */
/*             = 'S': X is skew-Hamiltonian; */
/*             = 'H': X is Hamiltonian. */

/*     COMPU   CHARACTER*1 */
/*             = 'U': update the matrices U1 and U2 containing the */
/*                    Schur vectors; */
/*             = 'N': do not update U1 and U2. */

/*     SELECT  (input/output) LOGICAL array, dimension (N) */
/*             SELECT specifies the eigenvalues in the selected cluster. */
/*             To select a real eigenvalue w(j), SELECT(j) must be set */
/*             to .TRUE.. To select a complex conjugate pair of */
/*             eigenvalues w(j) and w(j+1), corresponding to a 2-by-2 */
/*             diagonal block, both SELECT(j) and SELECT(j+1) must be set */
/*             to .TRUE.; a complex conjugate pair of eigenvalues must be */
/*             either both included in the cluster or both excluded. */

/*     LOWER   (input/output) LOGICAL array, dimension (N) */
/*             LOWER controls which copy of a selected eigenvalue is */
/*             included in the cluster. If SELECT(j) is set to .TRUE. */
/*             for a real eigenvalue w(j); then LOWER(j) must be set to */
/*             .TRUE. if the eigenvalue corresponding to the (n+j)-th */
/*             diagonal element of X is to be reordered to the leading */
/*             part; and LOWER(j) must be set to .FALSE. if the */
/*             eigenvalue corresponding to the j-th diagonal element of */
/*             X is to be reordered to the leading part. Similarly, for */
/*             a complex conjugate pair of eigenvalues w(j) and w(j+1), */
/*             both LOWER(j) and LOWER(j+1) must be set to .TRUE. if the */
/*             eigenvalues corresponding to the (n+j:n+j+1,n+j:n+j+1) */
/*             diagonal block of X are to be reordered to the leading */
/*             part. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A. N >= 0. */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the upper quasi-triangular matrix A in Schur */
/*             canonical form. */
/*             On exit, the leading N-by-N part of this array contains */
/*             the reordered matrix A, again in Schur canonical form, */
/*             with the selected eigenvalues in the diagonal blocks. */

/*     LDA     INTEGER */
/*             The leading dimension of the array A.  LDA >= MAX(1,N). */

/*     G       (input/output) DOUBLE PRECISION array, dimension (LDG,N) */
/*             On entry, if TYP = 'S', the leading N-by-N part of this */
/*             array must contain the strictly upper triangular part of */
/*             the skew-symmetric matrix G. The rest of this array is not */
/*             referenced. */
/*             On entry, if TYP = 'H', the leading N-by-N part of this */
/*             array must contain the upper triangular part of the */
/*             symmetric matrix G. The rest of this array is not */
/*             referenced. */
/*             On exit, if TYP = 'S', the leading N-by-N part of this */
/*             array contains the strictly upper triangular part of the */
/*             skew-symmetric matrix G, updated by the orthogonal */
/*             symplectic transformation which reorders X. */
/*             On exit, if TYP = 'H', the leading N-by-N part of this */
/*             array contains the upper triangular part of the symmetric */
/*             matrix G, updated by the orthogonal symplectic */
/*             transformation which reorders X. */

/*     LDG     INTEGER */
/*             The leading dimension of the array G.  LDG >= MAX(1,N). */

/*     U1      (input/output) DOUBLE PRECISION array, dimension (LDU1,N) */
/*             On entry, if COMPU = 'U', the leading N-by-N part of this */
/*             array must contain U1, the (1,1) block of an orthogonal */
/*             symplectic matrix U = [ U1, U2; -U2, U1 ]. */
/*             On exit, if COMPU = 'U', the leading N-by-N part of this */
/*             array contains the (1,1) block of the matrix U, */
/*             postmultiplied by the orthogonal symplectic transformation */
/*             which reorders X. The leading M columns of U form an */
/*             orthonormal basis for the specified invariant subspace. */
/*             If COMPU = 'N', this array is not referenced. */

/*     LDU1    INTEGER */
/*             The leading dimension of the array U1. */
/*             LDU1 >= MAX(1,N),  if COMPU = 'U'; */
/*             LDU1 >= 1,         otherwise. */

/*     U2      (input/output) DOUBLE PRECISION array, dimension (LDU2,N) */
/*             On entry, if COMPU = 'U', the leading N-by-N part of this */
/*             array must contain U2, the (1,2) block of an orthogonal */
/*             symplectic matrix U = [ U1, U2; -U2, U1 ]. */
/*             On exit, if COMPU = 'U', the leading N-by-N part of this */
/*             array contains the (1,2) block of the matrix U, */
/*             postmultiplied by the orthogonal symplectic transformation */
/*             which reorders X. */
/*             If COMPU = 'N', this array is not referenced. */

/*     LDU2    INTEGER */
/*             The leading dimension of the array U2. */
/*             LDU2 >= MAX(1,N),  if COMPU = 'U'; */
/*             LDU2 >= 1,         otherwise. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             The real and imaginary parts, respectively, of the */
/*             reordered eigenvalues of A. The eigenvalues are stored */
/*             in the same order as on the diagonal of A, with */
/*             WR(i) = A(i,i) and, if A(i:i+1,i:i+1) is a 2-by-2 diagonal */
/*             block, WI(i) > 0 and WI(i+1) = -WI(i). Note that if an */
/*             eigenvalue is sufficiently ill-conditioned, then its value */
/*             may differ significantly from its value before reordering. */

/*     M       (output) INTEGER */
/*             The dimension of the specified invariant subspace. */
/*             0 <= M <= N. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0,  DWORK(1)  returns the optimal */
/*             value of LDWORK. */
/*             On exit, if  INFO = -18,  DWORK(1)  returns the minimum */
/*             value of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= MAX(1,N). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */
/*             = 1:  reordering of X failed because some eigenvalue pairs */
/*                   are too close to separate (the problem is very */
/*                   ill-conditioned); X may have been partially */
/*                   reordered, and WR and WI contain the eigenvalues in */
/*                   the same order as in X. */

/*     REFERENCES */

/*     [1] Bai, Z. and Demmel, J.W. */
/*         On Swapping Diagonal Blocks in Real Schur Form. */
/*         Linear Algebra Appl., 186, pp. 73-95, 1993. */

/*     [2] Benner, P., Kressner, D., and Mehrmann, V. */
/*         Skew-Hamiltonian and Hamiltonian Eigenvalue Problems: Theory, */
/*         Algorithms and Applications. Techn. Report, TU Berlin, 2003. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, May 2008 (SLICOT version of the HAPACK routine DHAORD). */

/*     KEYWORDS */

/*     Hamiltonian matrix, skew-Hamiltonian matrix, invariant subspace. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

/*     Decode and check input parameters. */

#line 254 "MB03TD.f"
    /* Parameter adjustments */
#line 254 "MB03TD.f"
    --select;
#line 254 "MB03TD.f"
    --lower;
#line 254 "MB03TD.f"
    a_dim1 = *lda;
#line 254 "MB03TD.f"
    a_offset = 1 + a_dim1;
#line 254 "MB03TD.f"
    a -= a_offset;
#line 254 "MB03TD.f"
    g_dim1 = *ldg;
#line 254 "MB03TD.f"
    g_offset = 1 + g_dim1;
#line 254 "MB03TD.f"
    g -= g_offset;
#line 254 "MB03TD.f"
    u1_dim1 = *ldu1;
#line 254 "MB03TD.f"
    u1_offset = 1 + u1_dim1;
#line 254 "MB03TD.f"
    u1 -= u1_offset;
#line 254 "MB03TD.f"
    u2_dim1 = *ldu2;
#line 254 "MB03TD.f"
    u2_offset = 1 + u2_dim1;
#line 254 "MB03TD.f"
    u2 -= u2_offset;
#line 254 "MB03TD.f"
    --wr;
#line 254 "MB03TD.f"
    --wi;
#line 254 "MB03TD.f"
    --dwork;
#line 254 "MB03TD.f"

#line 254 "MB03TD.f"
    /* Function Body */
#line 254 "MB03TD.f"
    isham = lsame_(typ, "H", (ftnlen)1, (ftnlen)1);
#line 255 "MB03TD.f"
    wantu = lsame_(compu, "U", (ftnlen)1, (ftnlen)1);
#line 256 "MB03TD.f"
    wrkmin = max(1,*n);
#line 257 "MB03TD.f"
    *info = 0;
#line 258 "MB03TD.f"
    if (! isham && ! lsame_(typ, "S", (ftnlen)1, (ftnlen)1)) {
#line 259 "MB03TD.f"
	*info = -1;
#line 260 "MB03TD.f"
    } else if (! wantu && ! lsame_(compu, "N", (ftnlen)1, (ftnlen)1)) {
#line 261 "MB03TD.f"
	*info = -2;
#line 262 "MB03TD.f"
    } else if (*n < 0) {
#line 263 "MB03TD.f"
	*info = -5;
#line 264 "MB03TD.f"
    } else if (*lda < max(1,*n)) {
#line 265 "MB03TD.f"
	*info = -7;
#line 266 "MB03TD.f"
    } else if (*ldg < max(1,*n)) {
#line 267 "MB03TD.f"
	*info = -9;
#line 268 "MB03TD.f"
    } else if (*ldu1 < 1 || wantu && *ldu1 < *n) {
#line 269 "MB03TD.f"
	*info = -11;
#line 270 "MB03TD.f"
    } else if (*ldu2 < 1 || wantu && *ldu2 < *n) {
#line 271 "MB03TD.f"
	*info = -13;
#line 272 "MB03TD.f"
    } else if (*ldwork < wrkmin) {
#line 273 "MB03TD.f"
	*info = -18;
#line 274 "MB03TD.f"
	dwork[1] = (doublereal) wrkmin;
#line 275 "MB03TD.f"
    }

/*     Return if there were illegal values. */

#line 279 "MB03TD.f"
    if (*info != 0) {
#line 280 "MB03TD.f"
	i__1 = -(*info);
#line 280 "MB03TD.f"
	xerbla_("MB03TD", &i__1, (ftnlen)6);
#line 281 "MB03TD.f"
	return 0;
#line 282 "MB03TD.f"
    }

/*     Set M to the dimension of the specified invariant subspace. */

#line 286 "MB03TD.f"
    *m = 0;
#line 287 "MB03TD.f"
    pair = FALSE_;
#line 288 "MB03TD.f"
    i__1 = *n;
#line 288 "MB03TD.f"
    for (k = 1; k <= i__1; ++k) {
#line 289 "MB03TD.f"
	if (pair) {
#line 290 "MB03TD.f"
	    pair = FALSE_;
#line 291 "MB03TD.f"
	} else {
#line 292 "MB03TD.f"
	    if (k < *n) {
#line 293 "MB03TD.f"
		if (a[k + 1 + k * a_dim1] == 0.) {
#line 294 "MB03TD.f"
		    if (select[k]) {
#line 294 "MB03TD.f"
			++(*m);
#line 294 "MB03TD.f"
		    }
#line 296 "MB03TD.f"
		} else {
#line 297 "MB03TD.f"
		    pair = TRUE_;
#line 298 "MB03TD.f"
		    if (select[k] || select[k + 1]) {
#line 298 "MB03TD.f"
			*m += 2;
#line 298 "MB03TD.f"
		    }
#line 300 "MB03TD.f"
		}
#line 301 "MB03TD.f"
	    } else {
#line 302 "MB03TD.f"
		if (select[*n]) {
#line 302 "MB03TD.f"
		    ++(*m);
#line 302 "MB03TD.f"
		}
#line 304 "MB03TD.f"
	    }
#line 305 "MB03TD.f"
	}
#line 306 "MB03TD.f"
/* L10: */
#line 306 "MB03TD.f"
    }

/*     Quick return if possible. */

#line 310 "MB03TD.f"
    if (*n == 0) {
#line 311 "MB03TD.f"
	dwork[1] = 1.;
#line 312 "MB03TD.f"
	return 0;
#line 313 "MB03TD.f"
    }

/*     Collect the selected blocks at the top-left corner of X. */

#line 317 "MB03TD.f"
    ks = 0;
#line 318 "MB03TD.f"
    pair = FALSE_;
#line 319 "MB03TD.f"
    i__1 = *n;
#line 319 "MB03TD.f"
    for (k = 1; k <= i__1; ++k) {
#line 320 "MB03TD.f"
	if (pair) {
#line 321 "MB03TD.f"
	    pair = FALSE_;
#line 322 "MB03TD.f"
	} else {
#line 323 "MB03TD.f"
	    swap = select[k];
#line 324 "MB03TD.f"
	    flow = lower[k];
#line 325 "MB03TD.f"
	    if (k < *n) {
#line 326 "MB03TD.f"
		if (a[k + 1 + k * a_dim1] != 0.) {
#line 327 "MB03TD.f"
		    pair = TRUE_;
#line 328 "MB03TD.f"
		    swap = swap || select[k + 1];
#line 329 "MB03TD.f"
		    flow = flow || lower[k + 1];
#line 330 "MB03TD.f"
		}
#line 331 "MB03TD.f"
	    }

#line 333 "MB03TD.f"
	    if (pair) {
#line 334 "MB03TD.f"
		nbf = 2;
#line 335 "MB03TD.f"
	    } else {
#line 336 "MB03TD.f"
		nbf = 1;
#line 337 "MB03TD.f"
	    }

#line 339 "MB03TD.f"
	    if (swap) {
#line 340 "MB03TD.f"
		++ks;
#line 341 "MB03TD.f"
		if (flow) {

/*                 Step 1: Swap the K-th block to position N. */

#line 345 "MB03TD.f"
		    ifst = k;
#line 346 "MB03TD.f"
		    ilst = *n;
#line 347 "MB03TD.f"
		    nbl = 1;
#line 348 "MB03TD.f"
		    if (ilst > 1) {
#line 349 "MB03TD.f"
			if (a[ilst + (ilst - 1) * a_dim1] != 0.) {
#line 350 "MB03TD.f"
			    --ilst;
#line 351 "MB03TD.f"
			    nbl = 2;
#line 352 "MB03TD.f"
			}
#line 353 "MB03TD.f"
		    }

/*                 Update ILST. */

#line 357 "MB03TD.f"
		    if (nbf == 2 && nbl == 1) {
#line 357 "MB03TD.f"
			--ilst;
#line 357 "MB03TD.f"
		    }
#line 359 "MB03TD.f"
		    if (nbf == 1 && nbl == 2) {
#line 359 "MB03TD.f"
			++ilst;
#line 359 "MB03TD.f"
		    }

#line 362 "MB03TD.f"
		    if (ilst == ifst) {
#line 362 "MB03TD.f"
			goto L30;
#line 362 "MB03TD.f"
		    }

#line 365 "MB03TD.f"
		    here = ifst;

#line 367 "MB03TD.f"
L20:

/*                 Swap block with next one below. */

#line 371 "MB03TD.f"
		    if (nbf == 1 || nbf == 2) {

/*                    Current block is either 1-by-1 or 2-by-2. */

#line 375 "MB03TD.f"
			nbnext = 1;
#line 376 "MB03TD.f"
			if (here + nbf + 1 <= *n) {
#line 377 "MB03TD.f"
			    if (a[here + nbf + 1 + (here + nbf) * a_dim1] != 
				    0.) {
#line 377 "MB03TD.f"
				nbnext = 2;
#line 377 "MB03TD.f"
			    }
#line 379 "MB03TD.f"
			}
#line 380 "MB03TD.f"
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &here, &nbf, &nbnext, &
				dwork[1], &ierr);
#line 383 "MB03TD.f"
			if (ierr != 0) {
#line 384 "MB03TD.f"
			    *info = 1;
#line 385 "MB03TD.f"
			    goto L70;
#line 386 "MB03TD.f"
			}
#line 387 "MB03TD.f"
			here += nbnext;

/*                    Test if 2-by-2 block breaks into two 1-by-1 blocks. */

#line 391 "MB03TD.f"
			if (nbf == 2) {
#line 392 "MB03TD.f"
			    if (a[here + 1 + here * a_dim1] == 0.) {
#line 392 "MB03TD.f"
				nbf = 3;
#line 392 "MB03TD.f"
			    }
#line 394 "MB03TD.f"
			}

#line 396 "MB03TD.f"
		    } else {

/*                    Current block consists of two 1-by-1 blocks each of */
/*                    which must be swapped individually. */

#line 401 "MB03TD.f"
			nbnext = 1;
#line 402 "MB03TD.f"
			if (here + 3 <= *n) {
#line 403 "MB03TD.f"
			    if (a[here + 3 + (here + 2) * a_dim1] != 0.) {
#line 403 "MB03TD.f"
				nbnext = 2;
#line 403 "MB03TD.f"
			    }
#line 405 "MB03TD.f"
			}
#line 406 "MB03TD.f"
			i__2 = here + 1;
#line 406 "MB03TD.f"
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &i__2, &c__1, &nbnext, &
				dwork[1], &ierr);
#line 409 "MB03TD.f"
			if (ierr != 0) {
#line 410 "MB03TD.f"
			    *info = 1;
#line 411 "MB03TD.f"
			    goto L70;
#line 412 "MB03TD.f"
			}
#line 413 "MB03TD.f"
			if (nbnext == 1) {

/*                       Swap two 1-by-1 blocks, no problems possible. */

#line 417 "MB03TD.f"
			    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				    g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				    u2_offset], ldu2, &here, &c__1, &nbnext, &
				    dwork[1], &ierr);
#line 420 "MB03TD.f"
			    ++here;
#line 421 "MB03TD.f"
			} else {

/*                       Recompute NBNEXT in case 2 by 2 split. */

#line 425 "MB03TD.f"
			    if (a[here + 2 + (here + 1) * a_dim1] == 0.) {
#line 425 "MB03TD.f"
				nbnext = 1;
#line 425 "MB03TD.f"
			    }
#line 427 "MB03TD.f"
			    if (nbnext == 2) {

/*                          2-by-2 block did not split */

#line 431 "MB03TD.f"
				mb03ts_(&isham, &wantu, n, &a[a_offset], lda, 
					&g[g_offset], ldg, &u1[u1_offset], 
					ldu1, &u2[u2_offset], ldu2, &here, &
					c__1, &nbnext, &dwork[1], &ierr);
#line 434 "MB03TD.f"
				if (ierr != 0) {
#line 435 "MB03TD.f"
				    *info = 1;
#line 436 "MB03TD.f"
				    goto L70;
#line 437 "MB03TD.f"
				}
#line 438 "MB03TD.f"
				here += 2;
#line 439 "MB03TD.f"
			    } else {

/*                          2-by-2 block did split */

#line 443 "MB03TD.f"
				mb03ts_(&isham, &wantu, n, &a[a_offset], lda, 
					&g[g_offset], ldg, &u1[u1_offset], 
					ldu1, &u2[u2_offset], ldu2, &here, &
					c__1, &c__1, &dwork[1], &ierr);
#line 446 "MB03TD.f"
				i__2 = here + 1;
#line 446 "MB03TD.f"
				mb03ts_(&isham, &wantu, n, &a[a_offset], lda, 
					&g[g_offset], ldg, &u1[u1_offset], 
					ldu1, &u2[u2_offset], ldu2, &i__2, &
					c__1, &c__1, &dwork[1], &ierr);
#line 449 "MB03TD.f"
				here += 2;
#line 450 "MB03TD.f"
			    }
#line 451 "MB03TD.f"
			}
#line 452 "MB03TD.f"
		    }
#line 453 "MB03TD.f"
		    if (here < ilst) {
#line 453 "MB03TD.f"
			goto L20;
#line 453 "MB03TD.f"
		    }

#line 456 "MB03TD.f"
L30:

/*                 Step 2: Apply an orthogonal symplectic transformation */
/*                         to swap the last blocks in A and -A' (or A'). */

#line 461 "MB03TD.f"
		    if (nbf == 1) {

/*                    Exchange columns/rows N <-> 2*N. No problems */
/*                    possible. */

#line 466 "MB03TD.f"
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, n, &c__1, &c__1, &dwork[1], 
				&ierr);

#line 470 "MB03TD.f"
		    } else if (nbf == 2) {

/*                    Swap last block with its equivalent by an */
/*                    orthogonal symplectic transformation. */

#line 475 "MB03TD.f"
			i__2 = *n - 1;
#line 475 "MB03TD.f"
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &i__2, &c__2, &c__2, &dwork[
				1], &ierr);
#line 478 "MB03TD.f"
			if (ierr != 0) {
#line 479 "MB03TD.f"
			    *info = 1;
#line 480 "MB03TD.f"
			    goto L70;
#line 481 "MB03TD.f"
			}

/*                    Test if 2-by-2 block breaks into two 1-by-1 blocks. */

#line 485 "MB03TD.f"
			if (a[*n - 1 + *n * a_dim1] == 0.) {
#line 485 "MB03TD.f"
			    nbf = 3;
#line 485 "MB03TD.f"
			}
#line 487 "MB03TD.f"
		    } else {

/*                    Block did split. Swap (N-1)-th and N-th elements */
/*                    consecutively by symplectic generalized */
/*                    permutations and one rotation. */

#line 493 "MB03TD.f"
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, n, &c__1, &c__1, &dwork[1], 
				&ierr);
#line 495 "MB03TD.f"
			i__2 = *n - 1;
#line 495 "MB03TD.f"
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &i__2, &c__1, &c__1, &dwork[
				1], &ierr);
#line 498 "MB03TD.f"
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, n, &c__1, &c__1, &dwork[1], 
				&ierr);
#line 500 "MB03TD.f"
		    }
#line 501 "MB03TD.f"
		    ifst = *n;
#line 502 "MB03TD.f"
		    if (pair) {
#line 502 "MB03TD.f"
			ifst = *n - 1;
#line 502 "MB03TD.f"
		    }
#line 504 "MB03TD.f"
		} else {
#line 505 "MB03TD.f"
		    ifst = k;
#line 506 "MB03TD.f"
		}

/*              Step 3: Swap the K-th / N-th block to position KS. */

#line 510 "MB03TD.f"
		ilst = ks;
#line 511 "MB03TD.f"
		nbl = 1;
#line 512 "MB03TD.f"
		if (ilst > 1) {
#line 513 "MB03TD.f"
		    if (a[ilst + (ilst - 1) * a_dim1] != 0.) {
#line 514 "MB03TD.f"
			--ilst;
#line 515 "MB03TD.f"
			nbl = 2;
#line 516 "MB03TD.f"
		    }
#line 517 "MB03TD.f"
		}

#line 519 "MB03TD.f"
		if (ilst == ifst) {
#line 519 "MB03TD.f"
		    goto L50;
#line 519 "MB03TD.f"
		}

#line 522 "MB03TD.f"
		here = ifst;
#line 523 "MB03TD.f"
L40:

/*              Swap block with next one above. */

#line 527 "MB03TD.f"
		if (nbf == 1 || nbf == 2) {

/*                 Current block either 1 by 1 or 2 by 2. */

#line 531 "MB03TD.f"
		    nbnext = 1;
#line 532 "MB03TD.f"
		    if (here >= 3) {
#line 533 "MB03TD.f"
			if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
#line 533 "MB03TD.f"
			    nbnext = 2;
#line 533 "MB03TD.f"
			}
#line 535 "MB03TD.f"
		    }
#line 536 "MB03TD.f"
		    i__2 = here - nbnext;
#line 536 "MB03TD.f"
		    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[g_offset]
			    , ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2,
			     &i__2, &nbnext, &nbf, &dwork[1], &ierr);
#line 539 "MB03TD.f"
		    if (ierr != 0) {
#line 540 "MB03TD.f"
			*info = 1;
#line 541 "MB03TD.f"
			goto L70;
#line 542 "MB03TD.f"
		    }
#line 543 "MB03TD.f"
		    here -= nbnext;

/*                 Test if 2-by-2 block breaks into two 1-by-1 blocks. */

#line 547 "MB03TD.f"
		    if (nbf == 2) {
#line 548 "MB03TD.f"
			if (a[here + 1 + here * a_dim1] == 0.) {
#line 548 "MB03TD.f"
			    nbf = 3;
#line 548 "MB03TD.f"
			}
#line 550 "MB03TD.f"
		    }

#line 552 "MB03TD.f"
		} else {

/*                 Current block consists of two 1 by 1 blocks each of */
/*                 which must be swapped individually. */

#line 557 "MB03TD.f"
		    nbnext = 1;
#line 558 "MB03TD.f"
		    if (here >= 3) {
#line 559 "MB03TD.f"
			if (a[here - 1 + (here - 2) * a_dim1] != 0.) {
#line 559 "MB03TD.f"
			    nbnext = 2;
#line 559 "MB03TD.f"
			}
#line 561 "MB03TD.f"
		    }
#line 562 "MB03TD.f"
		    i__2 = here - nbnext;
#line 562 "MB03TD.f"
		    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[g_offset]
			    , ldg, &u1[u1_offset], ldu1, &u2[u2_offset], ldu2,
			     &i__2, &nbnext, &c__1, &dwork[1], &ierr);
#line 565 "MB03TD.f"
		    if (ierr != 0) {
#line 566 "MB03TD.f"
			*info = 1;
#line 567 "MB03TD.f"
			goto L70;
#line 568 "MB03TD.f"
		    }
#line 569 "MB03TD.f"
		    if (nbnext == 1) {

/*                    Swap two 1-by-1 blocks, no problems possible. */

#line 573 "MB03TD.f"
			mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				u2_offset], ldu2, &here, &nbnext, &c__1, &
				dwork[1], &ierr);
#line 577 "MB03TD.f"
			--here;
#line 578 "MB03TD.f"
		    } else {

/*                    Recompute NBNEXT in case 2-by-2 split. */

#line 582 "MB03TD.f"
			if (a[here + (here - 1) * a_dim1] == 0.) {
#line 582 "MB03TD.f"
			    nbnext = 1;
#line 582 "MB03TD.f"
			}
#line 584 "MB03TD.f"
			if (nbnext == 2) {

/*                       2-by-2 block did not split */

#line 588 "MB03TD.f"
			    i__2 = here - 1;
#line 588 "MB03TD.f"
			    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				    g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				    u2_offset], ldu2, &i__2, &c__2, &c__1, &
				    dwork[1], &ierr);
#line 591 "MB03TD.f"
			    if (ierr != 0) {
#line 592 "MB03TD.f"
				*info = 1;
#line 593 "MB03TD.f"
				goto L70;
#line 594 "MB03TD.f"
			    }
#line 595 "MB03TD.f"
			    here += -2;
#line 596 "MB03TD.f"
			} else {

/*                       2-by-2 block did split */

#line 600 "MB03TD.f"
			    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				    g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				    u2_offset], ldu2, &here, &c__1, &c__1, &
				    dwork[1], &ierr);
#line 603 "MB03TD.f"
			    i__2 = here - 1;
#line 603 "MB03TD.f"
			    mb03ts_(&isham, &wantu, n, &a[a_offset], lda, &g[
				    g_offset], ldg, &u1[u1_offset], ldu1, &u2[
				    u2_offset], ldu2, &i__2, &c__1, &c__1, &
				    dwork[1], &ierr);
#line 606 "MB03TD.f"
			    here += -2;
#line 607 "MB03TD.f"
			}
#line 608 "MB03TD.f"
		    }
#line 609 "MB03TD.f"
		}

#line 611 "MB03TD.f"
		if (here > ilst) {
#line 611 "MB03TD.f"
		    goto L40;
#line 611 "MB03TD.f"
		}

#line 614 "MB03TD.f"
L50:
#line 615 "MB03TD.f"
		if (pair) {
#line 615 "MB03TD.f"
		    ++ks;
#line 615 "MB03TD.f"
		}
#line 617 "MB03TD.f"
	    }
#line 618 "MB03TD.f"
	}
#line 619 "MB03TD.f"
/* L60: */
#line 619 "MB03TD.f"
    }

#line 621 "MB03TD.f"
L70:

/*     Store eigenvalues. */

#line 625 "MB03TD.f"
    i__1 = *n;
#line 625 "MB03TD.f"
    for (k = 1; k <= i__1; ++k) {
#line 626 "MB03TD.f"
	wr[k] = a[k + k * a_dim1];
#line 627 "MB03TD.f"
	wi[k] = 0.;
#line 628 "MB03TD.f"
/* L80: */
#line 628 "MB03TD.f"
    }
#line 629 "MB03TD.f"
    i__1 = *n - 1;
#line 629 "MB03TD.f"
    for (k = 1; k <= i__1; ++k) {
#line 630 "MB03TD.f"
	if (a[k + 1 + k * a_dim1] != 0.) {
#line 631 "MB03TD.f"
	    wi[k] = sqrt((d__1 = a[k + (k + 1) * a_dim1], abs(d__1))) * sqrt((
		    d__2 = a[k + 1 + k * a_dim1], abs(d__2)));
#line 633 "MB03TD.f"
	    wi[k + 1] = -wi[k];
#line 634 "MB03TD.f"
	}
#line 635 "MB03TD.f"
/* L90: */
#line 635 "MB03TD.f"
    }

#line 637 "MB03TD.f"
    dwork[1] = (doublereal) wrkmin;

#line 639 "MB03TD.f"
    return 0;
/* *** Last line of MB03TD *** */
} /* mb03td_ */

