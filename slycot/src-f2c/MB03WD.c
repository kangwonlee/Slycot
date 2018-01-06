#line 1 "MB03WD.f"
/* MB03WD.f -- translated by f2c (version 20100827).
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

#line 1 "MB03WD.f"
/* Table of constant values */

static doublereal c_b10 = 0.;
static doublereal c_b11 = 1.;
static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int mb03wd_(char *job, char *compz, integer *n, integer *p, 
	integer *ilo, integer *ihi, integer *iloz, integer *ihiz, doublereal *
	h__, integer *ldh1, integer *ldh2, doublereal *z__, integer *ldz1, 
	integer *ldz2, doublereal *wr, doublereal *wi, doublereal *dwork, 
	integer *ldwork, integer *info, ftnlen job_len, ftnlen compz_len)
{
    /* System generated locals */
    integer h_dim1, h_dim2, h_offset, z_dim1, z_dim2, z_offset, i__1, i__2, 
	    i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    static integer i__, j, k, l, m;
    static doublereal s, v[3];
    static integer i1, i2;
    static doublereal v1, v2, v3, h11, h12, h21, h22, h33, h44;
    static integer nh;
    static doublereal cs;
    static integer nr;
    static doublereal sn;
    static integer nz;
    static doublereal hh10, hh11, hh12, hh21, hh22, hp00, hp01, ave, hp02, 
	    hp11, hp12, hp22, h33s, h44s, tau;
    static integer itn, its;
    static doublereal ulp, tst1, h43h34, disc;
    static integer jmin, jmax;
    static doublereal unfl, ovfl;
    extern /* Subroutine */ int drot_(integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, doublereal *);
    static integer nrow;
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int mb04py_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     ftnlen), dcopy_(integer *, doublereal *, integer *, doublereal *,
	     integer *);
    static logical initz, wantt, wantz;
    extern /* Subroutine */ int dlanv2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *), dlabad_(
	    doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int dlarfg_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *);
    extern doublereal dlanhs_(char *, integer *, doublereal *, integer *, 
	    doublereal *, ftnlen);
    extern /* Subroutine */ int dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    dlartg_(doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), xerbla_(char *, integer *, ftnlen), dlarfx_(char *,
	     integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    integer *, doublereal *, ftnlen);
    extern doublereal dlantr_(char *, char *, char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, ftnlen, ftnlen, ftnlen);
    static doublereal smlnum;


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

/*     To compute the Schur decomposition and the eigenvalues of a */
/*     product of matrices, H = H_1*H_2*...*H_p, with H_1 an upper */
/*     Hessenberg matrix and H_2, ..., H_p upper triangular matrices, */
/*     without evaluating the product. Specifically, the matrices Z_i */
/*     are computed, such that */

/*             Z_1' * H_1 * Z_2 = T_1, */
/*             Z_2' * H_2 * Z_3 = T_2, */
/*                    ... */
/*             Z_p' * H_p * Z_1 = T_p, */

/*     where T_1 is in real Schur form, and T_2, ..., T_p are upper */
/*     triangular. */

/*     The routine works primarily with the Hessenberg and triangular */
/*     submatrices in rows and columns ILO to IHI, but optionally applies */
/*     the transformations to all the rows and columns of the matrices */
/*     H_i, i = 1,...,p. The transformations can be optionally */
/*     accumulated. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     JOB     CHARACTER*1 */
/*             Indicates whether the user wishes to compute the full */
/*             Schur form or the eigenvalues only, as follows: */
/*             = 'E':  Compute the eigenvalues only; */
/*             = 'S':  Compute the factors T_1, ..., T_p of the full */
/*                     Schur form, T = T_1*T_2*...*T_p. */

/*     COMPZ   CHARACTER*1 */
/*             Indicates whether or not the user wishes to accumulate */
/*             the matrices Z_1, ..., Z_p, as follows: */
/*             = 'N':  The matrices Z_1, ..., Z_p are not required; */
/*             = 'I':  Z_i is initialized to the unit matrix and the */
/*                     orthogonal transformation matrix Z_i is returned, */
/*                     i = 1, ..., p; */
/*             = 'V':  Z_i must contain an orthogonal matrix Q_i on */
/*                     entry, and the product Q_i*Z_i is returned, */
/*                     i = 1, ..., p. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix H.  N >= 0. */

/*     P       (input) INTEGER */
/*             The number of matrices in the product H_1*H_2*...*H_p. */
/*             P >= 1. */

/*     ILO     (input) INTEGER */
/*     IHI     (input) INTEGER */
/*             It is assumed that all matrices H_j, j = 2, ..., p, are */
/*             already upper triangular in rows and columns 1:ILO-1 and */
/*             IHI+1:N, and H_1 is upper quasi-triangular in rows and */
/*             columns 1:ILO-1 and IHI+1:N, with H_1(ILO,ILO-1) = 0 */
/*             (unless ILO = 1), and H_1(IHI+1,IHI) = 0 (unless IHI = N). */
/*             The routine works primarily with the Hessenberg submatrix */
/*             in rows and columns ILO to IHI, but applies the */
/*             transformations to all the rows and columns of the */
/*             matrices H_i, i = 1,...,p, if JOB = 'S'. */
/*             1 <= ILO <= max(1,N); min(ILO,N) <= IHI <= N. */

/*     ILOZ    (input) INTEGER */
/*     IHIZ    (input) INTEGER */
/*             Specify the rows of Z to which the transformations must be */
/*             applied if COMPZ = 'I' or COMPZ = 'V'. */
/*             1 <= ILOZ <= ILO; IHI <= IHIZ <= N. */

/*     H       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDH1,LDH2,P) */
/*             On entry, the leading N-by-N part of H(*,*,1) must contain */
/*             the upper Hessenberg matrix H_1 and the leading N-by-N */
/*             part of H(*,*,j) for j > 1 must contain the upper */
/*             triangular matrix H_j, j = 2, ..., p. */
/*             On exit, if JOB = 'S', the leading N-by-N part of H(*,*,1) */
/*             is upper quasi-triangular in rows and columns ILO:IHI, */
/*             with any 2-by-2 diagonal blocks corresponding to a pair of */
/*             complex conjugated eigenvalues, and the leading N-by-N */
/*             part of H(*,*,j) for j > 1 contains the resulting upper */
/*             triangular matrix T_j. */
/*             If JOB = 'E', the contents of H are unspecified on exit. */

/*     LDH1    INTEGER */
/*             The first leading dimension of the array H. */
/*             LDH1 >= max(1,N). */

/*     LDH2    INTEGER */
/*             The second leading dimension of the array H. */
/*             LDH2 >= max(1,N). */

/*     Z       (input/output) DOUBLE PRECISION array, dimension */
/*             (LDZ1,LDZ2,P) */
/*             On entry, if COMPZ = 'V', the leading N-by-N-by-P part of */
/*             this array must contain the current matrix Q of */
/*             transformations accumulated by SLICOT Library routine */
/*             MB03VY. */
/*             If COMPZ = 'I', Z need not be set on entry. */
/*             On exit, if COMPZ = 'V', or COMPZ = 'I', the leading */
/*             N-by-N-by-P part of this array contains the transformation */
/*             matrices which produced the Schur form; the */
/*             transformations are applied only to the submatrices */
/*             Z_j(ILOZ:IHIZ,ILO:IHI), j = 1, ..., P. */
/*             If COMPZ = 'N', Z is not referenced. */

/*     LDZ1    INTEGER */
/*             The first leading dimension of the array Z. */
/*             LDZ1 >= 1,        if COMPZ = 'N'; */
/*             LDZ1 >= max(1,N), if COMPZ = 'I' or COMPZ = 'V'. */

/*     LDZ2    INTEGER */
/*             The second leading dimension of the array Z. */
/*             LDZ2 >= 1,        if COMPZ = 'N'; */
/*             LDZ2 >= max(1,N), if COMPZ = 'I' or COMPZ = 'V'. */

/*     WR      (output) DOUBLE PRECISION array, dimension (N) */
/*     WI      (output) DOUBLE PRECISION array, dimension (N) */
/*             The real and imaginary parts, respectively, of the */
/*             computed eigenvalues ILO to IHI are stored in the */
/*             corresponding elements of WR and WI. If two eigenvalues */
/*             are computed as a complex conjugate pair, they are stored */
/*             in consecutive elements of WR and WI, say the i-th and */
/*             (i+1)th, with WI(i) > 0 and WI(i+1) < 0. If JOB = 'S', the */
/*             eigenvalues are stored in the same order as on the */
/*             diagonal of the Schur form returned in H. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION work array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= IHI-ILO+P-1. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             > 0:  if INFO = i, ILO <= i <= IHI, the QR algorithm */
/*                   failed to compute all the eigenvalues ILO to IHI */
/*                   in a total of 30*(IHI-ILO+1) iterations; */
/*                   the elements i+1:IHI of WR and WI contain those */
/*                   eigenvalues which have been successfully computed. */

/*     METHOD */

/*     A refined version of the QR algorithm proposed in [1] and [2] is */
/*     used. The elements of the subdiagonal, diagonal, and the first */
/*     supradiagonal of current principal submatrix of H are computed */
/*     in the process. */

/*     REFERENCES */

/*     [1] Bojanczyk, A.W., Golub, G. and Van Dooren, P. */
/*         The periodic Schur decomposition: algorithms and applications. */
/*         Proc. of the SPIE Conference (F.T. Luk, Ed.), 1770, pp. 31-42, */
/*         1992. */

/*     [2] Sreedhar, J. and Van Dooren, P. */
/*         Periodic Schur form and some matrix equations. */
/*         Proc. of the Symposium on the Mathematical Theory of Networks */
/*         and Systems (MTNS'93), Regensburg, Germany (U. Helmke, */
/*         R. Mennicken and J. Saurer, Eds.), Vol. 1, pp. 339-362, 1994. */

/*     NUMERICAL ASPECTS */

/*     The algorithm is numerically stable. */

/*     FURTHER COMMENTS */

/*     Note that for P = 1, the LAPACK Library routine DHSEQR could be */
/*     more efficient on some computer architectures than this routine, */
/*     because DHSEQR uses a block multishift QR algorithm. */
/*     When P is large and JOB = 'S', it could be more efficient to */
/*     compute the product matrix H, and use the LAPACK Library routines. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, and A. Varga, */
/*     German Aerospace Center, DLR Oberpfaffenhofen, February 1999. */
/*     Partly based on the routine PSHQR by A. Varga */
/*     (DLR Oberpfaffenhofen), January 22, 1996. */

/*     REVISIONS */

/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */

/*     KEYWORDS */

/*     Eigenvalue, eigenvalue decomposition, Hessenberg form, */
/*     orthogonal transformation, periodic systems, (periodic) Schur */
/*     form, real Schur form, similarity transformation, triangular form. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. */
/*     .. Scalar Arguments .. */
/*     .. */
/*     .. Array Arguments .. */
/*     .. */
/*     .. Local Scalars .. */
/*     .. */
/*     .. Local Arrays .. */
/*     .. */
/*     .. External Functions .. */
/*     .. */
/*     .. External Subroutines .. */
/*     .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

/*     Test the input scalar arguments. */

#line 265 "MB03WD.f"
    /* Parameter adjustments */
#line 265 "MB03WD.f"
    h_dim1 = *ldh1;
#line 265 "MB03WD.f"
    h_dim2 = *ldh2;
#line 265 "MB03WD.f"
    h_offset = 1 + h_dim1 * (1 + h_dim2);
#line 265 "MB03WD.f"
    h__ -= h_offset;
#line 265 "MB03WD.f"
    z_dim1 = *ldz1;
#line 265 "MB03WD.f"
    z_dim2 = *ldz2;
#line 265 "MB03WD.f"
    z_offset = 1 + z_dim1 * (1 + z_dim2);
#line 265 "MB03WD.f"
    z__ -= z_offset;
#line 265 "MB03WD.f"
    --wr;
#line 265 "MB03WD.f"
    --wi;
#line 265 "MB03WD.f"
    --dwork;
#line 265 "MB03WD.f"

#line 265 "MB03WD.f"
    /* Function Body */
#line 265 "MB03WD.f"
    wantt = lsame_(job, "S", (ftnlen)1, (ftnlen)1);
#line 266 "MB03WD.f"
    initz = lsame_(compz, "I", (ftnlen)1, (ftnlen)1);
#line 267 "MB03WD.f"
    wantz = lsame_(compz, "V", (ftnlen)1, (ftnlen)1) || initz;
#line 268 "MB03WD.f"
    *info = 0;
#line 269 "MB03WD.f"
    if (! (wantt || lsame_(job, "E", (ftnlen)1, (ftnlen)1))) {
#line 270 "MB03WD.f"
	*info = -1;
#line 271 "MB03WD.f"
    } else if (! (wantz || lsame_(compz, "N", (ftnlen)1, (ftnlen)1))) {
#line 272 "MB03WD.f"
	*info = -2;
#line 273 "MB03WD.f"
    } else if (*n < 0) {
#line 274 "MB03WD.f"
	*info = -3;
#line 275 "MB03WD.f"
    } else if (*p < 1) {
#line 276 "MB03WD.f"
	*info = -4;
#line 277 "MB03WD.f"
    } else if (*ilo < 1 || *ilo > max(1,*n)) {
#line 278 "MB03WD.f"
	*info = -5;
#line 279 "MB03WD.f"
    } else if (*ihi < min(*ilo,*n) || *ihi > *n) {
#line 280 "MB03WD.f"
	*info = -6;
#line 281 "MB03WD.f"
    } else if (*iloz < 1 || *iloz > *ilo) {
#line 282 "MB03WD.f"
	*info = -7;
#line 283 "MB03WD.f"
    } else if (*ihiz < *ihi || *ihiz > *n) {
#line 284 "MB03WD.f"
	*info = -8;
#line 285 "MB03WD.f"
    } else if (*ldh1 < max(1,*n)) {
#line 286 "MB03WD.f"
	*info = -10;
#line 287 "MB03WD.f"
    } else if (*ldh2 < max(1,*n)) {
#line 288 "MB03WD.f"
	*info = -11;
#line 289 "MB03WD.f"
    } else if (*ldz1 < 1 || wantz && *ldz1 < *n) {
#line 290 "MB03WD.f"
	*info = -13;
#line 291 "MB03WD.f"
    } else if (*ldz2 < 1 || wantz && *ldz2 < *n) {
#line 292 "MB03WD.f"
	*info = -14;
#line 293 "MB03WD.f"
    } else if (*ldwork < *ihi - *ilo + *p - 1) {
#line 294 "MB03WD.f"
	*info = -18;
#line 295 "MB03WD.f"
    }
#line 296 "MB03WD.f"
    if (*info == 0) {
#line 297 "MB03WD.f"
	if (*ilo > 1) {
#line 298 "MB03WD.f"
	    if (h__[*ilo + (*ilo - 1 + h_dim2) * h_dim1] != 0.) {
#line 298 "MB03WD.f"
		*info = -5;
#line 298 "MB03WD.f"
	    }
#line 300 "MB03WD.f"
	} else if (*ihi < *n) {
#line 301 "MB03WD.f"
	    if (h__[*ihi + 1 + (*ihi + h_dim2) * h_dim1] != 0.) {
#line 301 "MB03WD.f"
		*info = -6;
#line 301 "MB03WD.f"
	    }
#line 303 "MB03WD.f"
	}
#line 304 "MB03WD.f"
    }

#line 306 "MB03WD.f"
    if (*info != 0) {
#line 307 "MB03WD.f"
	i__1 = -(*info);
#line 307 "MB03WD.f"
	xerbla_("MB03WD", &i__1, (ftnlen)6);
#line 308 "MB03WD.f"
	return 0;
#line 309 "MB03WD.f"
    }

/*     Quick return if possible. */

#line 313 "MB03WD.f"
    if (*n == 0) {
#line 313 "MB03WD.f"
	return 0;
#line 313 "MB03WD.f"
    }

/*     Initialize Z, if necessary. */

#line 318 "MB03WD.f"
    if (initz) {

#line 320 "MB03WD.f"
	i__1 = *p;
#line 320 "MB03WD.f"
	for (j = 1; j <= i__1; ++j) {
#line 321 "MB03WD.f"
	    dlaset_("Full", n, n, &c_b10, &c_b11, &z__[(j * z_dim2 + 1) * 
		    z_dim1 + 1], ldz1, (ftnlen)4);
#line 322 "MB03WD.f"
/* L10: */
#line 322 "MB03WD.f"
	}

#line 324 "MB03WD.f"
    }

#line 326 "MB03WD.f"
    nh = *ihi - *ilo + 1;

#line 328 "MB03WD.f"
    if (nh == 1) {
#line 329 "MB03WD.f"
	hp00 = 1.;

#line 331 "MB03WD.f"
	i__1 = *p;
#line 331 "MB03WD.f"
	for (j = 1; j <= i__1; ++j) {
#line 332 "MB03WD.f"
	    hp00 *= h__[*ilo + (*ilo + j * h_dim2) * h_dim1];
#line 333 "MB03WD.f"
/* L20: */
#line 333 "MB03WD.f"
	}

#line 335 "MB03WD.f"
	wr[*ilo] = hp00;
#line 336 "MB03WD.f"
	wi[*ilo] = 0.;
#line 337 "MB03WD.f"
	return 0;
#line 338 "MB03WD.f"
    }

/*     Set machine-dependent constants for the stopping criterion. */
/*     If norm(H) <= sqrt(OVFL), overflow should not occur. */

#line 343 "MB03WD.f"
    unfl = dlamch_("Safe minimum", (ftnlen)12);
#line 344 "MB03WD.f"
    ovfl = 1. / unfl;
#line 345 "MB03WD.f"
    dlabad_(&unfl, &ovfl);
#line 346 "MB03WD.f"
    ulp = dlamch_("Precision", (ftnlen)9);
#line 347 "MB03WD.f"
    smlnum = unfl * ((doublereal) nh / ulp);

/*     Set the elements in rows and columns ILO to IHI to zero below the */
/*     first subdiagonal in H(*,*,1) and below the first diagonal in */
/*     H(*,*,j), j >= 2. In the same loop, compute and store in */
/*     DWORK(NH:NH+P-2) the 1-norms of the matrices H_2, ..., H_p, to be */
/*     used later. */

#line 355 "MB03WD.f"
    i__ = nh;
#line 356 "MB03WD.f"
    s = ulp * (doublereal) (*n);
#line 357 "MB03WD.f"
    if (nh > 2) {
#line 357 "MB03WD.f"
	i__1 = nh - 2;
#line 357 "MB03WD.f"
	i__2 = nh - 2;
#line 357 "MB03WD.f"
	dlaset_("Lower", &i__1, &i__2, &c_b10, &c_b10, &h__[*ilo + 2 + (*ilo 
		+ h_dim2) * h_dim1], ldh1, (ftnlen)5);
#line 357 "MB03WD.f"
    }

#line 361 "MB03WD.f"
    i__1 = *p;
#line 361 "MB03WD.f"
    for (j = 2; j <= i__1; ++j) {
#line 362 "MB03WD.f"
	i__2 = nh - 1;
#line 362 "MB03WD.f"
	i__3 = nh - 1;
#line 362 "MB03WD.f"
	dlaset_("Lower", &i__2, &i__3, &c_b10, &c_b10, &h__[*ilo + 1 + (*ilo 
		+ j * h_dim2) * h_dim1], ldh1, (ftnlen)5);
#line 364 "MB03WD.f"
	dwork[i__] = s * dlantr_("1-norm", "Upper", "NonUnit", &nh, &nh, &h__[
		*ilo + (*ilo + j * h_dim2) * h_dim1], ldh1, &dwork[1], (
		ftnlen)6, (ftnlen)5, (ftnlen)7);
#line 366 "MB03WD.f"
	++i__;
#line 367 "MB03WD.f"
/* L30: */
#line 367 "MB03WD.f"
    }

/*     I1 and I2 are the indices of the first row and last column of H */
/*     to which transformations must be applied. If eigenvalues only are */
/*     being computed, I1 and I2 are set inside the main loop. */

#line 373 "MB03WD.f"
    if (wantt) {
#line 374 "MB03WD.f"
	i1 = 1;
#line 375 "MB03WD.f"
	i2 = *n;
#line 376 "MB03WD.f"
    }

#line 378 "MB03WD.f"
    if (wantz) {
#line 378 "MB03WD.f"
	nz = *ihiz - *iloz + 1;
#line 378 "MB03WD.f"
    }

/*     ITN is the total number of QR iterations allowed. */

#line 383 "MB03WD.f"
    itn = nh * 30;

/*     The main loop begins here. I is the loop index and decreases from */
/*     IHI to ILO in steps of 1 or 2. Each iteration of the loop works */
/*     with the active submatrix in rows and columns L to I. */
/*     Eigenvalues I+1 to IHI have already converged. Either L = ILO or */
/*     H(L,L-1) is negligible so that the matrix splits. */

#line 391 "MB03WD.f"
    i__ = *ihi;

#line 393 "MB03WD.f"
L40:
#line 394 "MB03WD.f"
    l = *ilo;

/*     Perform QR iterations on rows and columns ILO to I until a */
/*     submatrix of order 1 or 2 splits off at the bottom because a */
/*     subdiagonal element has become negligible. */

/*     Let T = H_2*...*H_p, and H = H_1*T. Part of the currently */
/*     free locations of WR and WI are temporarily used as workspace. */

/*     WR(L:I):      the current diagonal elements of h = H(L:I,L:I); */
/*     WI(L+1:I):    the current elements of the first subdiagonal of h; */
/*     DWORK(NH-I+L:NH-1): the current elements of the first */
/*                   supradiagonal of h. */

#line 408 "MB03WD.f"
    i__1 = itn;
#line 408 "MB03WD.f"
    for (its = 0; its <= i__1; ++its) {

/*        Initialization: compute H(I,I) (and H(I,I-1) if I > L). */

#line 412 "MB03WD.f"
	hp22 = 1.;
#line 413 "MB03WD.f"
	if (i__ > l) {
#line 414 "MB03WD.f"
	    hp12 = 0.;
#line 415 "MB03WD.f"
	    hp11 = 1.;

#line 417 "MB03WD.f"
	    i__2 = *p;
#line 417 "MB03WD.f"
	    for (j = 2; j <= i__2; ++j) {
#line 418 "MB03WD.f"
		hp22 *= h__[i__ + (i__ + j * h_dim2) * h_dim1];
#line 419 "MB03WD.f"
		hp12 = hp11 * h__[i__ - 1 + (i__ + j * h_dim2) * h_dim1] + 
			hp12 * h__[i__ + (i__ + j * h_dim2) * h_dim1];
#line 420 "MB03WD.f"
		hp11 *= h__[i__ - 1 + (i__ - 1 + j * h_dim2) * h_dim1];
#line 421 "MB03WD.f"
/* L50: */
#line 421 "MB03WD.f"
	    }

#line 423 "MB03WD.f"
	    hh21 = h__[i__ + (i__ - 1 + h_dim2) * h_dim1] * hp11;
#line 424 "MB03WD.f"
	    hh22 = h__[i__ + (i__ - 1 + h_dim2) * h_dim1] * hp12 + h__[i__ + (
		    i__ + h_dim2) * h_dim1] * hp22;

#line 426 "MB03WD.f"
	    wr[i__] = hh22;
#line 427 "MB03WD.f"
	    wi[i__] = hh21;
#line 428 "MB03WD.f"
	} else {

#line 430 "MB03WD.f"
	    i__2 = *p;
#line 430 "MB03WD.f"
	    for (j = 1; j <= i__2; ++j) {
#line 431 "MB03WD.f"
		hp22 *= h__[i__ + (i__ + j * h_dim2) * h_dim1];
#line 432 "MB03WD.f"
/* L60: */
#line 432 "MB03WD.f"
	    }

#line 434 "MB03WD.f"
	    wr[i__] = hp22;
#line 435 "MB03WD.f"
	}

/*        Look for a single small subdiagonal element. */
/*        The loop also computes the needed current elements of the */
/*        diagonal and the first two supradiagonals of T, as well as */
/*        the current elements of the central tridiagonal of H. */

#line 442 "MB03WD.f"
	i__2 = l + 1;
#line 442 "MB03WD.f"
	for (k = i__; k >= i__2; --k) {

/*           Evaluate H(K-1,K-1), H(K-1,K) (and H(K-1,K-2) if K > L+1). */

#line 446 "MB03WD.f"
	    hp00 = 1.;
#line 447 "MB03WD.f"
	    hp01 = 0.;
#line 448 "MB03WD.f"
	    if (k > l + 1) {
#line 449 "MB03WD.f"
		hp02 = 0.;

#line 451 "MB03WD.f"
		i__3 = *p;
#line 451 "MB03WD.f"
		for (j = 2; j <= i__3; ++j) {
#line 452 "MB03WD.f"
		    hp02 = hp00 * h__[k - 2 + (k + j * h_dim2) * h_dim1] + 
			    hp01 * h__[k - 1 + (k + j * h_dim2) * h_dim1] + 
			    hp02 * h__[k + (k + j * h_dim2) * h_dim1];
#line 454 "MB03WD.f"
		    hp01 = hp00 * h__[k - 2 + (k - 1 + j * h_dim2) * h_dim1] 
			    + hp01 * h__[k - 1 + (k - 1 + j * h_dim2) * 
			    h_dim1];
#line 455 "MB03WD.f"
		    hp00 *= h__[k - 2 + (k - 2 + j * h_dim2) * h_dim1];
#line 456 "MB03WD.f"
/* L70: */
#line 456 "MB03WD.f"
		}

#line 458 "MB03WD.f"
		hh10 = h__[k - 1 + (k - 2 + h_dim2) * h_dim1] * hp00;
#line 459 "MB03WD.f"
		hh11 = h__[k - 1 + (k - 2 + h_dim2) * h_dim1] * hp01 + h__[k 
			- 1 + (k - 1 + h_dim2) * h_dim1] * hp11;
#line 460 "MB03WD.f"
		hh12 = h__[k - 1 + (k - 2 + h_dim2) * h_dim1] * hp02 + h__[k 
			- 1 + (k - 1 + h_dim2) * h_dim1] * hp12 + h__[k - 1 + 
			(k + h_dim2) * h_dim1] * hp22;
#line 462 "MB03WD.f"
		wi[k - 1] = hh10;
#line 463 "MB03WD.f"
	    } else {
#line 464 "MB03WD.f"
		hh10 = 0.;
#line 465 "MB03WD.f"
		hh11 = h__[k - 1 + (k - 1 + h_dim2) * h_dim1] * hp11;
#line 466 "MB03WD.f"
		hh12 = h__[k - 1 + (k - 1 + h_dim2) * h_dim1] * hp12 + h__[k 
			- 1 + (k + h_dim2) * h_dim1] * hp22;
#line 467 "MB03WD.f"
	    }
#line 468 "MB03WD.f"
	    wr[k - 1] = hh11;
#line 469 "MB03WD.f"
	    dwork[nh - i__ + k - 1] = hh12;

/*           Test for a negligible subdiagonal element. */

#line 473 "MB03WD.f"
	    tst1 = abs(hh11) + abs(hh22);
#line 474 "MB03WD.f"
	    if (tst1 == 0.) {
#line 474 "MB03WD.f"
		i__3 = i__ - l + 1;
#line 474 "MB03WD.f"
		tst1 = dlanhs_("1-norm", &i__3, &h__[l + (l + h_dim2) * 
			h_dim1], ldh1, &dwork[1], (ftnlen)6);
#line 474 "MB03WD.f"
	    }
/* Computing MAX */
#line 477 "MB03WD.f"
	    d__1 = ulp * tst1;
#line 477 "MB03WD.f"
	    if (abs(hh21) <= max(d__1,smlnum)) {
#line 477 "MB03WD.f"
		goto L90;
#line 477 "MB03WD.f"
	    }

/*           Update the values for the next cycle. */

#line 482 "MB03WD.f"
	    hp22 = hp11;
#line 483 "MB03WD.f"
	    hp11 = hp00;
#line 484 "MB03WD.f"
	    hp12 = hp01;
#line 485 "MB03WD.f"
	    hh22 = hh11;
#line 486 "MB03WD.f"
	    hh21 = hh10;
#line 487 "MB03WD.f"
/* L80: */
#line 487 "MB03WD.f"
	}

#line 489 "MB03WD.f"
L90:
#line 490 "MB03WD.f"
	l = k;

#line 492 "MB03WD.f"
	if (l > *ilo) {

/*           H(L,L-1) is negligible. */

#line 496 "MB03WD.f"
	    if (wantt) {

/*              If H(L,L-1,1) is also negligible, set it to 0; otherwise, */
/*              annihilate the subdiagonal elements bottom-up, and */
/*              restore the triangular form of H(*,*,j). Since H(L,L-1) */
/*              is negligible, the second case can only appear when the */
/*              product of H(L-1,L-1,j), j >= 2, is negligible. */

#line 504 "MB03WD.f"
		tst1 = (d__1 = h__[l - 1 + (l - 1 + h_dim2) * h_dim1], abs(
			d__1)) + (d__2 = h__[l + (l + h_dim2) * h_dim1], abs(
			d__2));
#line 505 "MB03WD.f"
		if (tst1 == 0.) {
#line 505 "MB03WD.f"
		    i__2 = i__ - l + 1;
#line 505 "MB03WD.f"
		    tst1 = dlanhs_("1-norm", &i__2, &h__[l + (l + h_dim2) * 
			    h_dim1], ldh1, &dwork[1], (ftnlen)6);
#line 505 "MB03WD.f"
		}
/* Computing MAX */
#line 508 "MB03WD.f"
		d__2 = ulp * tst1;
#line 508 "MB03WD.f"
		if ((d__1 = h__[l + (l - 1 + h_dim2) * h_dim1], abs(d__1)) > 
			max(d__2,smlnum)) {

#line 511 "MB03WD.f"
		    i__2 = l;
#line 511 "MB03WD.f"
		    for (k = i__; k >= i__2; --k) {

#line 513 "MB03WD.f"
			i__3 = *p - 1;
#line 513 "MB03WD.f"
			for (j = 1; j <= i__3; ++j) {

/*                       Compute G to annihilate from the right the */
/*                       (K,K-1) element of the matrix H_j. */

#line 518 "MB03WD.f"
			    v[0] = h__[k + (k - 1 + j * h_dim2) * h_dim1];
#line 519 "MB03WD.f"
			    dlarfg_(&c__2, &h__[k + (k + j * h_dim2) * h_dim1]
				    , v, &c__1, &tau);
#line 520 "MB03WD.f"
			    h__[k + (k - 1 + j * h_dim2) * h_dim1] = 0.;
#line 521 "MB03WD.f"
			    v[1] = 1.;

/*                       Apply G from the right to transform the columns */
/*                       of the matrix H_j in rows I1 to K-1. */

#line 526 "MB03WD.f"
			    i__4 = k - i1;
#line 526 "MB03WD.f"
			    dlarfx_("Right", &i__4, &c__2, v, &tau, &h__[i1 + 
				    (k - 1 + j * h_dim2) * h_dim1], ldh1, &
				    dwork[1], (ftnlen)5);

/*                       Apply G from the left to transform the rows of */
/*                       the matrix H_(j+1) in columns K-1 to I2. */

#line 532 "MB03WD.f"
			    i__4 = i2 - k + 2;
#line 532 "MB03WD.f"
			    dlarfx_("Left", &c__2, &i__4, v, &tau, &h__[k - 1 
				    + (k - 1 + (j + 1) * h_dim2) * h_dim1], 
				    ldh1, &dwork[1], (ftnlen)4);

#line 535 "MB03WD.f"
			    if (wantz) {

/*                          Accumulate transformations in the matrix */
/*                          Z_(j+1). */

#line 540 "MB03WD.f"
				dlarfx_("Right", &nz, &c__2, v, &tau, &z__[*
					iloz + (k - 1 + (j + 1) * z_dim2) * 
					z_dim1], ldz1, &dwork[1], (ftnlen)5);
#line 543 "MB03WD.f"
			    }
#line 544 "MB03WD.f"
/* L100: */
#line 544 "MB03WD.f"
			}

#line 546 "MB03WD.f"
			if (k < i__) {

/*                       Compute G to annihilate from the right the */
/*                       (K+1,K) element of the matrix H_p. */

#line 551 "MB03WD.f"
			    v[0] = h__[k + 1 + (k + *p * h_dim2) * h_dim1];
#line 552 "MB03WD.f"
			    dlarfg_(&c__2, &h__[k + 1 + (k + 1 + *p * h_dim2) 
				    * h_dim1], v, &c__1, &tau);
#line 553 "MB03WD.f"
			    h__[k + 1 + (k + *p * h_dim2) * h_dim1] = 0.;
#line 554 "MB03WD.f"
			    v[1] = 1.;

/*                       Apply G from the right to transform the columns */
/*                       of the matrix H_p in rows I1 to K. */

#line 559 "MB03WD.f"
			    i__3 = k - i1 + 1;
#line 559 "MB03WD.f"
			    dlarfx_("Right", &i__3, &c__2, v, &tau, &h__[i1 + 
				    (k + *p * h_dim2) * h_dim1], ldh1, &dwork[
				    1], (ftnlen)5);

/*                       Apply G from the left to transform the rows of */
/*                       the matrix H_1 in columns K to I2. */

#line 565 "MB03WD.f"
			    i__3 = i2 - k + 1;
#line 565 "MB03WD.f"
			    dlarfx_("Left", &c__2, &i__3, v, &tau, &h__[k + (
				    k + h_dim2) * h_dim1], ldh1, &dwork[1], (
				    ftnlen)4);

#line 568 "MB03WD.f"
			    if (wantz) {

/*                          Accumulate transformations in the matrix Z_1. */

#line 572 "MB03WD.f"
				dlarfx_("Right", &nz, &c__2, v, &tau, &z__[*
					iloz + (k + z_dim2) * z_dim1], ldz1, &
					dwork[1], (ftnlen)5);
#line 574 "MB03WD.f"
			    }
#line 575 "MB03WD.f"
			}
#line 576 "MB03WD.f"
/* L110: */
#line 576 "MB03WD.f"
		    }

#line 578 "MB03WD.f"
		    h__[l + (l - 1 + *p * h_dim2) * h_dim1] = 0.;
#line 579 "MB03WD.f"
		}
#line 580 "MB03WD.f"
		h__[l + (l - 1 + h_dim2) * h_dim1] = 0.;
#line 581 "MB03WD.f"
	    }
#line 582 "MB03WD.f"
	}

/*        Exit from loop if a submatrix of order 1 or 2 has split off. */

#line 586 "MB03WD.f"
	if (l >= i__ - 1) {
#line 586 "MB03WD.f"
	    goto L170;
#line 586 "MB03WD.f"
	}

/*        Now the active submatrix is in rows and columns L to I. If */
/*        eigenvalues only are being computed, only the active submatrix */
/*        need be transformed. */

#line 593 "MB03WD.f"
	if (! wantt) {
#line 594 "MB03WD.f"
	    i1 = l;
#line 595 "MB03WD.f"
	    i2 = i__;
#line 596 "MB03WD.f"
	}

#line 598 "MB03WD.f"
	if (its == 10 || its == 20) {

/*           Exceptional shift. */

#line 602 "MB03WD.f"
	    s = (d__1 = wi[i__], abs(d__1)) + (d__2 = wi[i__ - 1], abs(d__2));
#line 603 "MB03WD.f"
	    h44 = s * .75 + wr[i__];
#line 604 "MB03WD.f"
	    h33 = h44;
#line 605 "MB03WD.f"
	    h43h34 = s * -.4375 * s;
#line 606 "MB03WD.f"
	} else {

/*           Prepare to use Francis' double shift (i.e., second degree */
/*           generalized Rayleigh quotient). */

#line 611 "MB03WD.f"
	    h44 = wr[i__];
#line 612 "MB03WD.f"
	    h33 = wr[i__ - 1];
#line 613 "MB03WD.f"
	    h43h34 = wi[i__] * dwork[nh - 1];
#line 614 "MB03WD.f"
	    disc = (h33 - h44) * .5;
#line 615 "MB03WD.f"
	    disc = disc * disc + h43h34;
#line 616 "MB03WD.f"
	    if (disc > 0.) {

/*              Real roots: use Wilkinson's shift twice. */

#line 620 "MB03WD.f"
		disc = sqrt(disc);
#line 621 "MB03WD.f"
		ave = (h33 + h44) * .5;
#line 622 "MB03WD.f"
		if (abs(h33) - abs(h44) > 0.) {
#line 623 "MB03WD.f"
		    h33 = h33 * h44 - h43h34;
#line 624 "MB03WD.f"
		    h44 = h33 / (d_sign(&disc, &ave) + ave);
#line 625 "MB03WD.f"
		} else {
#line 626 "MB03WD.f"
		    h44 = d_sign(&disc, &ave) + ave;
#line 627 "MB03WD.f"
		}
#line 628 "MB03WD.f"
		h33 = h44;
#line 629 "MB03WD.f"
		h43h34 = 0.;
#line 630 "MB03WD.f"
	    }
#line 631 "MB03WD.f"
	}

/*        Look for two consecutive small subdiagonal elements. */

#line 635 "MB03WD.f"
	i__2 = l;
#line 635 "MB03WD.f"
	for (m = i__ - 2; m >= i__2; --m) {

/*           Determine the effect of starting the double-shift QR */
/*           iteration at row M, and see if this would make H(M,M-1) */
/*           negligible. */

#line 641 "MB03WD.f"
	    h11 = wr[m];
#line 642 "MB03WD.f"
	    h12 = dwork[nh - i__ + m];
#line 643 "MB03WD.f"
	    h21 = wi[m + 1];
#line 644 "MB03WD.f"
	    h22 = wr[m + 1];
#line 645 "MB03WD.f"
	    h44s = h44 - h11;
#line 646 "MB03WD.f"
	    h33s = h33 - h11;
#line 647 "MB03WD.f"
	    v1 = (h33s * h44s - h43h34) / h21 + h12;
#line 648 "MB03WD.f"
	    v2 = h22 - h11 - h33s - h44s;
#line 649 "MB03WD.f"
	    v3 = wi[m + 2];
#line 650 "MB03WD.f"
	    s = abs(v1) + abs(v2) + abs(v3);
#line 651 "MB03WD.f"
	    v1 /= s;
#line 652 "MB03WD.f"
	    v2 /= s;
#line 653 "MB03WD.f"
	    v3 /= s;
#line 654 "MB03WD.f"
	    v[0] = v1;
#line 655 "MB03WD.f"
	    v[1] = v2;
#line 656 "MB03WD.f"
	    v[2] = v3;
#line 657 "MB03WD.f"
	    if (m == l) {
#line 657 "MB03WD.f"
		goto L130;
#line 657 "MB03WD.f"
	    }
#line 659 "MB03WD.f"
	    tst1 = abs(v1) * ((d__1 = wr[m - 1], abs(d__1)) + abs(h11) + abs(
		    h22));
#line 661 "MB03WD.f"
	    if ((d__1 = wi[m], abs(d__1)) * (abs(v2) + abs(v3)) <= ulp * tst1)
		     {
#line 661 "MB03WD.f"
		goto L130;
#line 661 "MB03WD.f"
	    }
#line 663 "MB03WD.f"
/* L120: */
#line 663 "MB03WD.f"
	}

#line 665 "MB03WD.f"
L130:

/*        Double-shift QR step. */

#line 669 "MB03WD.f"
	i__2 = i__ - 1;
#line 669 "MB03WD.f"
	for (k = m; k <= i__2; ++k) {

/*           The first iteration of this loop determines a reflection G */
/*           from the vector V and applies it from left and right to H, */
/*           thus creating a nonzero bulge below the subdiagonal. */

/*           Each subsequent iteration determines a reflection G to */
/*           restore the Hessenberg form in the (K-1)th column, and thus */
/*           chases the bulge one step toward the bottom of the active */
/*           submatrix. NR is the order of G. */

/* Computing MIN */
#line 680 "MB03WD.f"
	    i__3 = 3, i__4 = i__ - k + 1;
#line 680 "MB03WD.f"
	    nr = min(i__3,i__4);
/* Computing MIN */
#line 681 "MB03WD.f"
	    i__3 = k + nr;
#line 681 "MB03WD.f"
	    nrow = min(i__3,i__) - i1 + 1;
#line 682 "MB03WD.f"
	    if (k > m) {
#line 682 "MB03WD.f"
		dcopy_(&nr, &h__[k + (k - 1 + h_dim2) * h_dim1], &c__1, v, &
			c__1);
#line 682 "MB03WD.f"
	    }
#line 684 "MB03WD.f"
	    dlarfg_(&nr, v, &v[1], &c__1, &tau);
#line 685 "MB03WD.f"
	    if (k > m) {
#line 686 "MB03WD.f"
		h__[k + (k - 1 + h_dim2) * h_dim1] = v[0];
#line 687 "MB03WD.f"
		h__[k + 1 + (k - 1 + h_dim2) * h_dim1] = 0.;
#line 688 "MB03WD.f"
		if (k < i__ - 1) {
#line 688 "MB03WD.f"
		    h__[k + 2 + (k - 1 + h_dim2) * h_dim1] = 0.;
#line 688 "MB03WD.f"
		}
#line 690 "MB03WD.f"
	    } else if (m > l) {
#line 691 "MB03WD.f"
		h__[k + (k - 1 + h_dim2) * h_dim1] = -h__[k + (k - 1 + h_dim2)
			 * h_dim1];
#line 692 "MB03WD.f"
	    }

/*           Apply G from the left to transform the rows of the matrix */
/*           H_1 in columns K to I2. */

#line 697 "MB03WD.f"
	    i__3 = i2 - k + 1;
#line 697 "MB03WD.f"
	    mb04py_("Left", &nr, &i__3, &v[1], &tau, &h__[k + (k + h_dim2) * 
		    h_dim1], ldh1, &dwork[1], (ftnlen)4);

/*           Apply G from the right to transform the columns of the */
/*           matrix H_p in rows I1 to min(K+NR,I). */

#line 703 "MB03WD.f"
	    mb04py_("Right", &nrow, &nr, &v[1], &tau, &h__[i1 + (k + *p * 
		    h_dim2) * h_dim1], ldh1, &dwork[1], (ftnlen)5);

#line 706 "MB03WD.f"
	    if (wantz) {

/*              Accumulate transformations in the matrix Z_1. */

#line 710 "MB03WD.f"
		mb04py_("Right", &nz, &nr, &v[1], &tau, &z__[*iloz + (k + 
			z_dim2) * z_dim1], ldz1, &dwork[1], (ftnlen)5);
#line 712 "MB03WD.f"
	    }

#line 714 "MB03WD.f"
	    for (j = *p; j >= 2; --j) {

/*              Apply G1 (and G2, if NR = 3) from the left to transform */
/*              the NR-by-NR submatrix of H_j in position (K,K) to upper */
/*              triangular form. */

/*              Compute G1. */

#line 722 "MB03WD.f"
		i__3 = nr - 1;
#line 722 "MB03WD.f"
		dcopy_(&i__3, &h__[k + 1 + (k + j * h_dim2) * h_dim1], &c__1, 
			v, &c__1);
#line 723 "MB03WD.f"
		dlarfg_(&nr, &h__[k + (k + j * h_dim2) * h_dim1], v, &c__1, &
			tau);
#line 724 "MB03WD.f"
		h__[k + 1 + (k + j * h_dim2) * h_dim1] = 0.;
#line 725 "MB03WD.f"
		if (nr == 3) {
#line 725 "MB03WD.f"
		    h__[k + 2 + (k + j * h_dim2) * h_dim1] = 0.;
#line 725 "MB03WD.f"
		}

/*              Apply G1 from the left to transform the rows of the */
/*              matrix H_j in columns K+1 to I2. */

#line 731 "MB03WD.f"
		i__3 = i2 - k;
#line 731 "MB03WD.f"
		mb04py_("Left", &nr, &i__3, v, &tau, &h__[k + (k + 1 + j * 
			h_dim2) * h_dim1], ldh1, &dwork[1], (ftnlen)4);

/*              Apply G1 from the right to transform the columns of the */
/*              matrix H_(j-1) in rows I1 to min(K+NR,I). */

#line 737 "MB03WD.f"
		mb04py_("Right", &nrow, &nr, v, &tau, &h__[i1 + (k + (j - 1) *
			 h_dim2) * h_dim1], ldh1, &dwork[1], (ftnlen)5);

#line 740 "MB03WD.f"
		if (wantz) {

/*                 Accumulate transformations in the matrix Z_j. */

#line 744 "MB03WD.f"
		    mb04py_("Right", &nz, &nr, v, &tau, &z__[*iloz + (k + j * 
			    z_dim2) * z_dim1], ldz1, &dwork[1], (ftnlen)5);
#line 746 "MB03WD.f"
		}

#line 748 "MB03WD.f"
		if (nr == 3) {

/*                 Compute G2. */

#line 752 "MB03WD.f"
		    v[0] = h__[k + 2 + (k + 1 + j * h_dim2) * h_dim1];
#line 753 "MB03WD.f"
		    dlarfg_(&c__2, &h__[k + 1 + (k + 1 + j * h_dim2) * h_dim1]
			    , v, &c__1, &tau);
#line 754 "MB03WD.f"
		    h__[k + 2 + (k + 1 + j * h_dim2) * h_dim1] = 0.;

/*                 Apply G2 from the left to transform the rows of the */
/*                 matrix H_j in columns K+2 to I2. */

#line 759 "MB03WD.f"
		    i__3 = i2 - k - 1;
#line 759 "MB03WD.f"
		    mb04py_("Left", &c__2, &i__3, v, &tau, &h__[k + 1 + (k + 
			    2 + j * h_dim2) * h_dim1], ldh1, &dwork[1], (
			    ftnlen)4);

/*                 Apply G2 from the right to transform the columns of */
/*                 the matrix H_(j-1) in rows I1 to min(K+3,I). */

#line 765 "MB03WD.f"
		    mb04py_("Right", &nrow, &c__2, v, &tau, &h__[i1 + (k + 1 
			    + (j - 1) * h_dim2) * h_dim1], ldh1, &dwork[1], (
			    ftnlen)5);

#line 768 "MB03WD.f"
		    if (wantz) {

/*                    Accumulate transformations in the matrix Z_j. */

#line 772 "MB03WD.f"
			mb04py_("Right", &nz, &c__2, v, &tau, &z__[*iloz + (k 
				+ 1 + j * z_dim2) * z_dim1], ldz1, &dwork[1], 
				(ftnlen)5);
#line 774 "MB03WD.f"
		    }
#line 775 "MB03WD.f"
		}
#line 776 "MB03WD.f"
/* L140: */
#line 776 "MB03WD.f"
	    }

#line 778 "MB03WD.f"
/* L150: */
#line 778 "MB03WD.f"
	}

#line 780 "MB03WD.f"
/* L160: */
#line 780 "MB03WD.f"
    }

/*     Failure to converge in remaining number of iterations. */

#line 784 "MB03WD.f"
    *info = i__;
#line 785 "MB03WD.f"
    return 0;

#line 787 "MB03WD.f"
L170:

#line 789 "MB03WD.f"
    if (l == i__) {

/*        H(I,I-1,1) is negligible: one eigenvalue has converged. */
/*        Note that WR(I) has already been set. */

#line 794 "MB03WD.f"
	wi[i__] = 0.;
#line 795 "MB03WD.f"
    } else if (l == i__ - 1) {

/*        H(I-1,I-2,1) is negligible: a pair of eigenvalues have */
/*        converged. */

/*        Transform the 2-by-2 submatrix of H_1*H_2*...*H_p in position */
/*        (I-1,I-1) to standard Schur form, and compute and store its */
/*        eigenvalues. If the Schur form is not required, then the */
/*        previously stored values of a similar submatrix are used. */
/*        For real eigenvalues, a Givens transformation is used to */
/*        triangularize the submatrix. */

#line 807 "MB03WD.f"
	if (wantt) {
#line 808 "MB03WD.f"
	    hp22 = 1.;
#line 809 "MB03WD.f"
	    hp12 = 0.;
#line 810 "MB03WD.f"
	    hp11 = 1.;

#line 812 "MB03WD.f"
	    i__1 = *p;
#line 812 "MB03WD.f"
	    for (j = 2; j <= i__1; ++j) {
#line 813 "MB03WD.f"
		hp22 *= h__[i__ + (i__ + j * h_dim2) * h_dim1];
#line 814 "MB03WD.f"
		hp12 = hp11 * h__[i__ - 1 + (i__ + j * h_dim2) * h_dim1] + 
			hp12 * h__[i__ + (i__ + j * h_dim2) * h_dim1];
#line 815 "MB03WD.f"
		hp11 *= h__[i__ - 1 + (i__ - 1 + j * h_dim2) * h_dim1];
#line 816 "MB03WD.f"
/* L180: */
#line 816 "MB03WD.f"
	    }

#line 818 "MB03WD.f"
	    hh21 = h__[i__ + (i__ - 1 + h_dim2) * h_dim1] * hp11;
#line 819 "MB03WD.f"
	    hh22 = h__[i__ + (i__ - 1 + h_dim2) * h_dim1] * hp12 + h__[i__ + (
		    i__ + h_dim2) * h_dim1] * hp22;
#line 820 "MB03WD.f"
	    hh11 = h__[i__ - 1 + (i__ - 1 + h_dim2) * h_dim1] * hp11;
#line 821 "MB03WD.f"
	    hh12 = h__[i__ - 1 + (i__ - 1 + h_dim2) * h_dim1] * hp12 + h__[
		    i__ - 1 + (i__ + h_dim2) * h_dim1] * hp22;
#line 822 "MB03WD.f"
	} else {
#line 823 "MB03WD.f"
	    hh11 = wr[i__ - 1];
#line 824 "MB03WD.f"
	    hh12 = dwork[nh - 1];
#line 825 "MB03WD.f"
	    hh21 = wi[i__];
#line 826 "MB03WD.f"
	    hh22 = wr[i__];
#line 827 "MB03WD.f"
	}

#line 829 "MB03WD.f"
	dlanv2_(&hh11, &hh12, &hh21, &hh22, &wr[i__ - 1], &wi[i__ - 1], &wr[
		i__], &wi[i__], &cs, &sn);

#line 832 "MB03WD.f"
	if (wantt) {

/*           Detect negligible diagonal elements in positions (I-1,I-1) */
/*           and (I,I) in H_j, J > 1. */

#line 837 "MB03WD.f"
	    jmin = 0;
#line 838 "MB03WD.f"
	    jmax = 0;

#line 840 "MB03WD.f"
	    i__1 = *p;
#line 840 "MB03WD.f"
	    for (j = 2; j <= i__1; ++j) {
#line 841 "MB03WD.f"
		if (jmin == 0) {
#line 842 "MB03WD.f"
		    if ((d__1 = h__[i__ - 1 + (i__ - 1 + j * h_dim2) * h_dim1]
			    , abs(d__1)) <= dwork[nh + j - 2]) {
#line 842 "MB03WD.f"
			jmin = j;
#line 842 "MB03WD.f"
		    }
#line 844 "MB03WD.f"
		}
#line 845 "MB03WD.f"
		if ((d__1 = h__[i__ + (i__ + j * h_dim2) * h_dim1], abs(d__1))
			 <= dwork[nh + j - 2]) {
#line 845 "MB03WD.f"
		    jmax = j;
#line 845 "MB03WD.f"
		}
#line 846 "MB03WD.f"
/* L190: */
#line 846 "MB03WD.f"
	    }

#line 848 "MB03WD.f"
	    if (jmin != 0 && jmax != 0) {

/*              Choose the shorter path if zero elements in both */
/*              (I-1,I-1) and (I,I) positions are present. */

#line 853 "MB03WD.f"
		if (jmin - 1 <= *p - jmax + 1) {
#line 854 "MB03WD.f"
		    jmax = 0;
#line 855 "MB03WD.f"
		} else {
#line 856 "MB03WD.f"
		    jmin = 0;
#line 857 "MB03WD.f"
		}
#line 858 "MB03WD.f"
	    }

#line 860 "MB03WD.f"
	    if (jmin != 0) {

#line 862 "MB03WD.f"
		i__1 = jmin - 1;
#line 862 "MB03WD.f"
		for (j = 1; j <= i__1; ++j) {

/*                 Compute G to annihilate from the right the (I,I-1) */
/*                 element of the matrix H_j. */

#line 867 "MB03WD.f"
		    v[0] = h__[i__ + (i__ - 1 + j * h_dim2) * h_dim1];
#line 868 "MB03WD.f"
		    dlarfg_(&c__2, &h__[i__ + (i__ + j * h_dim2) * h_dim1], v,
			     &c__1, &tau);
#line 869 "MB03WD.f"
		    h__[i__ + (i__ - 1 + j * h_dim2) * h_dim1] = 0.;
#line 870 "MB03WD.f"
		    v[1] = 1.;

/*                 Apply G from the right to transform the columns of the */
/*                 matrix H_j in rows I1 to I-1. */

#line 875 "MB03WD.f"
		    i__2 = i__ - i1;
#line 875 "MB03WD.f"
		    dlarfx_("Right", &i__2, &c__2, v, &tau, &h__[i1 + (i__ - 
			    1 + j * h_dim2) * h_dim1], ldh1, &dwork[1], (
			    ftnlen)5);

/*                 Apply G from the left to transform the rows of the */
/*                 matrix H_(j+1) in columns I-1 to I2. */

#line 881 "MB03WD.f"
		    i__2 = i2 - i__ + 2;
#line 881 "MB03WD.f"
		    dlarfx_("Left", &c__2, &i__2, v, &tau, &h__[i__ - 1 + (
			    i__ - 1 + (j + 1) * h_dim2) * h_dim1], ldh1, &
			    dwork[1], (ftnlen)4);

#line 884 "MB03WD.f"
		    if (wantz) {

/*                    Accumulate transformations in the matrix Z_(j+1). */

#line 888 "MB03WD.f"
			dlarfx_("Right", &nz, &c__2, v, &tau, &z__[*iloz + (
				i__ - 1 + (j + 1) * z_dim2) * z_dim1], ldz1, &
				dwork[1], (ftnlen)5);
#line 890 "MB03WD.f"
		    }
#line 891 "MB03WD.f"
/* L200: */
#line 891 "MB03WD.f"
		}

#line 893 "MB03WD.f"
		h__[i__ + (i__ - 1 + jmin * h_dim2) * h_dim1] = 0.;

#line 895 "MB03WD.f"
	    } else {
#line 896 "MB03WD.f"
		if (jmax > 0 && wi[i__ - 1] == 0.) {
#line 896 "MB03WD.f"
		    dlartg_(&h__[i__ - 1 + (i__ - 1 + h_dim2) * h_dim1], &h__[
			    i__ + (i__ - 1 + h_dim2) * h_dim1], &cs, &sn, &
			    tau);
#line 896 "MB03WD.f"
		}

/*              Apply the transformation to H. */

#line 902 "MB03WD.f"
		i__1 = i2 - i__ + 2;
#line 902 "MB03WD.f"
		drot_(&i__1, &h__[i__ - 1 + (i__ - 1 + h_dim2) * h_dim1], 
			ldh1, &h__[i__ + (i__ - 1 + h_dim2) * h_dim1], ldh1, &
			cs, &sn);
#line 904 "MB03WD.f"
		i__1 = i__ - i1 + 1;
#line 904 "MB03WD.f"
		drot_(&i__1, &h__[i1 + (i__ - 1 + *p * h_dim2) * h_dim1], &
			c__1, &h__[i1 + (i__ + *p * h_dim2) * h_dim1], &c__1, 
			&cs, &sn);
#line 906 "MB03WD.f"
		if (wantz) {

/*                 Apply transformation to Z_1. */

#line 910 "MB03WD.f"
		    drot_(&nz, &z__[*iloz + (i__ - 1 + z_dim2) * z_dim1], &
			    c__1, &z__[*iloz + (i__ + z_dim2) * z_dim1], &
			    c__1, &cs, &sn);
#line 912 "MB03WD.f"
		}

/* Computing MAX */
#line 914 "MB03WD.f"
		i__2 = 2, i__3 = jmax + 1;
#line 914 "MB03WD.f"
		i__1 = max(i__2,i__3);
#line 914 "MB03WD.f"
		for (j = *p; j >= i__1; --j) {

/*                 Compute G1 to annihilate from the left the (I,I-1) */
/*                 element of the matrix H_j. */

#line 919 "MB03WD.f"
		    v[0] = h__[i__ + (i__ - 1 + j * h_dim2) * h_dim1];
#line 920 "MB03WD.f"
		    dlarfg_(&c__2, &h__[i__ - 1 + (i__ - 1 + j * h_dim2) * 
			    h_dim1], v, &c__1, &tau);
#line 921 "MB03WD.f"
		    h__[i__ + (i__ - 1 + j * h_dim2) * h_dim1] = 0.;

/*                 Apply G1 from the left to transform the rows of the */
/*                 matrix H_j in columns I to I2. */

#line 926 "MB03WD.f"
		    i__2 = i2 - i__ + 1;
#line 926 "MB03WD.f"
		    mb04py_("Left", &c__2, &i__2, v, &tau, &h__[i__ - 1 + (
			    i__ + j * h_dim2) * h_dim1], ldh1, &dwork[1], (
			    ftnlen)4);

/*                 Apply G1 from the right to transform the columns of */
/*                 the matrix H_(j-1) in rows I1 to I. */

#line 932 "MB03WD.f"
		    i__2 = i__ - i1 + 1;
#line 932 "MB03WD.f"
		    mb04py_("Right", &i__2, &c__2, v, &tau, &h__[i1 + (i__ - 
			    1 + (j - 1) * h_dim2) * h_dim1], ldh1, &dwork[1], 
			    (ftnlen)5);

#line 935 "MB03WD.f"
		    if (wantz) {

/*                    Apply G1 to Z_j. */

#line 939 "MB03WD.f"
			mb04py_("Right", &nz, &c__2, v, &tau, &z__[*iloz + (
				i__ - 1 + j * z_dim2) * z_dim1], ldz1, &dwork[
				1], (ftnlen)5);
#line 941 "MB03WD.f"
		    }
#line 942 "MB03WD.f"
/* L210: */
#line 942 "MB03WD.f"
		}

#line 944 "MB03WD.f"
		if (jmax > 0) {
#line 945 "MB03WD.f"
		    h__[i__ + (i__ - 1 + h_dim2) * h_dim1] = 0.;
#line 946 "MB03WD.f"
		    h__[i__ + (i__ - 1 + jmax * h_dim2) * h_dim1] = 0.;
#line 947 "MB03WD.f"
		} else {
#line 948 "MB03WD.f"
		    if (hh21 == 0.) {
#line 948 "MB03WD.f"
			h__[i__ + (i__ - 1 + h_dim2) * h_dim1] = 0.;
#line 948 "MB03WD.f"
		    }
#line 950 "MB03WD.f"
		}
#line 951 "MB03WD.f"
	    }
#line 952 "MB03WD.f"
	}
#line 953 "MB03WD.f"
    }

/*     Decrement number of remaining iterations, and return to start of */
/*     the main loop with new value of I. */

#line 958 "MB03WD.f"
    itn -= its;
#line 959 "MB03WD.f"
    i__ = l - 1;
#line 960 "MB03WD.f"
    if (i__ >= *ilo) {
#line 960 "MB03WD.f"
	goto L40;
#line 960 "MB03WD.f"
    }

#line 963 "MB03WD.f"
    return 0;

/* *** Last line of MB03WD *** */
} /* mb03wd_ */

