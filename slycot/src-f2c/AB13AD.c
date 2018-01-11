#line 1 "AB13AD.f"
/* AB13AD.f -- translated by f2c (version 20100827).
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

#line 1 "AB13AD.f"
doublereal ab13ad_(char *dico, char *equil, integer *n, integer *m, integer *
	p, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, doublereal *c__, integer *ldc, integer *ns, doublereal *
	hsv, doublereal *dwork, integer *ldwork, integer *info, ftnlen 
	dico_len, ftnlen equil_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, i__1, i__2, 
	    i__3;
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer kt, kw, kw1, kw2, ierr;
    extern doublereal ab13ax_(char *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen);
    extern /* Subroutine */ int tb01id_(char *, integer *, integer *, integer 
	    *, doublereal *, doublereal *, integer *, doublereal *, integer *,
	     doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    tb01kd_(char *, char *, char *, integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr;
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal maxred, alpwrk, wrkopt;


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

/*     To compute the Hankel-norm of the ALPHA-stable projection of the */
/*     transfer-function matrix G of the state-space system (A,B,C). */

/*     FUNCTION VALUE */

/*     AB13AD  DOUBLE PRECISION */
/*             The Hankel-norm of the ALPHA-stable projection of G */
/*             (if INFO = 0). */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     EQUIL   CHARACTER*1 */
/*             Specifies whether the user wishes to preliminarily */
/*             equilibrate the triplet (A,B,C) as follows: */
/*             = 'S':  perform equilibration (scaling); */
/*             = 'N':  do not perform equilibration. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the state-space representation, i.e. the */
/*             order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     ALPHA   (input) DOUBLE PRECISION */
/*             Specifies the ALPHA-stability boundary for the eigenvalues */
/*             of the state dynamics matrix A. For a continuous-time */
/*             system (DICO = 'C'), ALPHA <= 0 is the boundary value for */
/*             the real parts of eigenvalues, while for a discrete-time */
/*             system (DICO = 'D'), 0 <= ALPHA <= 1 represents the */
/*             boundary value for the moduli of eigenvalues. */
/*             The ALPHA-stability domain does not include the boundary */
/*             (see the Note below). */

/*     A       (input/output) DOUBLE PRECISION array, dimension (LDA,N) */
/*             On entry, the leading N-by-N part of this array must */
/*             contain the state dynamics matrix A. */
/*             On exit, if INFO = 0, the leading N-by-N part of this */
/*             array contains the state dynamics matrix A in a block */
/*             diagonal real Schur form with its eigenvalues reordered */
/*             and separated. The resulting A has two diagonal blocks. */
/*             The leading NS-by-NS part of A has eigenvalues in the */
/*             ALPHA-stability domain and the trailing (N-NS) x (N-NS) */
/*             part has eigenvalues outside the ALPHA-stability domain. */
/*             Note: The ALPHA-stability domain is defined either */
/*                   as the open half complex plane left to ALPHA, */
/*                   for a continous-time system (DICO = 'C'), or the */
/*                   interior of the ALPHA-radius circle centered in the */
/*                   origin, for a discrete-time system (DICO = 'D'). */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input/output) DOUBLE PRECISION array, dimension (LDB,M) */
/*             On entry, the leading N-by-M part of this array must */
/*             contain the original input/state matrix B. */
/*             On exit, if INFO = 0, the leading N-by-M part of this */
/*             array contains the input/state matrix B of the transformed */
/*             system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input/output) DOUBLE PRECISION array, dimension (LDC,N) */
/*             On entry, the leading P-by-N part of this array must */
/*             contain the original state/output matrix C. */
/*             On exit, if INFO = 0, the leading P-by-N part of this */
/*             array contains the state/output matrix C of the */
/*             transformed system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     NS      (output) INTEGER */
/*             The dimension of the ALPHA-stable subsystem. */

/*     HSV     (output) DOUBLE PRECISION array, dimension (N) */
/*             If INFO = 0, the leading NS elements of HSV contain the */
/*             Hankel singular values of the ALPHA-stable part of the */
/*             original system ordered decreasingly. */
/*             HSV(1) is the Hankel norm of the ALPHA-stable subsystem. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1,N*(MAX(N,M,P)+5)+N*(N+1)/2). */
/*             For optimum performance LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  the computation of the ordered real Schur form of A */
/*                   failed; */
/*             = 2:  the separation of the ALPHA-stable/unstable diagonal */
/*                   blocks failed because of very close eigenvalues; */
/*             = 3:  the computed ALPHA-stable part is just stable, */
/*                   having stable eigenvalues very near to the imaginary */
/*                   axis (if DICO = 'C') or to the unit circle */
/*                   (if DICO = 'D'); */
/*             = 4:  the computation of Hankel singular values failed. */

/*     METHOD */

/*     Let be the following linear system */

/*          d[x(t)] = Ax(t) + Bu(t) */
/*          y(t)    = Cx(t)                               (1) */

/*     where d[x(t)] is dx(t)/dt for a continuous-time system and x(t+1) */
/*     for a discrete-time system, and let G be the corresponding */
/*     transfer-function matrix. The following procedure is used to */
/*     compute the Hankel-norm of the ALPHA-stable projection of G: */

/*     1) Decompose additively G as */

/*          G = G1 + G2 */

/*        such that G1 = (As,Bs,Cs) has only ALPHA-stable poles and */
/*        G2 = (Au,Bu,Cu) has only ALPHA-unstable poles. */
/*        For the computation of the additive decomposition, the */
/*        algorithm presented in [1] is used. */

/*     2) Compute the Hankel-norm of ALPHA-stable projection G1 as the */
/*        the maximum Hankel singular value of the system (As,Bs,Cs). */
/*        The computation of the Hankel singular values is performed */
/*        by using the square-root method of [2]. */

/*     REFERENCES */

/*     [1] Safonov, M.G., Jonckheere, E.A., Verma, M. and Limebeer, D.J. */
/*         Synthesis of positive real multivariable feedback systems, */
/*         Int. J. Control, Vol. 45, pp. 817-842, 1987. */

/*     [2] Tombs, M.S. and Postlethwaite, I. */
/*         Truncated balanced realization of stable, non-minimal */
/*         state-space systems. */
/*         Int. J. Control, Vol. 46, pp. 1319-1330, 1987. */

/*     NUMERICAL ASPECTS */

/*     The implemented method relies on a square-root technique. */
/*                                     3 */
/*     The algorithms require about 17N  floating point operations. */

/*     CONTRIBUTOR */

/*     C. Oara and A. Varga, German Aerospace Center, */
/*     DLR Oberpfaffenhofen, July 1998. */
/*     Based on the RASP routine SHANRM. */

/*     REVISIONS */

/*     Nov. 1998, V. Sima, Research Institute for Informatics, Bucharest. */
/*     Dec. 1998, V. Sima, Katholieke Univ. Leuven, Leuven. */
/*     Oct. 2001, V. Sima, Research Institute for Informatics, Bucharest. */

/*     KEYWORDS */

/*     Additive spectral decomposition, model reduction, */
/*     multivariable system, state-space model, system norms. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 231 "AB13AD.f"
    /* Parameter adjustments */
#line 231 "AB13AD.f"
    a_dim1 = *lda;
#line 231 "AB13AD.f"
    a_offset = 1 + a_dim1;
#line 231 "AB13AD.f"
    a -= a_offset;
#line 231 "AB13AD.f"
    b_dim1 = *ldb;
#line 231 "AB13AD.f"
    b_offset = 1 + b_dim1;
#line 231 "AB13AD.f"
    b -= b_offset;
#line 231 "AB13AD.f"
    c_dim1 = *ldc;
#line 231 "AB13AD.f"
    c_offset = 1 + c_dim1;
#line 231 "AB13AD.f"
    c__ -= c_offset;
#line 231 "AB13AD.f"
    --hsv;
#line 231 "AB13AD.f"
    --dwork;
#line 231 "AB13AD.f"

#line 231 "AB13AD.f"
    /* Function Body */
#line 231 "AB13AD.f"
    *info = 0;
#line 232 "AB13AD.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);

/*     Test the input scalar arguments. */

#line 236 "AB13AD.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 237 "AB13AD.f"
	*info = -1;
#line 238 "AB13AD.f"
    } else if (! (lsame_(equil, "S", (ftnlen)1, (ftnlen)1) || lsame_(equil, 
	    "N", (ftnlen)1, (ftnlen)1))) {
#line 240 "AB13AD.f"
	*info = -2;
#line 241 "AB13AD.f"
    } else if (*n < 0) {
#line 242 "AB13AD.f"
	*info = -3;
#line 243 "AB13AD.f"
    } else if (*m < 0) {
#line 244 "AB13AD.f"
	*info = -4;
#line 245 "AB13AD.f"
    } else if (*p < 0) {
#line 246 "AB13AD.f"
	*info = -5;
#line 247 "AB13AD.f"
    } else if (discr && (*alpha < 0. || *alpha > 1.) || ! discr && *alpha > 
	    0.) {
#line 249 "AB13AD.f"
	*info = -6;
#line 250 "AB13AD.f"
    } else if (*lda < max(1,*n)) {
#line 251 "AB13AD.f"
	*info = -8;
#line 252 "AB13AD.f"
    } else if (*ldb < max(1,*n)) {
#line 253 "AB13AD.f"
	*info = -10;
#line 254 "AB13AD.f"
    } else if (*ldc < max(1,*p)) {
#line 255 "AB13AD.f"
	*info = -12;
#line 256 "AB13AD.f"
    } else /* if(complicated condition) */ {
/* Computing MAX */
/* Computing MAX */
#line 256 "AB13AD.f"
	i__3 = max(*n,*m);
#line 256 "AB13AD.f"
	i__1 = 1, i__2 = *n * (max(i__3,*p) + 5) + *n * (*n + 1) / 2;
#line 256 "AB13AD.f"
	if (*ldwork < max(i__1,i__2)) {
#line 258 "AB13AD.f"
	    *info = -16;
#line 259 "AB13AD.f"
	}
#line 259 "AB13AD.f"
    }

#line 261 "AB13AD.f"
    if (*info != 0) {

/*        Error return. */

#line 265 "AB13AD.f"
	i__1 = -(*info);
#line 265 "AB13AD.f"
	xerbla_("AB13AD", &i__1, (ftnlen)6);
#line 266 "AB13AD.f"
	return ret_val;
#line 267 "AB13AD.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 271 "AB13AD.f"
    i__1 = min(*n,*m);
#line 271 "AB13AD.f"
    if (min(i__1,*p) == 0) {
#line 272 "AB13AD.f"
	*ns = 0;
#line 273 "AB13AD.f"
	ret_val = 0.;
#line 274 "AB13AD.f"
	dwork[1] = 1.;
#line 275 "AB13AD.f"
	return ret_val;
#line 276 "AB13AD.f"
    }

#line 278 "AB13AD.f"
    if (lsame_(equil, "S", (ftnlen)1, (ftnlen)1)) {

/*        Scale simultaneously the matrices A, B and C: */
/*        A <- inv(D)*A*D,  B <- inv(D)*B and C <- C*D, where D is a */
/*        diagonal matrix. */
/*        Workspace: N. */

#line 285 "AB13AD.f"
	maxred = 100.;
#line 286 "AB13AD.f"
	tb01id_("All", n, m, p, &maxred, &a[a_offset], lda, &b[b_offset], ldb,
		 &c__[c_offset], ldc, &dwork[1], info, (ftnlen)3);
#line 288 "AB13AD.f"
    }

/*     Correct the value of ALPHA to ensure stability. */

#line 292 "AB13AD.f"
    alpwrk = *alpha;
#line 293 "AB13AD.f"
    if (discr) {
#line 294 "AB13AD.f"
	if (*alpha == 1.) {
#line 294 "AB13AD.f"
	    alpwrk = 1. - sqrt(dlamch_("E", (ftnlen)1));
#line 294 "AB13AD.f"
	}
#line 295 "AB13AD.f"
    } else {
#line 296 "AB13AD.f"
	if (*alpha == 0.) {
#line 296 "AB13AD.f"
	    alpwrk = -sqrt(dlamch_("E", (ftnlen)1));
#line 296 "AB13AD.f"
	}
#line 297 "AB13AD.f"
    }

/*     Allocate working storage. */

#line 301 "AB13AD.f"
    kt = 1;
#line 302 "AB13AD.f"
    kw1 = *n * *n + 1;
#line 303 "AB13AD.f"
    kw2 = kw1 + *n;
#line 304 "AB13AD.f"
    kw = kw2 + *n;

/*     Reduce A to a block diagonal real Schur form, with the */
/*     ALPHA-stable part in the leading diagonal position, using a */
/*     non-orthogonal similarity transformation A <- inv(T)*A*T and */
/*     apply the transformation to B and C: B <- inv(T)*B and C <- C*T. */

/*     Workspace needed:      N*(N+2); */
/*     Additional workspace:  need   3*N; */
/*                            prefer larger. */

#line 315 "AB13AD.f"
    i__1 = *ldwork - kw + 1;
#line 315 "AB13AD.f"
    tb01kd_(dico, "Stable", "General", n, m, p, &alpwrk, &a[a_offset], lda, &
	    b[b_offset], ldb, &c__[c_offset], ldc, ns, &dwork[kt], n, &dwork[
	    kw1], &dwork[kw2], &dwork[kw], &i__1, &ierr, (ftnlen)1, (ftnlen)6,
	     (ftnlen)7);
#line 318 "AB13AD.f"
    if (ierr != 0) {
#line 319 "AB13AD.f"
	if (ierr != 3) {
#line 320 "AB13AD.f"
	    *info = 1;
#line 321 "AB13AD.f"
	} else {
#line 322 "AB13AD.f"
	    *info = 2;
#line 323 "AB13AD.f"
	}
#line 324 "AB13AD.f"
	return ret_val;
#line 325 "AB13AD.f"
    }

#line 327 "AB13AD.f"
    wrkopt = dwork[kw] + (doublereal) (kw - 1);

#line 329 "AB13AD.f"
    if (*ns == 0) {
#line 330 "AB13AD.f"
	ret_val = 0.;
#line 331 "AB13AD.f"
    } else {

/*        Workspace:  need   N*(MAX(N,M,P)+5)+N*(N+1)/2; */
/*                    prefer larger. */

#line 336 "AB13AD.f"
	ret_val = ab13ax_(dico, ns, m, p, &a[a_offset], lda, &b[b_offset], 
		ldb, &c__[c_offset], ldc, &hsv[1], &dwork[1], ldwork, &ierr, (
		ftnlen)1);

#line 339 "AB13AD.f"
	if (ierr != 0) {
#line 340 "AB13AD.f"
	    *info = ierr + 2;
#line 341 "AB13AD.f"
	    return ret_val;
#line 342 "AB13AD.f"
	}

#line 344 "AB13AD.f"
	dwork[1] = max(wrkopt,dwork[1]);
#line 345 "AB13AD.f"
    }

#line 347 "AB13AD.f"
    return ret_val;
/* *** Last line of AB13AD *** */
} /* ab13ad_ */

