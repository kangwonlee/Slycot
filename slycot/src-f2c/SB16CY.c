#line 1 "SB16CY.f"
/* SB16CY.f -- translated by f2c (version 20100827).
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

#line 1 "SB16CY.f"
/* Table of constant values */

static doublereal c_b10 = 1.;

/* Subroutine */ int sb16cy_(char *dico, char *jobcf, integer *n, integer *m, 
	integer *p, doublereal *a, integer *lda, doublereal *b, integer *ldb, 
	doublereal *c__, integer *ldc, doublereal *f, integer *ldf, 
	doublereal *g, integer *ldg, doublereal *scalec, doublereal *scaleo, 
	doublereal *s, integer *lds, doublereal *r__, integer *ldr, 
	doublereal *dwork, integer *ldwork, integer *info, ftnlen dico_len, 
	ftnlen jobcf_len)
{
    /* System generated locals */
    integer a_dim1, a_offset, b_dim1, b_offset, c_dim1, c_offset, f_dim1, 
	    f_offset, g_dim1, g_offset, r_dim1, r_offset, s_dim1, s_offset, 
	    i__1, i__2;

    /* Local variables */
    static integer me, mp, ku, kw, lw, kaw, ldu, kwi, kwr, ierr;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     sb03od_(char *, char *, char *, integer *, integer *, doublereal 
	    *, integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, ftnlen, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    static logical discr, leftw;
    extern /* Subroutine */ int dlacpy_(char *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
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

/*     To compute, for a given open-loop model (A,B,C,0), and for */
/*     given state feedback gain F and full observer gain G, */
/*     such that A+B*F and A+G*C are stable, the Cholesky factors */
/*     Su and Ru of a controllability Grammian P = Su*Su' and of */
/*     an observability Grammian Q = Ru'*Ru corresponding to a */
/*     frequency-weighted model reduction of the left or right coprime */
/*     factors of the state-feedback controller. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     DICO    CHARACTER*1 */
/*             Specifies the type of the open-loop system as follows: */
/*             = 'C':  continuous-time system; */
/*             = 'D':  discrete-time system. */

/*     JOBCF   CHARACTER*1 */
/*             Specifies whether a left or right coprime factorization */
/*             of the state-feedback controller is to be used as follows: */
/*             = 'L':  use a left coprime factorization; */
/*             = 'R':  use a right coprime factorization. */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the open-loop state-space representation, */
/*             i.e., the order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     A       (input) DOUBLE PRECISION array, dimension (LDA,N) */
/*             The leading N-by-N part of this array must contain the */
/*             state matrix A of the open-loop system. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= MAX(1,N). */

/*     B       (input) DOUBLE PRECISION array, dimension (LDB,M) */
/*             The leading N-by-M part of this array must contain the */
/*             input/state matrix B of the open-loop system. */

/*     LDB     INTEGER */
/*             The leading dimension of array B.  LDB >= MAX(1,N). */

/*     C       (input) DOUBLE PRECISION array, dimension (LDC,N) */
/*             The leading P-by-N part of this array must contain the */
/*             state/output matrix C of the open-loop system. */

/*     LDC     INTEGER */
/*             The leading dimension of array C.  LDC >= MAX(1,P). */

/*     F       (input) DOUBLE PRECISION array, dimension (LDF,N) */
/*             The leading M-by-N part of this array must contain a */
/*             stabilizing state feedback matrix. */

/*     LDF     INTEGER */
/*             The leading dimension of array F.  LDF >= MAX(1,M). */

/*     G       (input) DOUBLE PRECISION array, dimension (LDG,P) */
/*             The leading N-by-P part of this array must contain a */
/*             stabilizing observer gain matrix. */

/*     LDG     INTEGER */
/*             The leading dimension of array G.  LDG >= MAX(1,N). */

/*     SCALEC  (output) DOUBLE PRECISION */
/*             Scaling factor for the controllability Grammian. */
/*             See METHOD. */

/*     SCALEO  (output) DOUBLE PRECISION */
/*             Scaling factor for the observability Grammian. */
/*             See METHOD. */

/*     S       (output) DOUBLE PRECISION array, dimension (LDS,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the Cholesky factor Su of frequency-weighted */
/*             cotrollability Grammian P = Su*Su'. See METHOD. */

/*     LDS     INTEGER */
/*             The leading dimension of the array S.  LDS >= MAX(1,N). */

/*     R       (output) DOUBLE PRECISION array, dimension (LDR,N) */
/*             The leading N-by-N upper triangular part of this array */
/*             contains the Cholesky factor Ru of the frequency-weighted */
/*             observability Grammian Q = Ru'*Ru. See METHOD. */

/*     LDR     INTEGER */
/*             The leading dimension of the array R.  LDR >= MAX(1,N). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             On exit, if INFO = 0, DWORK(1) returns the optimal value */
/*             of LDWORK. */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= MAX(1, N*(N + MAX(N,M) + MIN(N,M) + 6)), */
/*                                                       if JOBCF = 'L'; */
/*             LDWORK >= MAX(1, N*(N + MAX(N,P) + MIN(N,P) + 6)), */
/*                                                       if JOBCF = 'R'. */
/*             For optimum performance LDWORK should be larger. */
/*             An upper bound for both cases is */
/*             LDWORK >= MAX(1, N*(N + MAX(N,M,P) + 7)). */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value; */
/*             = 1:  eigenvalue computation failure; */
/*             = 2:  the matrix A+G*C is not stable; */
/*             = 3:  the matrix A+B*F is not stable; */
/*             = 4:  the Lyapunov equation for computing the */
/*                   observability Grammian is (nearly) singular; */
/*             = 5:  the Lyapunov equation for computing the */
/*                   controllability Grammian is (nearly) singular. */

/*     METHOD */

/*     In accordance with the type of the coprime factorization */
/*     of the controller (left or right), the Cholesky factors Su and Ru */
/*     of the frequency-weighted controllability Grammian P = Su*Su' and */
/*     of the frequency-weighted observability Grammian Q = Ru'*Ru are */
/*     computed by solving appropriate Lyapunov or Stein equations [1]. */

/*     If JOBCF = 'L' and DICO = 'C', P and Q are computed as the */
/*     solutions of the following Lyapunov equations: */

/*            (A+B*F)*P + P*(A+B*F)' +  scalec^2*B*B' = 0,  (1) */

/*            (A+G*C)'*Q + Q*(A+G*C) +  scaleo^2*F'*F = 0.  (2) */

/*     If JOBCF = 'L' and DICO = 'D', P and Q are computed as the */
/*     solutions of the following Stein equations: */

/*            (A+B*F)*P*(A+B*F)' - P +  scalec^2*B*B' = 0,  (3) */

/*            (A+G*C)'*Q*(A+G*C) - Q +  scaleo^2*F'*F = 0.  (4) */

/*     If JOBCF = 'R' and DICO = 'C', P and Q are computed as the */
/*     solutions of the following Lyapunov equations: */

/*            (A+B*F)*P + P*(A+B*F)' +  scalec^2*G*G' = 0,  (5) */

/*            (A+G*C)'*Q + Q*(A+G*C) +  scaleo^2*C'*C = 0.  (6) */

/*     If JOBCF = 'R' and DICO = 'D', P and Q are computed as the */
/*     solutions of the following Stein equations: */

/*            (A+B*F)*P*(A+B*F)' - P +  scalec^2*G*G' = 0,  (7) */

/*            (A+G*C)'*Q*(A+G*C) - Q +  scaleo^2*C'*C = 0.  (8) */

/*     REFERENCES */

/*     [1] Liu, Y., Anderson, B.D.O. and Ly, O.L. */
/*         Coprime factorization controller reduction with Bezout */
/*         identity induced frequency weighting. */
/*         Automatica, vol. 26, pp. 233-249, 1990. */

/*     CONTRIBUTORS */

/*     A. Varga, German Aerospace Center, Oberpfaffenhofen, October 2000. */
/*     D. Sima, University of Bucharest, October 2000. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Oct. 2000. */

/*     REVISIONS */

/*     A. Varga, Australian National University, Canberra, November 2000. */

/*     KEYWORDS */

/*     Controller reduction, frequency weighting, multivariable system, */
/*     state-space model, state-space representation. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 233 "SB16CY.f"
    /* Parameter adjustments */
#line 233 "SB16CY.f"
    a_dim1 = *lda;
#line 233 "SB16CY.f"
    a_offset = 1 + a_dim1;
#line 233 "SB16CY.f"
    a -= a_offset;
#line 233 "SB16CY.f"
    b_dim1 = *ldb;
#line 233 "SB16CY.f"
    b_offset = 1 + b_dim1;
#line 233 "SB16CY.f"
    b -= b_offset;
#line 233 "SB16CY.f"
    c_dim1 = *ldc;
#line 233 "SB16CY.f"
    c_offset = 1 + c_dim1;
#line 233 "SB16CY.f"
    c__ -= c_offset;
#line 233 "SB16CY.f"
    f_dim1 = *ldf;
#line 233 "SB16CY.f"
    f_offset = 1 + f_dim1;
#line 233 "SB16CY.f"
    f -= f_offset;
#line 233 "SB16CY.f"
    g_dim1 = *ldg;
#line 233 "SB16CY.f"
    g_offset = 1 + g_dim1;
#line 233 "SB16CY.f"
    g -= g_offset;
#line 233 "SB16CY.f"
    s_dim1 = *lds;
#line 233 "SB16CY.f"
    s_offset = 1 + s_dim1;
#line 233 "SB16CY.f"
    s -= s_offset;
#line 233 "SB16CY.f"
    r_dim1 = *ldr;
#line 233 "SB16CY.f"
    r_offset = 1 + r_dim1;
#line 233 "SB16CY.f"
    r__ -= r_offset;
#line 233 "SB16CY.f"
    --dwork;
#line 233 "SB16CY.f"

#line 233 "SB16CY.f"
    /* Function Body */
#line 233 "SB16CY.f"
    discr = lsame_(dico, "D", (ftnlen)1, (ftnlen)1);
#line 234 "SB16CY.f"
    leftw = lsame_(jobcf, "L", (ftnlen)1, (ftnlen)1);

#line 236 "SB16CY.f"
    *info = 0;
#line 237 "SB16CY.f"
    if (leftw) {
#line 238 "SB16CY.f"
	mp = *m;
#line 239 "SB16CY.f"
    } else {
#line 240 "SB16CY.f"
	mp = *p;
#line 241 "SB16CY.f"
    }
#line 242 "SB16CY.f"
    lw = *n * (*n + max(*n,mp) + min(*n,mp) + 6);

#line 244 "SB16CY.f"
    if (! (lsame_(dico, "C", (ftnlen)1, (ftnlen)1) || discr)) {
#line 245 "SB16CY.f"
	*info = -1;
#line 246 "SB16CY.f"
    } else if (! (leftw || lsame_(jobcf, "R", (ftnlen)1, (ftnlen)1))) {
#line 247 "SB16CY.f"
	*info = -2;
#line 248 "SB16CY.f"
    } else if (*n < 0) {
#line 249 "SB16CY.f"
	*info = -3;
#line 250 "SB16CY.f"
    } else if (*m < 0) {
#line 251 "SB16CY.f"
	*info = -4;
#line 252 "SB16CY.f"
    } else if (*p < 0) {
#line 253 "SB16CY.f"
	*info = -5;
#line 254 "SB16CY.f"
    } else if (*lda < max(1,*n)) {
#line 255 "SB16CY.f"
	*info = -7;
#line 256 "SB16CY.f"
    } else if (*ldb < max(1,*n)) {
#line 257 "SB16CY.f"
	*info = -9;
#line 258 "SB16CY.f"
    } else if (*ldc < max(1,*p)) {
#line 259 "SB16CY.f"
	*info = -11;
#line 260 "SB16CY.f"
    } else if (*ldf < max(1,*m)) {
#line 261 "SB16CY.f"
	*info = -13;
#line 262 "SB16CY.f"
    } else if (*ldg < max(1,*n)) {
#line 263 "SB16CY.f"
	*info = -15;
#line 264 "SB16CY.f"
    } else if (*lds < max(1,*n)) {
#line 265 "SB16CY.f"
	*info = -19;
#line 266 "SB16CY.f"
    } else if (*ldr < max(1,*n)) {
#line 267 "SB16CY.f"
	*info = -21;
#line 268 "SB16CY.f"
    } else if (*ldwork < max(1,lw)) {
#line 269 "SB16CY.f"
	*info = -23;
#line 270 "SB16CY.f"
    }

#line 272 "SB16CY.f"
    if (*info != 0) {

/*        Error return. */

#line 276 "SB16CY.f"
	i__1 = -(*info);
#line 276 "SB16CY.f"
	xerbla_("SB16CY", &i__1, (ftnlen)6);
#line 277 "SB16CY.f"
	return 0;
#line 278 "SB16CY.f"
    }

/*     Quick return if possible. */

/* Computing MIN */
#line 282 "SB16CY.f"
    i__1 = min(*n,*m);
#line 282 "SB16CY.f"
    if (min(i__1,*p) == 0) {
#line 283 "SB16CY.f"
	*scalec = 1.;
#line 284 "SB16CY.f"
	*scaleo = 1.;
#line 285 "SB16CY.f"
	dwork[1] = 1.;
#line 286 "SB16CY.f"
	return 0;
#line 287 "SB16CY.f"
    }

/*     Allocate storage for work arrays. */

#line 291 "SB16CY.f"
    kaw = 1;
#line 292 "SB16CY.f"
    ku = kaw + *n * *n;
#line 293 "SB16CY.f"
    kwr = ku + *n * max(*n,mp);
#line 294 "SB16CY.f"
    kwi = kwr + *n;
#line 295 "SB16CY.f"
    kw = kwi + *n;

/*     Form A+G*C. */

#line 299 "SB16CY.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[kaw], n, (ftnlen)4);
#line 300 "SB16CY.f"
    dgemm_("No-transpose", "No-transpose", n, n, p, &c_b10, &g[g_offset], ldg,
	     &c__[c_offset], ldc, &c_b10, &dwork[kaw], n, (ftnlen)12, (ftnlen)
	    12);

/*     Form the factor H of the free term. */

#line 305 "SB16CY.f"
    if (leftw) {

/*        H = F. */

#line 309 "SB16CY.f"
	ldu = max(*n,*m);
#line 310 "SB16CY.f"
	me = *m;
#line 311 "SB16CY.f"
	dlacpy_("Full", m, n, &f[f_offset], ldf, &dwork[ku], &ldu, (ftnlen)4);
#line 312 "SB16CY.f"
    } else {

/*        H = C. */

#line 316 "SB16CY.f"
	ldu = max(*n,*p);
#line 317 "SB16CY.f"
	me = *p;
#line 318 "SB16CY.f"
	dlacpy_("Full", p, n, &c__[c_offset], ldc, &dwork[ku], &ldu, (ftnlen)
		4);
#line 319 "SB16CY.f"
    }

/*     Solve for the Cholesky factor Ru of Q, Q = Ru'*Ru, */
/*     the continuous-time Lyapunov equation (if DICO = 'C') */

/*        (A+G*C)'*Q + Q*(A+G*C) +  scaleo^2*H'*H = 0, */

/*     or the discrete-time Lyapunov equation (if DICO = 'D') */

/*        (A+G*C)'*Q*(A+G*C) - Q +  scaleo^2*H'*H = 0. */

/*     Workspace:  need   N*(N + MAX(N,M) + MIN(N,M) + 6) if JOBCF = 'L'; */
/*                        N*(N + MAX(N,P) + MIN(N,P) + 6) if JOBCF = 'R'. */
/*                 prefer larger. */

#line 334 "SB16CY.f"
    i__1 = *ldwork - kw + 1;
#line 334 "SB16CY.f"
    sb03od_(dico, "NoFact", "NoTransp", n, &me, &dwork[kaw], n, &r__[r_offset]
	    , ldr, &dwork[ku], &ldu, scaleo, &dwork[kwr], &dwork[kwi], &dwork[
	    kw], &i__1, &ierr, (ftnlen)1, (ftnlen)6, (ftnlen)8);
#line 337 "SB16CY.f"
    if (ierr != 0) {
#line 338 "SB16CY.f"
	if (ierr == 2) {
#line 339 "SB16CY.f"
	    *info = 2;
#line 340 "SB16CY.f"
	} else if (ierr == 1) {
#line 341 "SB16CY.f"
	    *info = 4;
#line 342 "SB16CY.f"
	} else if (ierr == 6) {
#line 343 "SB16CY.f"
	    *info = 1;
#line 344 "SB16CY.f"
	}
#line 345 "SB16CY.f"
	return 0;
#line 346 "SB16CY.f"
    }

#line 348 "SB16CY.f"
    wrkopt = (integer) dwork[kw] + kw - 1;
#line 349 "SB16CY.f"
    dlacpy_("Upper", n, n, &dwork[ku], &ldu, &r__[r_offset], ldr, (ftnlen)5);

/*     Form A+B*F. */

#line 353 "SB16CY.f"
    dlacpy_("Full", n, n, &a[a_offset], lda, &dwork[kaw], n, (ftnlen)4);
#line 354 "SB16CY.f"
    dgemm_("No-transpose", "No-transpose", n, n, m, &c_b10, &b[b_offset], ldb,
	     &f[f_offset], ldf, &c_b10, &dwork[kaw], n, (ftnlen)12, (ftnlen)
	    12);

/*     Form the factor K of the free term. */

#line 359 "SB16CY.f"
    ldu = *n;
#line 360 "SB16CY.f"
    if (leftw) {

/*        K = B. */

#line 364 "SB16CY.f"
	me = *m;
#line 365 "SB16CY.f"
	dlacpy_("Full", n, m, &b[b_offset], ldb, &dwork[ku], &ldu, (ftnlen)4);
#line 366 "SB16CY.f"
    } else {

/*        K = G. */

#line 370 "SB16CY.f"
	me = *p;
#line 371 "SB16CY.f"
	dlacpy_("Full", n, p, &g[g_offset], ldg, &dwork[ku], &ldu, (ftnlen)4);
#line 372 "SB16CY.f"
    }

/*     Solve for the Cholesky factor Su of P, P = Su*Su', */
/*     the continuous-time Lyapunov equation (if DICO = 'C') */

/*         (A+B*F)*P + P*(A+B*F)' +  scalec^2*K*K' = 0, */

/*     or the discrete-time Lyapunov equation (if DICO = 'D') */

/*         (A+B*F)*P*(A+B*F)' - P +  scalec^2*K*K' = 0. */

/*     Workspace:  need   N*(N + MAX(N,M) + MIN(N,M) + 6) if JOBCF = 'L'; */
/*                        N*(N + MAX(N,P) + MIN(N,P) + 6) if JOBCF = 'R'. */
/*                        prefer larger. */

#line 387 "SB16CY.f"
    i__1 = *ldwork - kw + 1;
#line 387 "SB16CY.f"
    sb03od_(dico, "NoFact", "Transp", n, &me, &dwork[kaw], n, &s[s_offset], 
	    lds, &dwork[ku], &ldu, scalec, &dwork[kwr], &dwork[kwi], &dwork[
	    kw], &i__1, &ierr, (ftnlen)1, (ftnlen)6, (ftnlen)6);
#line 390 "SB16CY.f"
    if (ierr != 0) {
#line 391 "SB16CY.f"
	if (ierr == 2) {
#line 392 "SB16CY.f"
	    *info = 3;
#line 393 "SB16CY.f"
	} else if (ierr == 1) {
#line 394 "SB16CY.f"
	    *info = 5;
#line 395 "SB16CY.f"
	} else if (ierr == 6) {
#line 396 "SB16CY.f"
	    *info = 1;
#line 397 "SB16CY.f"
	}
#line 398 "SB16CY.f"
	return 0;
#line 399 "SB16CY.f"
    }
/* Computing MAX */
#line 400 "SB16CY.f"
    i__1 = wrkopt, i__2 = (integer) dwork[kw] + kw - 1;
#line 400 "SB16CY.f"
    wrkopt = max(i__1,i__2);
#line 401 "SB16CY.f"
    dlacpy_("Upper", n, n, &dwork[ku], &ldu, &s[s_offset], lds, (ftnlen)5);

/*     Save the optimal workspace. */

#line 405 "SB16CY.f"
    dwork[1] = (doublereal) wrkopt;

#line 407 "SB16CY.f"
    return 0;
/* *** Last line of SB16CY *** */
} /* sb16cy_ */

