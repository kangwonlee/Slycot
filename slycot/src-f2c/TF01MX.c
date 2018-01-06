#line 1 "TF01MX.f"
/* TF01MX.f -- translated by f2c (version 20100827).
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

#line 1 "TF01MX.f"
/* Table of constant values */

static doublereal c_b4 = 0.;
static doublereal c_b8 = 1.;
static integer c__1 = 1;
static integer c_n1 = -1;

/* Subroutine */ int tf01mx_(integer *n, integer *m, integer *p, integer *ny, 
	doublereal *s, integer *lds, doublereal *u, integer *ldu, doublereal *
	x, doublereal *y, integer *ldy, doublereal *dwork, integer *ldwork, 
	integer *info)
{
    /* System generated locals */
    integer s_dim1, s_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3, i__4;

    /* Local variables */
    static integer i__, j, k, ic, nb, nf, nm, iu, np, iw, jw, iy, ns, n2m, 
	    n2p;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlaset_(char *, integer *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, ftnlen), 
	    xerbla_(char *, integer *, ftnlen);
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);


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

/*     To compute the output sequence of a linear time-invariant */
/*     open-loop system given by its discrete-time state-space model */
/*     with an (N+P)-by-(N+M) general system matrix S, */

/*            ( A  B ) */
/*        S = (      ) . */
/*            ( C  D ) */

/*     The initial state vector x(1) must be supplied by the user. */

/*     The input and output trajectories are stored as in the SLICOT */
/*     Library routine TF01MY. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     N       (input) INTEGER */
/*             The order of the matrix A.  N >= 0. */

/*     M       (input) INTEGER */
/*             The number of system inputs.  M >= 0. */

/*     P       (input) INTEGER */
/*             The number of system outputs.  P >= 0. */

/*     NY      (input) INTEGER */
/*             The number of output vectors y(k) to be computed. */
/*             NY >= 0. */

/*     S       (input) DOUBLE PRECISION array, dimension (LDS,N+M) */
/*             The leading (N+P)-by-(N+M) part of this array must contain */
/*             the system matrix S. */

/*     LDS     INTEGER */
/*             The leading dimension of array S.  LDS >= MAX(1,N+P). */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU,M) */
/*             The leading NY-by-M part of this array must contain the */
/*             input vector sequence u(k), for k = 1,2,...,NY. */
/*             Specifically, the k-th row of U must contain u(k)'. */

/*     LDU     INTEGER */
/*             The leading dimension of array U.  LDU >= MAX(1,NY). */

/*     X       (input/output) DOUBLE PRECISION array, dimension (N) */
/*             On entry, this array must contain the initial state vector */
/*             x(1) which consists of the N initial states of the system. */
/*             On exit, this array contains the final state vector */
/*             x(NY+1) of the N states of the system at instant NY+1. */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY,P) */
/*             The leading NY-by-P part of this array contains the output */
/*             vector sequence y(1),y(2),...,y(NY) such that the k-th */
/*             row of Y contains y(k)' (the outputs at instant k), */
/*             for k = 1,2,...,NY. */

/*     LDY     INTEGER */
/*             The leading dimension of array Y.  LDY >= MAX(1,NY). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= 0,        if MIN(N,P,NY) = 0;  otherwise, */
/*             LDWORK >= N+P,      if M = 0; */
/*             LDWORK >= 2*N+M+P,  if M > 0. */
/*             For better performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     Given an initial state vector x(1), the output vector sequence */
/*     y(1), y(2),..., y(NY) is obtained via the formulae */

/*        ( x(k+1) )     ( x(k) ) */
/*        (        ) = S (      ) , */
/*        (  y(k)  )     ( u(k) ) */

/*     where each element y(k) is a vector of length P containing the */
/*     outputs at instant k, and k = 1,2,...,NY. */

/*     REFERENCES */

/*     [1] Luenberger, D.G. */
/*         Introduction to Dynamic Systems: Theory, Models and */
/*         Applications. */
/*         John Wiley & Sons, New York, 1979. */

/*     NUMERICAL ASPECTS */

/*     The algorithm requires approximately (N + M) x (N + P) x NY */
/*     multiplications and additions. */

/*     FURTHER COMMENTS */

/*     The implementation exploits data locality as much as possible, */
/*     given the workspace length. */

/*     CONTRIBUTOR */

/*     V. Sima, Katholieke Univ. Leuven, Belgium, Feb. 2002. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Discrete-time system, multivariable system, state-space model, */
/*     state-space representation, time response. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 165 "TF01MX.f"
    /* Parameter adjustments */
#line 165 "TF01MX.f"
    s_dim1 = *lds;
#line 165 "TF01MX.f"
    s_offset = 1 + s_dim1;
#line 165 "TF01MX.f"
    s -= s_offset;
#line 165 "TF01MX.f"
    u_dim1 = *ldu;
#line 165 "TF01MX.f"
    u_offset = 1 + u_dim1;
#line 165 "TF01MX.f"
    u -= u_offset;
#line 165 "TF01MX.f"
    --x;
#line 165 "TF01MX.f"
    y_dim1 = *ldy;
#line 165 "TF01MX.f"
    y_offset = 1 + y_dim1;
#line 165 "TF01MX.f"
    y -= y_offset;
#line 165 "TF01MX.f"
    --dwork;
#line 165 "TF01MX.f"

#line 165 "TF01MX.f"
    /* Function Body */
#line 165 "TF01MX.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 169 "TF01MX.f"
    np = *n + *p;
#line 170 "TF01MX.f"
    nm = *n + *m;
#line 171 "TF01MX.f"
    iw = nm + np;
#line 172 "TF01MX.f"
    if (*n < 0) {
#line 173 "TF01MX.f"
	*info = -1;
#line 174 "TF01MX.f"
    } else if (*m < 0) {
#line 175 "TF01MX.f"
	*info = -2;
#line 176 "TF01MX.f"
    } else if (*p < 0) {
#line 177 "TF01MX.f"
	*info = -3;
#line 178 "TF01MX.f"
    } else if (*ny < 0) {
#line 179 "TF01MX.f"
	*info = -4;
#line 180 "TF01MX.f"
    } else if (*lds < max(1,np)) {
#line 181 "TF01MX.f"
	*info = -6;
#line 182 "TF01MX.f"
    } else if (*ldu < max(1,*ny)) {
#line 183 "TF01MX.f"
	*info = -8;
#line 184 "TF01MX.f"
    } else if (*ldy < max(1,*ny)) {
#line 185 "TF01MX.f"
	*info = -11;
#line 186 "TF01MX.f"
    } else {
/* Computing MIN */
#line 187 "TF01MX.f"
	i__1 = min(*n,*p);
#line 187 "TF01MX.f"
	if (min(i__1,*ny) == 0) {
#line 188 "TF01MX.f"
	    jw = 0;
#line 189 "TF01MX.f"
	} else if (*m == 0) {
#line 190 "TF01MX.f"
	    jw = np;
#line 191 "TF01MX.f"
	} else {
#line 192 "TF01MX.f"
	    jw = iw;
#line 193 "TF01MX.f"
	}
#line 194 "TF01MX.f"
	if (*ldwork < jw) {
#line 194 "TF01MX.f"
	    *info = -13;
#line 194 "TF01MX.f"
	}
#line 196 "TF01MX.f"
    }

#line 198 "TF01MX.f"
    if (*info != 0) {

/*        Error return. */

#line 202 "TF01MX.f"
	i__1 = -(*info);
#line 202 "TF01MX.f"
	xerbla_("TF01MX", &i__1, (ftnlen)6);
#line 203 "TF01MX.f"
	return 0;
#line 204 "TF01MX.f"
    }

/*     Quick return if possible. */

#line 208 "TF01MX.f"
    if (min(*ny,*p) == 0) {
#line 209 "TF01MX.f"
	return 0;
#line 210 "TF01MX.f"
    } else if (*n == 0) {

/*        Non-dynamic system: compute the output vectors. */

#line 214 "TF01MX.f"
	if (*m == 0) {
#line 215 "TF01MX.f"
	    dlaset_("Full", ny, p, &c_b4, &c_b4, &y[y_offset], ldy, (ftnlen)4)
		    ;
#line 216 "TF01MX.f"
	} else {
#line 217 "TF01MX.f"
	    dgemm_("No transpose", "Transpose", ny, p, m, &c_b8, &u[u_offset],
		     ldu, &s[s_offset], lds, &c_b4, &y[y_offset], ldy, (
		    ftnlen)12, (ftnlen)9);
#line 219 "TF01MX.f"
	}
#line 220 "TF01MX.f"
	return 0;
#line 221 "TF01MX.f"
    }

/*     Determine the block size (taken as for LAPACK routine DGETRF). */

#line 225 "TF01MX.f"
    i__1 = max(*m,*p);
#line 225 "TF01MX.f"
    nb = ilaenv_(&c__1, "DGETRF", " ", ny, &i__1, &c_n1, &c_n1, (ftnlen)6, (
	    ftnlen)1);

/*     Find the number of state vectors, extended with inputs (if M > 0) */
/*     and outputs, that can be accommodated in the provided workspace. */

/* Computing MIN */
#line 230 "TF01MX.f"
    i__1 = *ldwork / jw, i__2 = nb * nb / jw, i__1 = min(i__1,i__2);
#line 230 "TF01MX.f"
    ns = min(i__1,*ny);
#line 231 "TF01MX.f"
    n2p = *n + np;

#line 233 "TF01MX.f"
    if (*m == 0) {

/*        System with no inputs. */
/*        Workspace: need   N + P; */
/*                   prefer larger. */

#line 239 "TF01MX.f"
	if (ns <= 1 || *ny * *p <= nb * nb) {
#line 240 "TF01MX.f"
	    iy = *n + 1;

/*           LDWORK < 2*(N+P), or small problem. */
/*           One row of array Y is computed for each loop index value. */

#line 245 "TF01MX.f"
	    i__1 = *ny;
#line 245 "TF01MX.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Compute */

/*              /x(i+1)\    /A\ */
/*              |      | =  | | * x(i). */
/*              \ y(i) /    \C/ */

#line 253 "TF01MX.f"
		dgemv_("NoTranspose", &np, n, &c_b8, &s[s_offset], lds, &x[1],
			 &c__1, &c_b4, &dwork[1], &c__1, (ftnlen)11);
#line 255 "TF01MX.f"
		dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);
#line 256 "TF01MX.f"
		dcopy_(p, &dwork[iy], &c__1, &y[i__ + y_dim1], ldy);
#line 257 "TF01MX.f"
/* L10: */
#line 257 "TF01MX.f"
	    }

#line 259 "TF01MX.f"
	} else {

/*           LDWORK >= 2*(N+P), and large problem. */
/*           NS rows of array Y are computed before being saved. */

#line 264 "TF01MX.f"
	    nf = *ny / ns * ns;
#line 265 "TF01MX.f"
	    dcopy_(n, &x[1], &c__1, &dwork[1], &c__1);

#line 267 "TF01MX.f"
	    i__1 = nf;
#line 267 "TF01MX.f"
	    i__2 = ns;
#line 267 "TF01MX.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {

/*              Compute the current NS extended state vectors in the */
/*              workspace: */

/*              /x(i+1)\    /A\ */
/*              |      | =  | | * x(i),  i = 1 : ns - 1. */
/*              \ y(i) /    \C/ */

#line 276 "TF01MX.f"
		i__3 = (ns - 1) * np;
#line 276 "TF01MX.f"
		i__4 = np;
#line 276 "TF01MX.f"
		for (ic = 1; i__4 < 0 ? ic >= i__3 : ic <= i__3; ic += i__4) {
#line 277 "TF01MX.f"
		    dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			    dwork[ic], &c__1, &c_b4, &dwork[ic + np], &c__1, (
			    ftnlen)12);
#line 279 "TF01MX.f"
/* L20: */
#line 279 "TF01MX.f"
		}

/*              Prepare the next iteration. */

#line 283 "TF01MX.f"
		dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			dwork[(ns - 1) * np + 1], &c__1, &c_b4, &dwork[1], &
			c__1, (ftnlen)12);

/*              Transpose the NS output vectors in the corresponding part */
/*              of Y (column-wise). */

#line 289 "TF01MX.f"
		i__4 = *p;
#line 289 "TF01MX.f"
		for (j = 1; j <= i__4; ++j) {
#line 290 "TF01MX.f"
		    i__3 = ns - 1;
#line 290 "TF01MX.f"
		    dcopy_(&i__3, &dwork[n2p + j], &np, &y[i__ + j * y_dim1], 
			    &c__1);
#line 291 "TF01MX.f"
		    y[i__ + ns - 1 + j * y_dim1] = dwork[*n + j];
#line 292 "TF01MX.f"
/* L30: */
#line 292 "TF01MX.f"
		}

#line 294 "TF01MX.f"
/* L40: */
#line 294 "TF01MX.f"
	    }

#line 296 "TF01MX.f"
	    ns = *ny - nf;

#line 298 "TF01MX.f"
	    if (ns > 1) {

/*              Compute similarly the last NS output vectors. */

#line 302 "TF01MX.f"
		i__2 = (ns - 1) * np;
#line 302 "TF01MX.f"
		i__1 = np;
#line 302 "TF01MX.f"
		for (ic = 1; i__1 < 0 ? ic >= i__2 : ic <= i__2; ic += i__1) {
#line 303 "TF01MX.f"
		    dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			    dwork[ic], &c__1, &c_b4, &dwork[ic + np], &c__1, (
			    ftnlen)12);
#line 305 "TF01MX.f"
/* L50: */
#line 305 "TF01MX.f"
		}

#line 307 "TF01MX.f"
		dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			dwork[(ns - 1) * np + 1], &c__1, &c_b4, &dwork[1], &
			c__1, (ftnlen)12);

#line 310 "TF01MX.f"
		i__1 = *p;
#line 310 "TF01MX.f"
		for (j = 1; j <= i__1; ++j) {
#line 311 "TF01MX.f"
		    i__2 = ns - 1;
#line 311 "TF01MX.f"
		    dcopy_(&i__2, &dwork[n2p + j], &np, &y[nf + 1 + j * 
			    y_dim1], &c__1);
#line 312 "TF01MX.f"
		    y[nf + ns + j * y_dim1] = dwork[*n + j];
#line 313 "TF01MX.f"
/* L60: */
#line 313 "TF01MX.f"
		}

#line 315 "TF01MX.f"
	    } else if (ns == 1) {

/*              Compute similarly the last NS = 1 output vectors. */

#line 319 "TF01MX.f"
		dcopy_(n, &dwork[1], &c__1, &dwork[np + 1], &c__1);
#line 320 "TF01MX.f"
		dgemv_("No transpose", &np, n, &c_b8, &s[s_offset], lds, &
			dwork[np + 1], &c__1, &c_b4, &dwork[1], &c__1, (
			ftnlen)12);
#line 322 "TF01MX.f"
		dcopy_(p, &dwork[*n + 1], &c__1, &y[nf + 1 + y_dim1], ldy);

#line 324 "TF01MX.f"
	    }

/*           Set the final state vector. */

#line 328 "TF01MX.f"
	    dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);

#line 330 "TF01MX.f"
	}

#line 332 "TF01MX.f"
    } else {

/*        General case. */
/*        Workspace: need   2*N + M + P; */
/*                   prefer larger. */

#line 338 "TF01MX.f"
	dcopy_(n, &x[1], &c__1, &dwork[1], &c__1);

#line 340 "TF01MX.f"
	if (ns <= 1 || *ny * (*m + *p) <= nb * nb) {
#line 341 "TF01MX.f"
	    iu = *n + 1;
#line 342 "TF01MX.f"
	    jw = iu + *m;
#line 343 "TF01MX.f"
	    iy = jw + *n;

/*           LDWORK < 2*(2*N+M+P), or small problem. */
/*           One row of array Y is computed for each loop index value. */

#line 348 "TF01MX.f"
	    i__1 = *ny;
#line 348 "TF01MX.f"
	    for (i__ = 1; i__ <= i__1; ++i__) {

/*              Compute */

/*              /x(i+1)\    /A, B\   /x(i)\ */
/*              |      | =  |    | * |    | . */
/*              \ y(i) /    \C, D/   \u(i)/ */

#line 356 "TF01MX.f"
		dcopy_(m, &u[i__ + u_dim1], ldu, &dwork[iu], &c__1);
#line 357 "TF01MX.f"
		dgemv_("NoTranspose", &np, &nm, &c_b8, &s[s_offset], lds, &
			dwork[1], &c__1, &c_b4, &dwork[jw], &c__1, (ftnlen)11)
			;
#line 359 "TF01MX.f"
		dcopy_(n, &dwork[jw], &c__1, &dwork[1], &c__1);
#line 360 "TF01MX.f"
		dcopy_(p, &dwork[iy], &c__1, &y[i__ + y_dim1], ldy);
#line 361 "TF01MX.f"
/* L70: */
#line 361 "TF01MX.f"
	    }

#line 363 "TF01MX.f"
	} else {

/*           LDWORK >= 2*(2*N+M+P), and large problem. */
/*           NS rows of array Y are computed before being saved. */

#line 368 "TF01MX.f"
	    nf = *ny / ns * ns;
#line 369 "TF01MX.f"
	    n2m = *n + nm;

#line 371 "TF01MX.f"
	    i__1 = nf;
#line 371 "TF01MX.f"
	    i__2 = ns;
#line 371 "TF01MX.f"
	    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
#line 372 "TF01MX.f"
		jw = 1;

/*              Compute the current NS extended state vectors in the */
/*              workspace: */

/*              /x(i+1)\    /A, B\   /x(i)\ */
/*              |      | =  |    | * |    | ,  i = 1 : ns - 1. */
/*              \ y(i) /    \C, D/   \u(i)/ */

#line 381 "TF01MX.f"
		i__4 = *m;
#line 381 "TF01MX.f"
		for (j = 1; j <= i__4; ++j) {
#line 382 "TF01MX.f"
		    dcopy_(&ns, &u[i__ + j * u_dim1], &c__1, &dwork[*n + j], &
			    iw);
#line 383 "TF01MX.f"
/* L80: */
#line 383 "TF01MX.f"
		}

#line 385 "TF01MX.f"
		i__4 = ns - 1;
#line 385 "TF01MX.f"
		for (k = 1; k <= i__4; ++k) {
#line 386 "TF01MX.f"
		    dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds,
			     &dwork[jw], &c__1, &c_b4, &dwork[jw + nm], &c__1,
			     (ftnlen)12);
#line 388 "TF01MX.f"
		    jw += nm;
#line 389 "TF01MX.f"
		    dcopy_(n, &dwork[jw], &c__1, &dwork[jw + np], &c__1);
#line 390 "TF01MX.f"
		    jw += np;
#line 391 "TF01MX.f"
/* L90: */
#line 391 "TF01MX.f"
		}

/*              Prepare the next iteration. */

#line 395 "TF01MX.f"
		dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds, &
			dwork[jw], &c__1, &c_b4, &dwork[jw + nm], &c__1, (
			ftnlen)12);
#line 397 "TF01MX.f"
		dcopy_(n, &dwork[jw + nm], &c__1, &dwork[1], &c__1);

/*              Transpose the NS output vectors in the corresponding part */
/*              of Y (column-wise). */

#line 402 "TF01MX.f"
		i__4 = *p;
#line 402 "TF01MX.f"
		for (j = 1; j <= i__4; ++j) {
#line 403 "TF01MX.f"
		    dcopy_(&ns, &dwork[n2m + j], &iw, &y[i__ + j * y_dim1], &
			    c__1);
#line 404 "TF01MX.f"
/* L100: */
#line 404 "TF01MX.f"
		}

#line 406 "TF01MX.f"
/* L110: */
#line 406 "TF01MX.f"
	    }

#line 408 "TF01MX.f"
	    ns = *ny - nf;

#line 410 "TF01MX.f"
	    if (ns > 1) {
#line 411 "TF01MX.f"
		jw = 1;

/*              Compute similarly the last NS output vectors. */

#line 415 "TF01MX.f"
		i__2 = *m;
#line 415 "TF01MX.f"
		for (j = 1; j <= i__2; ++j) {
#line 416 "TF01MX.f"
		    dcopy_(&ns, &u[nf + 1 + j * u_dim1], &c__1, &dwork[*n + j]
			    , &iw);
#line 417 "TF01MX.f"
/* L120: */
#line 417 "TF01MX.f"
		}

#line 419 "TF01MX.f"
		i__2 = ns - 1;
#line 419 "TF01MX.f"
		for (k = 1; k <= i__2; ++k) {
#line 420 "TF01MX.f"
		    dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds,
			     &dwork[jw], &c__1, &c_b4, &dwork[jw + nm], &c__1,
			     (ftnlen)12);
#line 422 "TF01MX.f"
		    jw += nm;
#line 423 "TF01MX.f"
		    dcopy_(n, &dwork[jw], &c__1, &dwork[jw + np], &c__1);
#line 424 "TF01MX.f"
		    jw += np;
#line 425 "TF01MX.f"
/* L130: */
#line 425 "TF01MX.f"
		}

#line 427 "TF01MX.f"
		dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds, &
			dwork[jw], &c__1, &c_b4, &dwork[jw + nm], &c__1, (
			ftnlen)12);
#line 429 "TF01MX.f"
		dcopy_(n, &dwork[jw + nm], &c__1, &dwork[1], &c__1);

#line 431 "TF01MX.f"
		i__2 = *p;
#line 431 "TF01MX.f"
		for (j = 1; j <= i__2; ++j) {
#line 432 "TF01MX.f"
		    dcopy_(&ns, &dwork[n2m + j], &iw, &y[nf + 1 + j * y_dim1],
			     &c__1);
#line 433 "TF01MX.f"
/* L140: */
#line 433 "TF01MX.f"
		}

#line 435 "TF01MX.f"
	    } else if (ns == 1) {

/*              Compute similarly the last NS = 1 output vectors. */

#line 439 "TF01MX.f"
		dcopy_(n, &dwork[1], &c__1, &dwork[np + 1], &c__1);
#line 440 "TF01MX.f"
		dcopy_(m, &u[nf + 1 + u_dim1], ldu, &dwork[n2p + 1], &c__1);
#line 441 "TF01MX.f"
		dgemv_("No transpose", &np, &nm, &c_b8, &s[s_offset], lds, &
			dwork[np + 1], &c__1, &c_b4, &dwork[1], &c__1, (
			ftnlen)12);
#line 443 "TF01MX.f"
		dcopy_(p, &dwork[*n + 1], &c__1, &y[nf + 1 + y_dim1], ldy);

#line 445 "TF01MX.f"
	    }

#line 447 "TF01MX.f"
	}

/*        Set the final state vector. */

#line 451 "TF01MX.f"
	dcopy_(n, &dwork[1], &c__1, &x[1], &c__1);

#line 453 "TF01MX.f"
    }

#line 455 "TF01MX.f"
    return 0;
/* *** Last line of TF01MX *** */
} /* tf01mx_ */

