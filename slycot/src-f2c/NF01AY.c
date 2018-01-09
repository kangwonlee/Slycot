#line 1 "NF01AY.f"
/* NF01AY.f -- translated by f2c (version 20100827).
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

#line 1 "NF01AY.f"
/* Table of constant values */

static doublereal c_b10 = -2.;
static doublereal c_b11 = 0.;
static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b17 = 1.;

/* Subroutine */ int nf01ay_(integer *nsmp, integer *nz, integer *l, integer *
	ipar, integer *lipar, doublereal *wb, integer *lwb, doublereal *z__, 
	integer *ldz, doublereal *y, integer *ldy, doublereal *dwork, integer 
	*ldwork, integer *info)
{
    /* System generated locals */
    integer y_dim1, y_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, j, k, m;
    static doublereal df;
    static integer ib, mf, lj, lk, nn, nv, ws;
    static doublereal tmp;
    static integer ldwb;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    static logical last;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen),
	     dgemv_(char *, integer *, integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, doublereal *, 
	    integer *, ftnlen), dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dlabad_(doublereal *, doublereal *);
    extern doublereal dlamch_(char *, ftnlen);
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static doublereal bignum, smlnum;


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

/*     To calculate the output of a set of neural networks with the */
/*     structure */

/*             - tanh(w1'*z+b1) - */
/*           /      :             \ */
/*         z ---    :           --- sum(ws(i)*...)+ b(n+1)  --- y, */
/*           \      :             / */
/*             - tanh(wn'*z+bn) - */

/*     given the input z and the parameter vectors wi, ws, and b, */
/*     where z, w1, ..., wn are vectors of length NZ, ws is a vector */
/*     of length n, b(1), ..., b(n+1) are scalars, and n is called the */
/*     number of neurons in the hidden layer, or just number of neurons. */
/*     Such a network is used for each L output variables. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NSMP    (input) INTEGER */
/*             The number of training samples.  NSMP >= 0. */

/*     NZ      (input) INTEGER */
/*             The length of each input sample.  NZ >= 0. */

/*     L       (input) INTEGER */
/*             The length of each output sample.  L >= 0. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters needed. */
/*             IPAR(1) must contain the number of neurons, n, per output */
/*             variable, denoted NN in the sequel.  NN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of the vector IPAR.  LIPAR >= 1. */

/*     WB      (input) DOUBLE PRECISION array, dimension (LWB) */
/*             The leading (NN*(NZ+2)+1)*L part of this array must */
/*             contain the weights and biases of the network. This vector */
/*             is partitioned into L vectors of length NN*(NZ+2)+1, */
/*             WB = [ wb(1), ..., wb(L) ]. Each wb(k), k = 1, ..., L, */
/*             corresponds to one output variable, and has the structure */
/*             wb(k) = [ w1(1), ..., w1(NZ), ..., wn(1), ..., wn(NZ), */
/*                       ws(1), ..., ws(n), b(1), ..., b(n+1) ], */
/*             where wi(j) are the weights of the hidden layer, */
/*             ws(i) are the weights of the linear output layer, and */
/*             b(i) are the biases, as in the scheme above. */

/*     LWB     (input) INTEGER */
/*             The length of the array WB. */
/*             LWB >= ( NN*(NZ + 2) + 1 )*L. */

/*     Z       (input) DOUBLE PRECISION array, dimension (LDZ, NZ) */
/*             The leading NSMP-by-NZ part of this array must contain the */
/*             set of input samples, */
/*             Z = ( Z(1,1),...,Z(1,NZ); ...; Z(NSMP,1),...,Z(NSMP,NZ) ). */

/*     LDZ     INTEGER */
/*             The leading dimension of the array Z.  LDZ >= MAX(1,NSMP). */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY, L) */
/*             The leading NSMP-by-L part of this array contains the set */
/*             of output samples, */
/*             Y = ( Y(1,1),...,Y(1,L); ...; Y(NSMP,1),...,Y(NSMP,L) ). */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= MAX(1,NSMP). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK.  LDWORK >= 2*NN. */
/*             For better performance, LDWORK should be larger. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     BLAS routines are used to compute the matrix-vector products. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000, during a stay at University of Twente, NL. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Input output description, neural network, nonlinear system, */
/*     simulation, system response. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 150 "NF01AY.f"
    /* Parameter adjustments */
#line 150 "NF01AY.f"
    --ipar;
#line 150 "NF01AY.f"
    --wb;
#line 150 "NF01AY.f"
    z_dim1 = *ldz;
#line 150 "NF01AY.f"
    z_offset = 1 + z_dim1;
#line 150 "NF01AY.f"
    z__ -= z_offset;
#line 150 "NF01AY.f"
    y_dim1 = *ldy;
#line 150 "NF01AY.f"
    y_offset = 1 + y_dim1;
#line 150 "NF01AY.f"
    y -= y_offset;
#line 150 "NF01AY.f"
    --dwork;
#line 150 "NF01AY.f"

#line 150 "NF01AY.f"
    /* Function Body */
#line 150 "NF01AY.f"
    *info = 0;
#line 151 "NF01AY.f"
    nn = ipar[1];
#line 152 "NF01AY.f"
    ldwb = nn * (*nz + 2) + 1;
#line 153 "NF01AY.f"
    if (*nsmp < 0) {
#line 154 "NF01AY.f"
	*info = -1;
#line 155 "NF01AY.f"
    } else if (*nz < 0) {
#line 156 "NF01AY.f"
	*info = -2;
#line 157 "NF01AY.f"
    } else if (*l < 0) {
#line 158 "NF01AY.f"
	*info = -3;
#line 159 "NF01AY.f"
    } else if (nn < 0) {
#line 160 "NF01AY.f"
	*info = -4;
#line 161 "NF01AY.f"
    } else if (*lipar < 1) {
#line 162 "NF01AY.f"
	*info = -5;
#line 163 "NF01AY.f"
    } else if (*lwb < ldwb * *l) {
#line 164 "NF01AY.f"
	*info = -7;
#line 165 "NF01AY.f"
    } else if (*ldz < max(1,*nsmp)) {
#line 166 "NF01AY.f"
	*info = -9;
#line 167 "NF01AY.f"
    } else if (*ldy < max(1,*nsmp)) {
#line 168 "NF01AY.f"
	*info = -11;
#line 169 "NF01AY.f"
    } else if (*ldwork < nn << 1) {
#line 170 "NF01AY.f"
	*info = -13;
#line 171 "NF01AY.f"
    }

/*     Return if there are illegal arguments. */

#line 175 "NF01AY.f"
    if (*info != 0) {
#line 176 "NF01AY.f"
	i__1 = -(*info);
#line 176 "NF01AY.f"
	xerbla_("NF01AY", &i__1, (ftnlen)6);
#line 177 "NF01AY.f"
	return 0;
#line 178 "NF01AY.f"
    }

/*     Quick return if possible. */

#line 182 "NF01AY.f"
    if (min(*nsmp,*l) == 0) {
#line 182 "NF01AY.f"
	return 0;
#line 182 "NF01AY.f"
    }

/*     Set parameters to avoid overflows and increase accuracy for */
/*     extreme values. */

#line 188 "NF01AY.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 189 "NF01AY.f"
    bignum = 1. / smlnum;
#line 190 "NF01AY.f"
    dlabad_(&smlnum, &bignum);
#line 191 "NF01AY.f"
    smlnum = log(smlnum);
#line 192 "NF01AY.f"
    bignum = log(bignum);

#line 194 "NF01AY.f"
    ws = *nz * nn + 1;
#line 195 "NF01AY.f"
    ib = ws + nn - 1;
#line 196 "NF01AY.f"
    lk = 0;
#line 197 "NF01AY.f"
    if (min(*nz,nn) == 0) {
#line 198 "NF01AY.f"
	nv = 2;
#line 199 "NF01AY.f"
    } else {
#line 200 "NF01AY.f"
	nv = (*ldwork - nn) / nn;
#line 201 "NF01AY.f"
    }

#line 203 "NF01AY.f"
    if (nv > 2) {
#line 204 "NF01AY.f"
	mf = *nsmp / nv * nv;
#line 205 "NF01AY.f"
	last = *nsmp % nv != 0;

/*        Some BLAS 3 calculations can be used. */

#line 209 "NF01AY.f"
	i__1 = *l - 1;
#line 209 "NF01AY.f"
	for (k = 0; k <= i__1; ++k) {
#line 210 "NF01AY.f"
	    tmp = wb[ib + nn + 1 + lk];

#line 212 "NF01AY.f"
	    i__2 = nn;
#line 212 "NF01AY.f"
	    for (j = 1; j <= i__2; ++j) {
#line 213 "NF01AY.f"
		dwork[j] = wb[ib + j + lk] * 2.;
#line 214 "NF01AY.f"
/* L10: */
#line 214 "NF01AY.f"
	    }

#line 216 "NF01AY.f"
	    i__2 = mf;
#line 216 "NF01AY.f"
	    i__3 = nv;
#line 216 "NF01AY.f"
	    for (i__ = 1; i__3 < 0 ? i__ >= i__2 : i__ <= i__2; i__ += i__3) {

/*              Compute -2*[w1 w2 ... wn]'*Z', where */
/*              Z = [z(i)';...; z(i+NV-1)']. */

#line 221 "NF01AY.f"
		dgemm_("Transpose", "Transpose", &nn, &nv, nz, &c_b10, &wb[lk 
			+ 1], nz, &z__[i__ + z_dim1], ldz, &c_b11, &dwork[nn 
			+ 1], &nn, (ftnlen)9, (ftnlen)9);
#line 224 "NF01AY.f"
		lj = nn;

#line 226 "NF01AY.f"
		i__4 = nv;
#line 226 "NF01AY.f"
		for (m = 1; m <= i__4; ++m) {
#line 227 "NF01AY.f"
		    i__5 = nn;
#line 227 "NF01AY.f"
		    for (j = 1; j <= i__5; ++j) {

/*                    Compute tanh(wj'*z(i) + bj), j = 1:n. */

#line 231 "NF01AY.f"
			++lj;
#line 232 "NF01AY.f"
			df = dwork[lj] - dwork[j];
#line 233 "NF01AY.f"
			if (abs(df) >= bignum) {
#line 234 "NF01AY.f"
			    if (df > 0.) {
#line 235 "NF01AY.f"
				dwork[lj] = -1.;
#line 236 "NF01AY.f"
			    } else {
#line 237 "NF01AY.f"
				dwork[lj] = 1.;
#line 238 "NF01AY.f"
			    }
#line 239 "NF01AY.f"
			} else if (abs(df) <= smlnum) {
#line 240 "NF01AY.f"
			    dwork[lj] = 0.;
#line 241 "NF01AY.f"
			} else {
#line 242 "NF01AY.f"
			    dwork[lj] = 2. / (exp(df) + 1.) - 1.;
#line 243 "NF01AY.f"
			}
#line 244 "NF01AY.f"
/* L20: */
#line 244 "NF01AY.f"
		    }

#line 246 "NF01AY.f"
/* L30: */
#line 246 "NF01AY.f"
		}

#line 248 "NF01AY.f"
		y[i__ + (k + 1) * y_dim1] = tmp;
#line 249 "NF01AY.f"
		i__4 = nv - 1;
#line 249 "NF01AY.f"
		dcopy_(&i__4, &y[i__ + (k + 1) * y_dim1], &c__0, &y[i__ + 1 + 
			(k + 1) * y_dim1], &c__1);
#line 250 "NF01AY.f"
		dgemv_("Transpose", &nn, &nv, &c_b17, &dwork[nn + 1], &nn, &
			wb[ws + lk], &c__1, &c_b17, &y[i__ + (k + 1) * y_dim1]
			, &c__1, (ftnlen)9);
#line 252 "NF01AY.f"
/* L40: */
#line 252 "NF01AY.f"
	    }

#line 254 "NF01AY.f"
	    if (last) {

/*              Process the last samples. */

#line 258 "NF01AY.f"
		nv = *nsmp - mf;
#line 259 "NF01AY.f"
		i__ = mf + 1;

/*              Compute -2*[w1 w2 ... wn]'*Z', where */
/*              Z = [z(i)';...; z(NSMP)']. */

#line 264 "NF01AY.f"
		dgemm_("Transpose", "Transpose", &nn, &nv, nz, &c_b10, &wb[lk 
			+ 1], nz, &z__[i__ + z_dim1], ldz, &c_b11, &dwork[nn 
			+ 1], &nn, (ftnlen)9, (ftnlen)9);
#line 267 "NF01AY.f"
		lj = nn;

#line 269 "NF01AY.f"
		i__3 = nv;
#line 269 "NF01AY.f"
		for (m = 1; m <= i__3; ++m) {
#line 270 "NF01AY.f"
		    i__2 = nn;
#line 270 "NF01AY.f"
		    for (j = 1; j <= i__2; ++j) {

/*                    Compute tanh(wj'*z(i) + bj), j = 1:n. */

#line 274 "NF01AY.f"
			++lj;
#line 275 "NF01AY.f"
			df = dwork[lj] - dwork[j];
#line 276 "NF01AY.f"
			if (abs(df) >= bignum) {
#line 277 "NF01AY.f"
			    if (df > 0.) {
#line 278 "NF01AY.f"
				dwork[lj] = -1.;
#line 279 "NF01AY.f"
			    } else {
#line 280 "NF01AY.f"
				dwork[lj] = 1.;
#line 281 "NF01AY.f"
			    }
#line 282 "NF01AY.f"
			} else if (abs(df) <= smlnum) {
#line 283 "NF01AY.f"
			    dwork[lj] = 0.;
#line 284 "NF01AY.f"
			} else {
#line 285 "NF01AY.f"
			    dwork[lj] = 2. / (exp(df) + 1.) - 1.;
#line 286 "NF01AY.f"
			}
#line 287 "NF01AY.f"
/* L50: */
#line 287 "NF01AY.f"
		    }

#line 289 "NF01AY.f"
/* L60: */
#line 289 "NF01AY.f"
		}

#line 291 "NF01AY.f"
		y[i__ + (k + 1) * y_dim1] = tmp;
#line 292 "NF01AY.f"
		if (nv > 1) {
#line 292 "NF01AY.f"
		    i__3 = nv - 1;
#line 292 "NF01AY.f"
		    dcopy_(&i__3, &y[i__ + (k + 1) * y_dim1], &c__0, &y[i__ + 
			    1 + (k + 1) * y_dim1], &c__1);
#line 292 "NF01AY.f"
		}
#line 294 "NF01AY.f"
		dgemv_("Transpose", &nn, &nv, &c_b17, &dwork[nn + 1], &nn, &
			wb[ws + lk], &c__1, &c_b17, &y[i__ + (k + 1) * y_dim1]
			, &c__1, (ftnlen)9);
#line 296 "NF01AY.f"
	    }

#line 298 "NF01AY.f"
	    lk += ldwb;
#line 299 "NF01AY.f"
/* L70: */
#line 299 "NF01AY.f"
	}

#line 301 "NF01AY.f"
    } else {

/*        BLAS 2 calculations only can be used. */

#line 305 "NF01AY.f"
	i__1 = *l - 1;
#line 305 "NF01AY.f"
	for (k = 0; k <= i__1; ++k) {
#line 306 "NF01AY.f"
	    tmp = wb[ib + nn + 1 + lk];

#line 308 "NF01AY.f"
	    i__3 = nn;
#line 308 "NF01AY.f"
	    for (j = 1; j <= i__3; ++j) {
#line 309 "NF01AY.f"
		dwork[j] = wb[ib + j + lk] * 2.;
#line 310 "NF01AY.f"
/* L80: */
#line 310 "NF01AY.f"
	    }

#line 312 "NF01AY.f"
	    i__3 = *nsmp;
#line 312 "NF01AY.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {

/*              Compute -2*[w1 w2 ... wn]'*z(i). */

#line 316 "NF01AY.f"
		if (*nz == 0) {
#line 317 "NF01AY.f"
		    dwork[nn + 1] = 0.;
#line 318 "NF01AY.f"
		    dcopy_(&nn, &dwork[nn + 1], &c__0, &dwork[nn + 1], &c__1);
#line 319 "NF01AY.f"
		} else {
#line 320 "NF01AY.f"
		    dgemv_("Transpose", nz, &nn, &c_b10, &wb[lk + 1], nz, &
			    z__[i__ + z_dim1], ldz, &c_b11, &dwork[nn + 1], &
			    c__1, (ftnlen)9);
#line 322 "NF01AY.f"
		}

#line 324 "NF01AY.f"
		i__2 = nn << 1;
#line 324 "NF01AY.f"
		for (j = nn + 1; j <= i__2; ++j) {

/*                 Compute tanh(wj'*z(i) + bj), j = 1:n. */

#line 328 "NF01AY.f"
		    df = dwork[j] - dwork[j - nn];
#line 329 "NF01AY.f"
		    if (abs(df) >= bignum) {
#line 330 "NF01AY.f"
			if (df > 0.) {
#line 331 "NF01AY.f"
			    dwork[j] = -1.;
#line 332 "NF01AY.f"
			} else {
#line 333 "NF01AY.f"
			    dwork[j] = 1.;
#line 334 "NF01AY.f"
			}
#line 335 "NF01AY.f"
		    } else if (abs(df) <= smlnum) {
#line 336 "NF01AY.f"
			dwork[j] = 0.;
#line 337 "NF01AY.f"
		    } else {
#line 338 "NF01AY.f"
			dwork[j] = 2. / (exp(df) + 1.) - 1.;
#line 339 "NF01AY.f"
		    }
#line 340 "NF01AY.f"
/* L90: */
#line 340 "NF01AY.f"
		}

#line 342 "NF01AY.f"
		y[i__ + (k + 1) * y_dim1] = ddot_(&nn, &wb[ws + lk], &c__1, &
			dwork[nn + 1], &c__1) + tmp;
#line 344 "NF01AY.f"
/* L100: */
#line 344 "NF01AY.f"
	    }

#line 346 "NF01AY.f"
	    lk += ldwb;
#line 347 "NF01AY.f"
/* L110: */
#line 347 "NF01AY.f"
	}

#line 349 "NF01AY.f"
    }
#line 350 "NF01AY.f"
    return 0;

/* *** Last line of NF01AY *** */
} /* nf01ay_ */

