#line 1 "NF01BY.f"
/* NF01BY.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BY.f"
/* Table of constant values */

static integer c__0 = 0;
static integer c__1 = 1;
static doublereal c_b15 = -2.;
static doublereal c_b22 = 1.;
static doublereal c_b24 = 0.;

/* Subroutine */ int nf01by_(char *cjte, integer *nsmp, integer *nz, integer *
	l, integer *ipar, integer *lipar, doublereal *wb, integer *lwb, 
	doublereal *z__, integer *ldz, doublereal *e, doublereal *j, integer *
	ldj, doublereal *jte, doublereal *dwork, integer *ldwork, integer *
	info, ftnlen cjte_len)
{
    /* System generated locals */
    integer j_dim1, j_offset, z_dim1, z_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    double log(doublereal), exp(doublereal);

    /* Local variables */
    static integer i__, k, m, ib, di, nn, ws, bp1, nwb;
    static doublereal tmp;
    static logical wjte;
    extern /* Subroutine */ int dgemm_(char *, char *, integer *, integer *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, ftnlen, ftnlen);
    extern logical lsame_(char *, char *, ftnlen, ftnlen);
    extern /* Subroutine */ int dgemv_(char *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, integer *, ftnlen), dcopy_(integer *, 
	    doublereal *, integer *, doublereal *, integer *), dlabad_(
	    doublereal *, doublereal *);
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

/*     To compute the Jacobian of the error function for a neural network */
/*     of the structure */

/*             - tanh(w1*z+b1) - */
/*           /      :            \ */
/*         z ---    :          --- sum(ws(i)*...)+ b(n+1)  --- y, */
/*           \      :            / */
/*             - tanh(wn*z+bn) - */

/*     for the single-output case. The Jacobian has the form */

/*                d e(1)  / d WB(1)   ...    d e(1)  / d WB(NWB) */
/*         J =            :                          :           , */
/*              d e(NSMP) / d WB(1)   ...  d e(NSMP) / d WB(NWB) */

/*     where e(z) is the error function, WB is the set of weights and */
/*     biases of the network (for the considered output), and NWB is */
/*     the number of elements of this set, NWB = IPAR(1)*(NZ+2)+1 */
/*     (see below). */

/*     In the multi-output case, this routine should be called for each */
/*     output. */

/*     NOTE: this routine must have the same arguments as SLICOT Library */
/*     routine NF01BD. */

/*     ARGUMENTS */

/*     Mode Parameters */

/*     CJTE    CHARACTER*1 */
/*             Specifies whether the matrix-vector product J'*e should be */
/*             computed or not, as follows: */
/*             = 'C' :  compute J'*e; */
/*             = 'N' :  do not compute J'*e. */

/*     Input/Output Parameters */

/*     NSMP    (input) INTEGER */
/*             The number of training samples.  NSMP >= 0. */

/*     NZ      (input) INTEGER */
/*             The length of each input sample.  NZ >= 0. */

/*     L       (input) INTEGER */
/*             The length of each output sample. */
/*             Currently, L must be 1. */

/*     IPAR    (input/output) INTEGER array, dimension (LIPAR) */
/*             The integer parameters needed. */
/*             On entry, the first element of this array must contain */
/*             a value related to the number of neurons, n; specifically, */
/*             n = abs(IPAR(1)), since setting IPAR(1) < 0 has a special */
/*             meaning (see below). */
/*             On exit, if IPAR(1) < 0 on entry, then no computations are */
/*             performed, except the needed tests on input parameters, */
/*             but the following values are returned: */
/*             IPAR(1) contains the length of the array J, LJ; */
/*             LDJ     contains the leading dimension of array J. */
/*             Otherwise, IPAR(1) and LDJ are unchanged on exit. */

/*     LIPAR   (input) INTEGER */
/*             The length of the vector IPAR.  LIPAR >= 1. */

/*     WB      (input) DOUBLE PRECISION array, dimension (LWB) */
/*             The leading NWB = IPAR(1)*(NZ+2)+1 part of this array */
/*             must contain the weights and biases of the network, */
/*             WB = ( w(1,1), ..., w(1,NZ), ..., w(n,1), ...,  w(n,NZ), */
/*                    ws(1), ..., ws(n), b(1), ..., b(n+1) ), */
/*             where w(i,j) are the weights of the hidden layer, */
/*             ws(i) are the weights of the linear output layer and */
/*             b(i) are the biases. */

/*     LWB     (input) INTEGER */
/*             The length of array WB.  LWB >= NWB. */

/*     Z       (input) DOUBLE PRECISION array, dimension (LDZ, NZ) */
/*             The leading NSMP-by-NZ part of this array must contain the */
/*             set of input samples, */
/*             Z = ( Z(1,1),...,Z(1,NZ); ...; Z(NSMP,1),...,Z(NSMP,NZ) ). */

/*     LDZ     INTEGER */
/*             The leading dimension of array Z.  LDZ >= MAX(1,NSMP). */

/*     E       (input) DOUBLE PRECISION array, dimension (NSMP) */
/*             If CJTE = 'C', this array must contain the error vector e. */
/*             If CJTE = 'N', this array is not referenced. */

/*     J       (output) DOUBLE PRECISION array, dimension (LDJ, NWB) */
/*             The leading NSMP-by-NWB part of this array contains the */
/*             Jacobian of the error function. */

/*     LDJ     INTEGER */
/*             The leading dimension of array J.  LDJ >= MAX(1,NSMP). */
/*             Note that LDJ is an input parameter, except for */
/*             IPAR(1) < 0 on entry, when it is an output parameter. */

/*     JTE     (output) DOUBLE PRECISION array, dimension (NWB) */
/*             If CJTE = 'C', this array contains the matrix-vector */
/*             product J'*e. */
/*             If CJTE = 'N', this array is not referenced. */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */
/*             This argument is included for combatibility with SLICOT */
/*             Library routine NF01BD. */

/*     LDWORK  INTEGER */
/*             Normally, the length of the array DWORK.  LDWORK >= 0. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The Jacobian is computed analytically. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Oct. 2000, during a stay at University of Twente, NL. */
/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001. */

/*     REVISIONS */

/*     - */

/*     KEYWORDS */

/*     Input output description, neural network, nonlinear system, */
/*     optimization, system response. */

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

#line 188 "NF01BY.f"
    /* Parameter adjustments */
#line 188 "NF01BY.f"
    --ipar;
#line 188 "NF01BY.f"
    --wb;
#line 188 "NF01BY.f"
    z_dim1 = *ldz;
#line 188 "NF01BY.f"
    z_offset = 1 + z_dim1;
#line 188 "NF01BY.f"
    z__ -= z_offset;
#line 188 "NF01BY.f"
    --e;
#line 188 "NF01BY.f"
    j_dim1 = *ldj;
#line 188 "NF01BY.f"
    j_offset = 1 + j_dim1;
#line 188 "NF01BY.f"
    j -= j_offset;
#line 188 "NF01BY.f"
    --jte;
#line 188 "NF01BY.f"
    --dwork;
#line 188 "NF01BY.f"

#line 188 "NF01BY.f"
    /* Function Body */
#line 188 "NF01BY.f"
    wjte = lsame_(cjte, "C", (ftnlen)1, (ftnlen)1);
#line 189 "NF01BY.f"
    *info = 0;
#line 190 "NF01BY.f"
    nn = ipar[1];
#line 191 "NF01BY.f"
    nwb = nn * (*nz + 2) + 1;
#line 192 "NF01BY.f"
    if (! (wjte || lsame_(cjte, "N", (ftnlen)1, (ftnlen)1))) {
#line 193 "NF01BY.f"
	*info = -1;
#line 194 "NF01BY.f"
    } else if (*nsmp < 0) {
#line 195 "NF01BY.f"
	*info = -2;
#line 196 "NF01BY.f"
    } else if (*nz < 0) {
#line 197 "NF01BY.f"
	*info = -3;
#line 198 "NF01BY.f"
    } else if (*l != 1) {
#line 199 "NF01BY.f"
	*info = -4;
#line 200 "NF01BY.f"
    } else if (*lipar < 1) {
#line 201 "NF01BY.f"
	*info = -6;
#line 202 "NF01BY.f"
    } else if (ipar[1] < 0) {
#line 203 "NF01BY.f"
	if (*info != 0) {
#line 204 "NF01BY.f"
	    i__1 = -(*info);
#line 204 "NF01BY.f"
	    xerbla_("NF01BY", &i__1, (ftnlen)6);
#line 205 "NF01BY.f"
	} else {
#line 206 "NF01BY.f"
	    ipar[1] = *nsmp * (abs(nn) * (*nz + 2) + 1);
#line 207 "NF01BY.f"
	    *ldj = *nsmp;
#line 208 "NF01BY.f"
	}
#line 209 "NF01BY.f"
	return 0;
#line 210 "NF01BY.f"
    } else if (*lwb < nwb) {
#line 211 "NF01BY.f"
	*info = -8;
#line 212 "NF01BY.f"
    } else if (*ldz < max(1,*nsmp)) {
#line 213 "NF01BY.f"
	*info = -10;
#line 214 "NF01BY.f"
    } else if (*ldj < max(1,*nsmp)) {
#line 215 "NF01BY.f"
	*info = -13;
#line 216 "NF01BY.f"
    }

/*     Return if there are illegal arguments. */

#line 220 "NF01BY.f"
    if (*info != 0) {
#line 221 "NF01BY.f"
	i__1 = -(*info);
#line 221 "NF01BY.f"
	xerbla_("NF01BY", &i__1, (ftnlen)6);
#line 222 "NF01BY.f"
	return 0;
#line 223 "NF01BY.f"
    }

/*     Quick return if possible. */

#line 227 "NF01BY.f"
    if (min(*nsmp,*nz) == 0) {
#line 227 "NF01BY.f"
	return 0;
#line 227 "NF01BY.f"
    }

/*     Set parameters to avoid overflows and increase accuracy for */
/*     extreme values. */

#line 233 "NF01BY.f"
    smlnum = dlamch_("Safe minimum", (ftnlen)12) / dlamch_("Precision", (
	    ftnlen)9);
#line 234 "NF01BY.f"
    bignum = 1. / smlnum;
#line 235 "NF01BY.f"
    dlabad_(&smlnum, &bignum);
#line 236 "NF01BY.f"
    smlnum = log(smlnum);
#line 237 "NF01BY.f"
    bignum = log(bignum);

#line 239 "NF01BY.f"
    ws = *nz * nn + 1;
#line 240 "NF01BY.f"
    ib = ws + nn;
#line 241 "NF01BY.f"
    bp1 = ib + nn;

#line 243 "NF01BY.f"
    j[bp1 * j_dim1 + 1] = 1.;
#line 244 "NF01BY.f"
    dcopy_(nsmp, &j[bp1 * j_dim1 + 1], &c__0, &j[bp1 * j_dim1 + 1], &c__1);

#line 246 "NF01BY.f"
    i__1 = nn - 1;
#line 246 "NF01BY.f"
    for (i__ = 0; i__ <= i__1; ++i__) {
#line 247 "NF01BY.f"
	dcopy_(nsmp, &wb[ib + i__], &c__0, &j[(ws + i__) * j_dim1 + 1], &c__1)
		;
#line 248 "NF01BY.f"
/* L10: */
#line 248 "NF01BY.f"
    }

#line 250 "NF01BY.f"
    dgemm_("NoTranspose", "NoTranspose", nsmp, &nn, nz, &c_b15, &z__[z_offset]
	    , ldz, &wb[1], nz, &c_b15, &j[ws * j_dim1 + 1], ldj, (ftnlen)11, (
	    ftnlen)11);
#line 252 "NF01BY.f"
    di = 1;

#line 254 "NF01BY.f"
    i__1 = nn - 1;
#line 254 "NF01BY.f"
    for (i__ = 0; i__ <= i__1; ++i__) {

#line 256 "NF01BY.f"
	i__2 = *nsmp;
#line 256 "NF01BY.f"
	for (k = 1; k <= i__2; ++k) {
#line 257 "NF01BY.f"
	    tmp = j[k + (ws + i__) * j_dim1];
#line 258 "NF01BY.f"
	    if (abs(tmp) >= bignum) {
#line 259 "NF01BY.f"
		if (tmp > 0.) {
#line 260 "NF01BY.f"
		    j[k + (ws + i__) * j_dim1] = -1.;
#line 261 "NF01BY.f"
		} else {
#line 262 "NF01BY.f"
		    j[k + (ws + i__) * j_dim1] = 1.;
#line 263 "NF01BY.f"
		}
#line 264 "NF01BY.f"
	    } else if (abs(tmp) <= smlnum) {
#line 265 "NF01BY.f"
		j[k + (ws + i__) * j_dim1] = 0.;
#line 266 "NF01BY.f"
	    } else {
#line 267 "NF01BY.f"
		j[k + (ws + i__) * j_dim1] = 2. / (exp(tmp) + 1.) - 1.;
#line 268 "NF01BY.f"
	    }
/* Computing 2nd power */
#line 269 "NF01BY.f"
	    d__1 = j[k + (ws + i__) * j_dim1];
#line 269 "NF01BY.f"
	    j[k + (ib + i__) * j_dim1] = wb[ws + i__] * (1. - d__1 * d__1);
#line 270 "NF01BY.f"
/* L20: */
#line 270 "NF01BY.f"
	}

#line 272 "NF01BY.f"
	i__2 = *nz - 1;
#line 272 "NF01BY.f"
	for (k = 0; k <= i__2; ++k) {

#line 274 "NF01BY.f"
	    i__3 = *nsmp;
#line 274 "NF01BY.f"
	    for (m = 1; m <= i__3; ++m) {
#line 275 "NF01BY.f"
		j[m + (di + k) * j_dim1] = j[m + (ib + i__) * j_dim1] * z__[m 
			+ (k + 1) * z_dim1];
#line 276 "NF01BY.f"
/* L30: */
#line 276 "NF01BY.f"
	    }

#line 278 "NF01BY.f"
/* L40: */
#line 278 "NF01BY.f"
	}

#line 280 "NF01BY.f"
	di += *nz;
#line 281 "NF01BY.f"
/* L50: */
#line 281 "NF01BY.f"
    }

#line 283 "NF01BY.f"
    if (wjte) {

/*        Compute J'e. */

#line 287 "NF01BY.f"
	dgemv_("Transpose", nsmp, &nwb, &c_b22, &j[j_offset], ldj, &e[1], &
		c__1, &c_b24, &jte[1], &c__1, (ftnlen)9);
#line 289 "NF01BY.f"
    }

#line 291 "NF01BY.f"
    return 0;

/* *** Last line of NF01BY *** */
} /* nf01by_ */

