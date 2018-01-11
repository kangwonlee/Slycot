#line 1 "NF01AD.f"
/* NF01AD.f -- translated by f2c (version 20100827).
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

#line 1 "NF01AD.f"
/* Subroutine */ int nf01ad_(integer *nsmp, integer *m, integer *l, integer *
	ipar, integer *lipar, doublereal *x, integer *lx, doublereal *u, 
	integer *ldu, doublereal *y, integer *ldy, doublereal *dwork, integer 
	*ldwork, integer *info)
{
    /* System generated locals */
    integer u_dim1, u_offset, y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    static integer n, z__, ac, bd, nn, ix, jw, ldac, lths, nths;
    extern /* Subroutine */ int nf01ay_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), tf01mx_(integer *, integer *, integer *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *, integer *), 
	    tb01vy_(char *, integer *, integer *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, integer *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), xerbla_(char *, 
	    integer *, ftnlen);


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

/*     To calculate the output y of the Wiener system */

/*        x(t+1) = A*x(t) + B*u(t) */
/*        z(t)   = C*x(t) + D*u(t), */

/*        y(t)   = f(z(t),wb(1:L)), */

/*     where t = 1, 2, ..., NSMP, and f is a nonlinear function, */
/*     evaluated by the SLICOT Library routine NF01AY. The parameter */
/*     vector X is partitioned as X = ( wb(1), ..., wb(L), theta ), */
/*     where wb(i), i = 1:L, correspond to the nonlinear part, theta */
/*     corresponds to the linear part, and the notation is fully */
/*     described below. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     NSMP    (input) INTEGER */
/*             The number of training samples.  NSMP >= 0. */

/*     M       (input) INTEGER */
/*             The length of each input sample.  M >= 0. */

/*     L       (input) INTEGER */
/*             The length of each output sample.  L >= 0. */

/*     IPAR    (input) INTEGER array, dimension (LIPAR) */
/*             The integer parameters needed. */
/*             IPAR(1)  must contain the order of the linear part, */
/*                      referred to as N below.  N >= 0. */
/*             IPAR(2)  must contain the number of neurons for the */
/*                      nonlinear part, referred to as NN below. */
/*                      NN >= 0. */

/*     LIPAR   (input) INTEGER */
/*             The length of IPAR.  LIPAR >= 2. */

/*     X       (input) DOUBLE PRECISION array, dimension (LX) */
/*             The parameter vector, partitioned as */
/*             X = (wb(1), ..., wb(L), theta), where the vectors */
/*             wb(i), of length NN*(L+2)+1, are parameters for the */
/*             static nonlinearity, which is simulated by the */
/*             SLICOT Library routine NF01AY. See the documentation of */
/*             NF01AY for further details. The vector theta, of length */
/*             N*(M + L + 1) + L*M, represents the matrices A, B, C, */
/*             D and x(1), and it can be retrieved from these matrices */
/*             by SLICOT Library routine TB01VD and retranslated by */
/*             TB01VY. */

/*     LX      (input) INTEGER */
/*             The length of the array X. */
/*             LX >= ( NN*(L+2)+1 )*L + N*(M + L + 1) + L*M. */

/*     U       (input) DOUBLE PRECISION array, dimension (LDU, M) */
/*             The leading NSMP-by-M part of this array must contain the */
/*             set of input samples, */
/*             U = ( U(1,1),...,U(1,M); ...; U(NSMP,1),...,U(NSMP,M) ). */

/*     LDU     INTEGER */
/*             The leading dimension of the array U.  LDU >= MAX(1,NSMP). */

/*     Y       (output) DOUBLE PRECISION array, dimension (LDY, L) */
/*             The leading NSMP-by-L part of this array contains the */
/*             simulated output. */

/*     LDY     INTEGER */
/*             The leading dimension of the array Y.  LDY >= MAX(1,NSMP). */

/*     Workspace */

/*     DWORK   DOUBLE PRECISION array, dimension (LDWORK) */

/*     LDWORK  INTEGER */
/*             The length of the array DWORK. */
/*             LDWORK >= NSMP*L + MAX( 2*NN, (N + L)*(N + M) + 2*N + */
/*                                     MAX( N*(N + L), N + M + L ) ) */
/*                                                              if M > 0; */
/*             LDWORK >= NSMP*L + MAX( 2*NN, (N + L)*N + 2*N + */
/*                                     MAX( N*(N + L), L ) ),   if M = 0. */
/*             A larger value of LDWORK could improve the efficiency. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */
/*     METHOD */

/*     BLAS routines are used for the matrix-vector multiplications and */
/*     the routine NF01AY is called for the calculation of the nonlinear */
/*     function. */

/*     CONTRIBUTORS */

/*     A. Riedel, R. Schneider, Chemnitz University of Technology, */
/*     Mar. 2001, during a stay at University of Twente, NL. */

/*     REVISIONS */

/*     V. Sima, Research Institute for Informatics, Bucharest, Mar. 2001, */
/*     Dec. 2001. */

/*     KEYWORDS */

/*     Nonlinear system, output normal form, simulation, state-space */
/*     representation, Wiener system. */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. */
/*     .. Executable Statements .. */

#line 149 "NF01AD.f"
    /* Parameter adjustments */
#line 149 "NF01AD.f"
    --ipar;
#line 149 "NF01AD.f"
    --x;
#line 149 "NF01AD.f"
    u_dim1 = *ldu;
#line 149 "NF01AD.f"
    u_offset = 1 + u_dim1;
#line 149 "NF01AD.f"
    u -= u_offset;
#line 149 "NF01AD.f"
    y_dim1 = *ldy;
#line 149 "NF01AD.f"
    y_offset = 1 + y_dim1;
#line 149 "NF01AD.f"
    y -= y_offset;
#line 149 "NF01AD.f"
    --dwork;
#line 149 "NF01AD.f"

#line 149 "NF01AD.f"
    /* Function Body */
#line 149 "NF01AD.f"
    *info = 0;
#line 150 "NF01AD.f"
    if (*nsmp < 0) {
#line 151 "NF01AD.f"
	*info = -1;
#line 152 "NF01AD.f"
    } else if (*m < 0) {
#line 153 "NF01AD.f"
	*info = -2;
#line 154 "NF01AD.f"
    } else if (*l < 0) {
#line 155 "NF01AD.f"
	*info = -3;
#line 156 "NF01AD.f"
    } else if (*lipar < 2) {
#line 157 "NF01AD.f"
	*info = -5;
#line 158 "NF01AD.f"
    } else {

#line 160 "NF01AD.f"
	n = ipar[1];
#line 161 "NF01AD.f"
	nn = ipar[2];
#line 162 "NF01AD.f"
	ldac = n + *l;
#line 163 "NF01AD.f"
	nths = (nn * (*l + 2) + 1) * *l;
#line 164 "NF01AD.f"
	lths = n * (*m + *l + 1) + *l * *m;

#line 166 "NF01AD.f"
	if (n < 0 || nn < 0) {
#line 167 "NF01AD.f"
	    *info = -4;
#line 168 "NF01AD.f"
	} else if (*lx < nths + lths) {
#line 169 "NF01AD.f"
	    *info = -7;
#line 170 "NF01AD.f"
	} else if (*ldu < max(1,*nsmp)) {
#line 171 "NF01AD.f"
	    *info = -9;
#line 172 "NF01AD.f"
	} else if (*ldy < max(1,*nsmp)) {
#line 173 "NF01AD.f"
	    *info = -11;
#line 174 "NF01AD.f"
	} else {
#line 175 "NF01AD.f"
	    if (*m > 0) {
/* Computing MAX */
#line 176 "NF01AD.f"
		i__1 = n * ldac, i__2 = n + *m + *l;
#line 176 "NF01AD.f"
		jw = max(i__1,i__2);
#line 177 "NF01AD.f"
	    } else {
/* Computing MAX */
#line 178 "NF01AD.f"
		i__1 = n * ldac;
#line 178 "NF01AD.f"
		jw = max(i__1,*l);
#line 179 "NF01AD.f"
	    }
/* Computing MAX */
#line 180 "NF01AD.f"
	    i__1 = nn << 1, i__2 = ldac * (n + *m) + (n << 1) + jw;
#line 180 "NF01AD.f"
	    if (*ldwork < *nsmp * *l + max(i__1,i__2)) {
#line 180 "NF01AD.f"
		*info = -13;
#line 180 "NF01AD.f"
	    }
#line 183 "NF01AD.f"
	}
#line 184 "NF01AD.f"
    }

/*     Return if there are illegal arguments. */

#line 188 "NF01AD.f"
    if (*info != 0) {
#line 189 "NF01AD.f"
	i__1 = -(*info);
#line 189 "NF01AD.f"
	xerbla_("NF01AD", &i__1, (ftnlen)6);
#line 190 "NF01AD.f"
	return 0;
#line 191 "NF01AD.f"
    }

/*     Quick return if possible. */

#line 195 "NF01AD.f"
    if (min(*nsmp,*l) == 0) {
#line 195 "NF01AD.f"
	return 0;
#line 195 "NF01AD.f"
    }

/*     Compute the output of the linear part. */
/*     Workspace: need   NSMP*L + (N + L)*(N + M) + N + N*(N + L + 1). */
/*     (NSMP*L locations are reserved for the output of the linear part.) */

#line 202 "NF01AD.f"
    z__ = 1;
#line 203 "NF01AD.f"
    ac = z__ + *nsmp * *l;
#line 204 "NF01AD.f"
    bd = ac + ldac * n;
#line 205 "NF01AD.f"
    ix = bd + ldac * *m;
#line 206 "NF01AD.f"
    jw = ix + n;

#line 208 "NF01AD.f"
    i__1 = *ldwork - jw + 1;
#line 208 "NF01AD.f"
    tb01vy_("Apply", &n, m, l, &x[nths + 1], &lths, &dwork[ac], &ldac, &dwork[
	    bd], &ldac, &dwork[ac + n], &ldac, &dwork[bd + n], &ldac, &dwork[
	    ix], &dwork[jw], &i__1, info, (ftnlen)5);

/*     Workspace: need   NSMP*L + (N + L)*(N + M) + 3*N + M + L, if M>0; */
/*                       NSMP*L + (N + L)*N + 2*N + L,           if M=0; */
/*                prefer larger. */

#line 216 "NF01AD.f"
    i__1 = *ldwork - jw + 1;
#line 216 "NF01AD.f"
    tf01mx_(&n, m, l, nsmp, &dwork[ac], &ldac, &u[u_offset], ldu, &dwork[ix], 
	    &dwork[z__], nsmp, &dwork[jw], &i__1, info);

/*     Simulate the static nonlinearity. */
/*     Workspace: need   NSMP*L + 2*NN; */
/*                prefer larger. */

#line 223 "NF01AD.f"
    jw = ac;
#line 224 "NF01AD.f"
    i__1 = *lipar - 1;
#line 224 "NF01AD.f"
    i__2 = *ldwork - jw + 1;
#line 224 "NF01AD.f"
    nf01ay_(nsmp, l, l, &ipar[2], &i__1, &x[1], &nths, &dwork[z__], nsmp, &y[
	    y_offset], ldy, &dwork[jw], &i__2, info);

#line 227 "NF01AD.f"
    return 0;

/* *** Last line of NF01AD *** */
} /* nf01ad_ */

