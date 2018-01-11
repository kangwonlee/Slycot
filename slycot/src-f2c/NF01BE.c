#line 1 "NF01BE.f"
/* NF01BE.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BE.f"
/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b3 = -1.;

/* Subroutine */ int nf01be_(integer *iflag, integer *nsmp, integer *n, 
	integer *ipar, integer *lipar, doublereal *z__, integer *ldz, 
	doublereal *y, integer *ldy, doublereal *x, integer *nfevl, 
	doublereal *e, doublereal *j, integer *ldj, doublereal *dwork, 
	integer *ldwork, integer *info)
{
    /* System generated locals */
    integer j_dim1, j_offset, y_dim1, y_offset, z_dim1, z_offset, i__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static doublereal err;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int nf01ay_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), nf01by_(char *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), daxpy_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);

    /* Fortran I/O blocks */
    static cilist io___2 = { 0, 6, 0, "(' Norm of current error = ', D15.6)", 
	    0 };



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

/*     This is the FCN routine for optimizing the parameters of the */
/*     nonlinear part of a Wiener system (initialization phase), using */
/*     SLICOT Library routine MD03BD. See the argument FCN in the */
/*     routine MD03BD for the description of parameters. Note that */
/*     NF01BE is called for each output of the Wiener system. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. CJTE is initialized to avoid the calculation of J'*e .. */
/*     .. NOUT is the unit number for printing intermediate results .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

#line 56 "NF01BE.f"
    /* Parameter adjustments */
#line 56 "NF01BE.f"
    --ipar;
#line 56 "NF01BE.f"
    z_dim1 = *ldz;
#line 56 "NF01BE.f"
    z_offset = 1 + z_dim1;
#line 56 "NF01BE.f"
    z__ -= z_offset;
#line 56 "NF01BE.f"
    y_dim1 = *ldy;
#line 56 "NF01BE.f"
    y_offset = 1 + y_dim1;
#line 56 "NF01BE.f"
    y -= y_offset;
#line 56 "NF01BE.f"
    --x;
#line 56 "NF01BE.f"
    --e;
#line 56 "NF01BE.f"
    j_dim1 = *ldj;
#line 56 "NF01BE.f"
    j_offset = 1 + j_dim1;
#line 56 "NF01BE.f"
    j -= j_offset;
#line 56 "NF01BE.f"
    --dwork;
#line 56 "NF01BE.f"

#line 56 "NF01BE.f"
    /* Function Body */
#line 56 "NF01BE.f"
    *info = 0;
#line 57 "NF01BE.f"
    if (*iflag == 1) {

/*        Call NF01AY to compute the output y of the Wiener system (in E) */
/*        and then the error functions (also in E). The array Z must */
/*        contain the output of the linear part of the Wiener system, and */
/*        Y must contain the original output Y of the Wiener system. */
/*        IPAR(2) must contain the number of outputs. */
/*        Workspace: need:    2*NN, NN = IPAR(3) (number of neurons); */
/*                   prefer:  larger. */

#line 67 "NF01BE.f"
	i__1 = *lipar - 2;
#line 67 "NF01BE.f"
	nf01ay_(nsmp, &ipar[2], &c__1, &ipar[3], &i__1, &x[1], n, &z__[
		z_offset], ldz, &e[1], nsmp, &dwork[1], ldwork, info);
#line 69 "NF01BE.f"
	daxpy_(nsmp, &c_b3, &y[y_offset], &c__1, &e[1], &c__1);
#line 70 "NF01BE.f"
	dwork[1] = (doublereal) (ipar[3] << 1);

#line 72 "NF01BE.f"
    } else if (*iflag == 2) {

/*        Call NF01BY to compute the Jacobian in a compressed form. */
/*        IPAR(2), IPAR(3) must have the same content as for IFLAG = 1. */
/*        Workspace: need:    0. */

#line 78 "NF01BE.f"
	i__1 = *lipar - 2;
#line 78 "NF01BE.f"
	nf01by_("N", nsmp, &ipar[2], &c__1, &ipar[3], &i__1, &x[1], n, &z__[
		z_offset], ldz, &e[1], &j[j_offset], ldj, &dwork[1], &dwork[1]
		, ldwork, info, (ftnlen)1);
#line 80 "NF01BE.f"
	*nfevl = 0;
#line 81 "NF01BE.f"
	dwork[1] = 0.;

#line 83 "NF01BE.f"
    } else if (*iflag == 3) {

/*        Set the parameter LDJ, the length of the array J, and the sizes */
/*        of the workspace for FCN (IFLAG = 1 or 2), QRFACT and LMPARM. */

#line 88 "NF01BE.f"
	*ldj = *nsmp;
#line 89 "NF01BE.f"
	ipar[1] = *nsmp * *n;
#line 90 "NF01BE.f"
	ipar[2] = ipar[3] << 1;
#line 91 "NF01BE.f"
	ipar[3] = 0;
#line 92 "NF01BE.f"
	ipar[4] = (*n << 2) + 1;
#line 93 "NF01BE.f"
	ipar[5] = *n << 2;

#line 95 "NF01BE.f"
    } else if (*iflag == 0) {

/*        Special call for printing intermediate results. */

#line 99 "NF01BE.f"
	err = dnrm2_(nsmp, &e[1], &c__1);
#line 100 "NF01BE.f"
	s_wsfe(&io___2);
#line 100 "NF01BE.f"
	do_fio(&c__1, (char *)&err, (ftnlen)sizeof(doublereal));
#line 100 "NF01BE.f"
	e_wsfe();
#line 101 "NF01BE.f"
    }
#line 102 "NF01BE.f"
    return 0;

/* *** Last line of NF01BE *** */
} /* nf01be_ */

