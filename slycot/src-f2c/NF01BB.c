#line 1 "NF01BB.f"
/* NF01BB.f -- translated by f2c (version 20100827).
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

#line 1 "NF01BB.f"
/* Table of constant values */

static doublereal c_b3 = -1.;
static integer c__1 = 1;

/* Subroutine */ int nf01bb_(integer *iflag, integer *nfun, integer *lx, 
	integer *ipar, integer *lipar, doublereal *u, integer *ldu, 
	doublereal *y, integer *ldy, doublereal *x, integer *nfevl, 
	doublereal *e, doublereal *j, integer *ldj, doublereal *jte, 
	doublereal *dwork, integer *ldwork, integer *info)
{
    /* System generated locals */
    integer j_dim1, j_offset, u_dim1, u_offset, y_dim1, y_offset, i__1, i__2, 
	    i__3, i__4;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, l, m, n, nn, st, bsn;
    static doublereal err;
    static integer nsmp;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int nf01ad_(integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *, 
	    integer *), nf01bd_(char *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, integer *, doublereal *, 
	    integer *, doublereal *, doublereal *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, ftnlen), daxpy_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    static integer jwork;

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 6, 0, "(' Norm of current error = ', D15.6)",
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

/*     This is the FCN routine for optimizing all parameters of a Wiener */
/*     system using SLICOT Library routine MD03AD. See the argument FCN */
/*     in the routine MD03AD for the description of parameters. */

/*     ****************************************************************** */

/*     .. Parameters .. */
/*     .. CJTE is initialized to activate the calculation of J'*e .. */
/*     .. NOUT is the unit number for printing intermediate results .. */
/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. External Subroutines .. */

/*     .. Executable Statements .. */

#line 55 "NF01BB.f"
    /* Parameter adjustments */
#line 55 "NF01BB.f"
    --ipar;
#line 55 "NF01BB.f"
    u_dim1 = *ldu;
#line 55 "NF01BB.f"
    u_offset = 1 + u_dim1;
#line 55 "NF01BB.f"
    u -= u_offset;
#line 55 "NF01BB.f"
    y_dim1 = *ldy;
#line 55 "NF01BB.f"
    y_offset = 1 + y_dim1;
#line 55 "NF01BB.f"
    y -= y_offset;
#line 55 "NF01BB.f"
    --x;
#line 55 "NF01BB.f"
    --e;
#line 55 "NF01BB.f"
    j_dim1 = *ldj;
#line 55 "NF01BB.f"
    j_offset = 1 + j_dim1;
#line 55 "NF01BB.f"
    j -= j_offset;
#line 55 "NF01BB.f"
    --jte;
#line 55 "NF01BB.f"
    --dwork;
#line 55 "NF01BB.f"

#line 55 "NF01BB.f"
    /* Function Body */
#line 55 "NF01BB.f"
    l = ipar[2];
#line 56 "NF01BB.f"
    m = ipar[5];
#line 57 "NF01BB.f"
    n = ipar[6];
#line 58 "NF01BB.f"
    if (l == 0) {
#line 59 "NF01BB.f"
	nsmp = *nfun;
#line 60 "NF01BB.f"
    } else {
#line 61 "NF01BB.f"
	nsmp = *nfun / l;
#line 62 "NF01BB.f"
    }

#line 64 "NF01BB.f"
    *info = 0;
#line 65 "NF01BB.f"
    if (*iflag == 1) {

/*        Call NF01AD to compute the output y of the Wiener system (in E) */
/*        and then the error functions (also in E). The array U must */
/*        contain the input to the linear part of the Wiener system, and */
/*        Y must contain the original output Y of the Wiener system. */
/*        IPAR(6) must contain the number of states of the linear part, n. */
/*        Workspace: need:    NFUN + MAX( 2*NN, (N + L)*(N + M) + 2*N + */
/*                                        MAX( N*(N + L), N + M + L ) ), */
/*                                                               if M>0, */
/*                            NFUN + MAX( 2*NN, (N + L)*N + 2*N + */
/*                                        MAX( N*(N + L), L ) ), if M=0, */
/*                            where NN = IPAR(7) (number of neurons); */
/*                   prefer:  larger. */

#line 80 "NF01BB.f"
	i__1 = *lipar - 2;
#line 80 "NF01BB.f"
	nf01ad_(&nsmp, &m, &l, &ipar[6], &i__1, &x[1], lx, &u[u_offset], ldu, 
		&e[1], &nsmp, &dwork[1], ldwork, info);

#line 83 "NF01BB.f"
	i__1 = l;
#line 83 "NF01BB.f"
	for (i__ = 1; i__ <= i__1; ++i__) {
#line 84 "NF01BB.f"
	    daxpy_(&nsmp, &c_b3, &y[i__ * y_dim1 + 1], &c__1, &e[(i__ - 1) * 
		    nsmp + 1], &c__1);
#line 85 "NF01BB.f"
/* L10: */
#line 85 "NF01BB.f"
	}

/* Computing MAX */
/* Computing MAX */
#line 87 "NF01BB.f"
	i__3 = n * (n + l), i__4 = n + m + l;
#line 87 "NF01BB.f"
	i__1 = ipar[7] << 1, i__2 = (n + l) * (n + m) + (n << 1) + max(i__3,
		i__4);
#line 87 "NF01BB.f"
	dwork[1] = (doublereal) (*nfun + max(i__1,i__2));

#line 90 "NF01BB.f"
    } else if (*iflag == 2) {

/*        Call NF01BD to compute the Jacobian in a compressed form. */
/*        Workspace: need:    2*NFUN + MAX( 2*NN, (N + L)*(N + M) + 2*N + */
/*                                          MAX( N*(N + L), N + M + L )), */
/*                                                              if M > 0, */
/*                            2*NFUN + MAX( 2*NN, (N + L)*N + 2*N + */
/*                                          MAX( N*(N + L), L ) ), */
/*                                                              if M = 0; */
/*                   prefer:  larger. */

#line 101 "NF01BB.f"
	i__1 = *lipar - 2;
#line 101 "NF01BB.f"
	nf01bd_("C", &nsmp, &m, &l, &ipar[6], &i__1, &x[1], lx, &u[u_offset], 
		ldu, &e[1], &j[j_offset], ldj, &jte[1], &dwork[1], ldwork, 
		info, (ftnlen)1);
#line 103 "NF01BB.f"
	*nfevl = ipar[6] * (m + l + 1) + l * m;
/* Computing MAX */
/* Computing MAX */
#line 104 "NF01BB.f"
	i__3 = n * (n + l), i__4 = n + m + l;
#line 104 "NF01BB.f"
	i__1 = ipar[7] << 1, i__2 = (n + l) * (n + m) + (n << 1) + max(i__3,
		i__4);
#line 104 "NF01BB.f"
	dwork[1] = (doublereal) ((*nfun << 1) + max(i__1,i__2));

#line 107 "NF01BB.f"
    } else if (*iflag == 3) {

/*        Set the parameter LDJ, the length of the array J, and the sizes */
/*        of the workspace for FCN (IFLAG = 1 or 2), and JTJ. */

#line 112 "NF01BB.f"
	st = ipar[1];
#line 113 "NF01BB.f"
	bsn = ipar[4];
#line 114 "NF01BB.f"
	nn = ipar[7];

#line 116 "NF01BB.f"
	*ldj = *nfun;
#line 117 "NF01BB.f"
	ipar[1] = *nfun * (bsn + st);
#line 118 "NF01BB.f"
	if (m > 0) {
/* Computing MAX */
#line 119 "NF01BB.f"
	    i__1 = n * (n + l), i__2 = n + m + l;
#line 119 "NF01BB.f"
	    jwork = max(i__1,i__2);
#line 120 "NF01BB.f"
	} else {
/* Computing MAX */
#line 121 "NF01BB.f"
	    i__1 = n * (n + l);
#line 121 "NF01BB.f"
	    jwork = max(i__1,l);
#line 122 "NF01BB.f"
	}
/* Computing MAX */
#line 123 "NF01BB.f"
	i__1 = (n + l) * (n + m) + (n << 1) + jwork, i__2 = nn << 1;
#line 123 "NF01BB.f"
	ipar[2] = *ldj + max(i__1,i__2);
#line 124 "NF01BB.f"
	ipar[3] = *ldj + ipar[2];
#line 125 "NF01BB.f"
	ipar[4] = 0;
#line 126 "NF01BB.f"
	ipar[5] = *nfun;

#line 128 "NF01BB.f"
    } else if (*iflag == 0) {

/*        Special call for printing intermediate results. */

#line 132 "NF01BB.f"
	err = dnrm2_(nfun, &e[1], &c__1);
#line 133 "NF01BB.f"
	s_wsfe(&io___11);
#line 133 "NF01BB.f"
	do_fio(&c__1, (char *)&err, (ftnlen)sizeof(doublereal));
#line 133 "NF01BB.f"
	e_wsfe();
#line 134 "NF01BB.f"
    }
#line 135 "NF01BB.f"
    return 0;

/* *** Last line of NF01BB *** */
} /* nf01bb_ */

