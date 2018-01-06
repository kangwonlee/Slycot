#line 1 "UD01MZ.f"
/* UD01MZ.f -- translated by f2c (version 20100827).
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

#line 1 "UD01MZ.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int ud01mz_(integer *m, integer *n, integer *l, integer *
	nout, doublecomplex *a, integer *lda, char *text, integer *info, 
	ftnlen text_len)
{
    /* Format strings */
    static char fmt_99996[] = "(1x,a,\002 (\002,i5,\002X\002,i5,\002)\002,/)";
    static char fmt_99999[] = "(7x,5(13x,i5,14x))";
    static char fmt_99997[] = "(1x,i5,2x,3(d15.7,sp,d15.7,s,\002i \002))";
    static char fmt_99998[] = "(\002 \002)";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Builtin functions */
    integer i_len(char *, ftnlen), s_wsfe(cilist *), do_fio(integer *, char *,
	     ftnlen), e_wsfe(void);

    /* Local variables */
    static integer i__, j, j1, j2, n1, jj, ltext;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lentxt;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 0, 0, fmt_99996, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_99999, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_99997, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_99998, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_99999, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_99997, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_99998, 0 };



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

/*     To print an M-by-N real matrix A row by row. The elements of A */
/*     are output to 7 significant figures. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     M       (input) INTEGER */
/*             The number of rows of matrix A to be printed.  M >= 1. */

/*     N       (input) INTEGER */
/*             The number of columns of matrix A to be printed.  N >= 1. */

/*     L       (input) INTEGER */
/*             The number of elements of matrix A to be printed per line. */
/*             1 <= L <= 3. */

/*     NOUT    (input) INTEGER */
/*             The output channel to which the results are sent. */
/*             NOUT >= 0. */

/*     A       (input) COMPLEX*16 array, dimension (LDA,N) */
/*             The leading M-by-N part of this array must contain the */
/*             matrix to be printed. */

/*     LDA     INTEGER */
/*             The leading dimension of array A.  LDA >= M. */

/*     TEXT    (input) CHARACTER*72. */
/*             Title caption of the matrix to be printed (up to a */
/*             maximum of 72 characters). For example, TEXT = 'Matrix A'. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     The routine first prints the contents of TEXT as a title, followed */
/*     by the elements of the matrix A such that */

/*     (i)  if N <= L, the leading M-by-N part is printed; */
/*     (ii) if N = k*L + p (where k,p > 0), then k M-by-L blocks of */
/*          consecutive columns of A are printed one after another */
/*          followed by one M-by-p block containing the last p columns */
/*          of A. */

/*     Row numbers are printed on the left of each row and a column */
/*     number appears on top of each complex column. */
/*     The routine uses 2 + (k + 1)*(m + 1) lines and 7 + 32*c positions */
/*     per line where c is the actual number of columns, (i.e. c = L */
/*     or c = p). */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, Oct. 1997. */
/*     Complex version: V. Sima, Research Institute for Informatics, */
/*     Bucharest, Dec. 2008. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */
/*     .. Executable Statements .. */

#line 112 "UD01MZ.f"
    /* Parameter adjustments */
#line 112 "UD01MZ.f"
    a_dim1 = *lda;
#line 112 "UD01MZ.f"
    a_offset = 1 + a_dim1;
#line 112 "UD01MZ.f"
    a -= a_offset;
#line 112 "UD01MZ.f"

#line 112 "UD01MZ.f"
    /* Function Body */
#line 112 "UD01MZ.f"
    *info = 0;

/*     Test the input scalar arguments. */

#line 116 "UD01MZ.f"
    if (*m < 1) {
#line 117 "UD01MZ.f"
	*info = -1;
#line 118 "UD01MZ.f"
    } else if (*n < 1) {
#line 119 "UD01MZ.f"
	*info = -2;
#line 120 "UD01MZ.f"
    } else if (*l < 1 || *l > 3) {
#line 121 "UD01MZ.f"
	*info = -3;
#line 122 "UD01MZ.f"
    } else if (*nout < 0) {
#line 123 "UD01MZ.f"
	*info = -4;
#line 124 "UD01MZ.f"
    } else if (*lda < *m) {
#line 125 "UD01MZ.f"
	*info = -6;
#line 126 "UD01MZ.f"
    }

#line 128 "UD01MZ.f"
    if (*info != 0) {

/*        Error return. */

#line 132 "UD01MZ.f"
	i__1 = -(*info);
#line 132 "UD01MZ.f"
	xerbla_("UD01MZ", &i__1, (ftnlen)6);
#line 133 "UD01MZ.f"
	return 0;
#line 134 "UD01MZ.f"
    }

#line 136 "UD01MZ.f"
    lentxt = i_len(text, text_len);

#line 138 "UD01MZ.f"
    for (ltext = min(72,lentxt); ltext >= 2; --ltext) {
#line 139 "UD01MZ.f"
	if (*(unsigned char *)&text[ltext - 1] != ' ') {
#line 139 "UD01MZ.f"
	    goto L40;
#line 139 "UD01MZ.f"
	}
#line 140 "UD01MZ.f"
/* L20: */
#line 140 "UD01MZ.f"
    }

#line 142 "UD01MZ.f"
L40:
#line 143 "UD01MZ.f"
    io___3.ciunit = *nout;
#line 143 "UD01MZ.f"
    s_wsfe(&io___3);
#line 143 "UD01MZ.f"
    do_fio(&c__1, text, ltext);
#line 143 "UD01MZ.f"
    do_fio(&c__1, (char *)&(*m), (ftnlen)sizeof(integer));
#line 143 "UD01MZ.f"
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
#line 143 "UD01MZ.f"
    e_wsfe();
#line 144 "UD01MZ.f"
    n1 = (*n - 1) / *l;
#line 145 "UD01MZ.f"
    j1 = 1;
#line 146 "UD01MZ.f"
    j2 = *l;

#line 148 "UD01MZ.f"
    i__1 = n1;
#line 148 "UD01MZ.f"
    for (j = 1; j <= i__1; ++j) {
#line 149 "UD01MZ.f"
	io___8.ciunit = *nout;
#line 149 "UD01MZ.f"
	s_wsfe(&io___8);
#line 149 "UD01MZ.f"
	i__2 = j2;
#line 149 "UD01MZ.f"
	for (jj = j1; jj <= i__2; ++jj) {
#line 149 "UD01MZ.f"
	    do_fio(&c__1, (char *)&jj, (ftnlen)sizeof(integer));
#line 149 "UD01MZ.f"
	}
#line 149 "UD01MZ.f"
	e_wsfe();

#line 151 "UD01MZ.f"
	i__2 = *m;
#line 151 "UD01MZ.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 152 "UD01MZ.f"
	    io___11.ciunit = *nout;
#line 152 "UD01MZ.f"
	    s_wsfe(&io___11);
#line 152 "UD01MZ.f"
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
#line 152 "UD01MZ.f"
	    i__3 = j2;
#line 152 "UD01MZ.f"
	    for (jj = j1; jj <= i__3; ++jj) {
#line 152 "UD01MZ.f"
		do_fio(&c__2, (char *)&a[i__ + jj * a_dim1], (ftnlen)sizeof(
			doublereal));
#line 152 "UD01MZ.f"
	    }
#line 152 "UD01MZ.f"
	    e_wsfe();
#line 153 "UD01MZ.f"
/* L60: */
#line 153 "UD01MZ.f"
	}

#line 155 "UD01MZ.f"
	io___12.ciunit = *nout;
#line 155 "UD01MZ.f"
	s_wsfe(&io___12);
#line 155 "UD01MZ.f"
	e_wsfe();
#line 156 "UD01MZ.f"
	j1 += *l;
#line 157 "UD01MZ.f"
	j2 += *l;
#line 158 "UD01MZ.f"
/* L80: */
#line 158 "UD01MZ.f"
    }

#line 160 "UD01MZ.f"
    io___13.ciunit = *nout;
#line 160 "UD01MZ.f"
    s_wsfe(&io___13);
#line 160 "UD01MZ.f"
    i__1 = *n;
#line 160 "UD01MZ.f"
    for (j = j1; j <= i__1; ++j) {
#line 160 "UD01MZ.f"
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
#line 160 "UD01MZ.f"
    }
#line 160 "UD01MZ.f"
    e_wsfe();

#line 162 "UD01MZ.f"
    i__1 = *m;
#line 162 "UD01MZ.f"
    for (i__ = 1; i__ <= i__1; ++i__) {
#line 163 "UD01MZ.f"
	io___14.ciunit = *nout;
#line 163 "UD01MZ.f"
	s_wsfe(&io___14);
#line 163 "UD01MZ.f"
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
#line 163 "UD01MZ.f"
	i__2 = *n;
#line 163 "UD01MZ.f"
	for (jj = j1; jj <= i__2; ++jj) {
#line 163 "UD01MZ.f"
	    do_fio(&c__2, (char *)&a[i__ + jj * a_dim1], (ftnlen)sizeof(
		    doublereal));
#line 163 "UD01MZ.f"
	}
#line 163 "UD01MZ.f"
	e_wsfe();
#line 164 "UD01MZ.f"
/* L100: */
#line 164 "UD01MZ.f"
    }

#line 166 "UD01MZ.f"
    io___15.ciunit = *nout;
#line 166 "UD01MZ.f"
    s_wsfe(&io___15);
#line 166 "UD01MZ.f"
    e_wsfe();

#line 168 "UD01MZ.f"
    return 0;

/* *** Last line of UD01MZ *** */
} /* ud01mz_ */

