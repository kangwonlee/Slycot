#line 1 "UD01ND.f"
/* UD01ND.f -- translated by f2c (version 20100827).
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

#line 1 "UD01ND.f"
/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int ud01nd_(integer *mp, integer *np, integer *dp, integer *
	l, integer *nout, doublereal *p, integer *ldp1, integer *ldp2, char *
	text, integer *info, ftnlen text_len)
{
    /* Format strings */
    static char fmt_99999[] = "(\002 \002)";
    static char fmt_99998[] = "(/,1x,a,\002(\002,i2,\002)\002,\002 (\002,i2"
	    ",\002X\002,i2,\002)\002)";
    static char fmt_99997[] = "(5x,5(6x,i2,7x))";
    static char fmt_99996[] = "(1x,i2,2x,5d15.7)";

    /* System generated locals */
    integer p_dim1, p_dim2, p_offset, i__1, i__2, i__3, i__4;

    /* Builtin functions */
    integer i_len(char *, ftnlen), s_wsfe(cilist *), e_wsfe(void), do_fio(
	    integer *, char *, ftnlen);

    /* Local variables */
    static integer i__, j, k, j1, j2, n1, jj, ltext;
    extern /* Subroutine */ int xerbla_(char *, integer *, ftnlen);
    static integer lentxt;

    /* Fortran I/O blocks */
    static cilist io___4 = { 0, 0, 0, fmt_99999, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_99998, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_99997, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_99996, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_99997, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_99996, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_99999, 0 };



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

/*     To print the MP-by-NP coefficient matrices of a matrix polynomial */
/*                                                    dp-1           dp */
/*        P(s) = P(0) + P(1) * s + . . . + P(dp-1) * s    + P(dp) * s  . */

/*     The elements of the matrices are output to 7 significant figures. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     MP      (input) INTEGER */
/*             The number of rows of the matrix polynomial P(s). */
/*             MP >= 1. */

/*     NP      (input) INTEGER */
/*             The number of columns of the matrix polynomial P(s). */
/*             NP >= 1. */

/*     DP      (input) INTEGER */
/*             The degree of the matrix polynomial P(s).  DP >= 0. */

/*     L       (input) INTEGER */
/*             The number of elements of the coefficient matrices to be */
/*             printed per line.  1 <= L <= 5. */

/*     NOUT    (input) INTEGER */
/*             The output channel to which the results are sent. */
/*             NOUT >= 0. */

/*     P       (input) DOUBLE PRECISION array, dimension (LDP1,LDP2,DP+1) */
/*             The leading MP-by-NP-by-(DP+1) part of this array must */
/*             contain the coefficients of the matrix polynomial P(s). */
/*             Specifically, P(i,j,k) must contain the coefficient of */
/*             s**(k-1) of the polynomial which is the (i,j)-th element */
/*             of P(s), where i = 1,2,...,MP, j = 1,2,...,NP and */
/*             k = 1,2,...,DP+1. */

/*     LDP1    INTEGER */
/*             The leading dimension of array P.  LDP1 >= MP. */

/*     LDP2    INTEGER */
/*             The second dimension of array P.  LDP2 >= NP. */

/*     TEXT    (input) CHARACTER*72 */
/*             Title caption of the coefficient matrices to be printed. */
/*             TEXT is followed by the degree of the coefficient matrix, */
/*             within brackets. If TEXT = ' ', then the coefficient */
/*             matrices are separated by an empty line. */

/*     Error Indicator */

/*     INFO    INTEGER */
/*             = 0:  successful exit; */
/*             < 0:  if INFO = -i, the i-th argument had an illegal */
/*                   value. */

/*     METHOD */

/*     For i = 1, 2, ..., DP + 1 the routine first prints the contents of */
/*     TEXT followed by (i-1) as a title, followed by the elements of the */
/*     MP-by-NP coefficient matrix P(i) such that */
/*     (i)  if NP < L, then the leading MP-by-NP part is printed; */
/*     (ii) if NP = k*L + p (where k, p > 0), then k MP-by-L blocks of */
/*          consecutive columns of P(i) are printed one after another */
/*          followed by one MP-by-p block containing the last p columns */
/*          of P(i). */
/*     Row numbers are printed on the left of each row and a column */
/*     number on top of each column. */

/*     REFERENCES */

/*     None. */

/*     NUMERICAL ASPECTS */

/*     None. */

/*     CONTRIBUTOR */

/*     Release 3.0: V. Sima, Katholieke Univ. Leuven, Belgium, June 1998. */
/*     Based on routine PRMAPO by A.J. Geurts, Eindhoven University of */
/*     Technology, Holland. */

/*     REVISIONS */

/*     - */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */
/*     .. Array Arguments .. */
/*     .. Local Scalars .. */
/*     .. External Subroutines .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

#line 127 "UD01ND.f"
    /* Parameter adjustments */
#line 127 "UD01ND.f"
    p_dim1 = *ldp1;
#line 127 "UD01ND.f"
    p_dim2 = *ldp2;
#line 127 "UD01ND.f"
    p_offset = 1 + p_dim1 * (1 + p_dim2);
#line 127 "UD01ND.f"
    p -= p_offset;
#line 127 "UD01ND.f"

#line 127 "UD01ND.f"
    /* Function Body */
#line 127 "UD01ND.f"
    *info = 0;

/*     Check the input scalar arguments. */

#line 131 "UD01ND.f"
    if (*mp < 1) {
#line 132 "UD01ND.f"
	*info = -1;
#line 133 "UD01ND.f"
    } else if (*np < 1) {
#line 134 "UD01ND.f"
	*info = -2;
#line 135 "UD01ND.f"
    } else if (*dp < 0) {
#line 136 "UD01ND.f"
	*info = -3;
#line 137 "UD01ND.f"
    } else if (*l < 1 || *l > 5) {
#line 138 "UD01ND.f"
	*info = -4;
#line 139 "UD01ND.f"
    } else if (*nout < 0) {
#line 140 "UD01ND.f"
	*info = -5;
#line 141 "UD01ND.f"
    } else if (*ldp1 < *mp) {
#line 142 "UD01ND.f"
	*info = -7;
#line 143 "UD01ND.f"
    } else if (*ldp2 < *np) {
#line 144 "UD01ND.f"
	*info = -8;
#line 145 "UD01ND.f"
    }

#line 147 "UD01ND.f"
    if (*info != 0) {

/*        Error return. */

#line 151 "UD01ND.f"
	i__1 = -(*info);
#line 151 "UD01ND.f"
	xerbla_("UD01ND", &i__1, (ftnlen)6);
#line 152 "UD01ND.f"
	return 0;
#line 153 "UD01ND.f"
    }

#line 155 "UD01ND.f"
    lentxt = i_len(text, text_len);
#line 156 "UD01ND.f"
    ltext = min(72,lentxt);
/*     WHILE ( TEXT(LTEXT:LTEXT) =  ' ' ) DO */
#line 158 "UD01ND.f"
L10:
#line 158 "UD01ND.f"
    if (*(unsigned char *)&text[ltext - 1] == ' ') {
#line 159 "UD01ND.f"
	--ltext;
#line 160 "UD01ND.f"
	goto L10;
#line 161 "UD01ND.f"
    }
/*     END WHILE 10 */

#line 164 "UD01ND.f"
    i__1 = *dp + 1;
#line 164 "UD01ND.f"
    for (k = 1; k <= i__1; ++k) {
#line 165 "UD01ND.f"
	if (ltext == 0) {
#line 166 "UD01ND.f"
	    io___4.ciunit = *nout;
#line 166 "UD01ND.f"
	    s_wsfe(&io___4);
#line 166 "UD01ND.f"
	    e_wsfe();
#line 167 "UD01ND.f"
	} else {
#line 168 "UD01ND.f"
	    io___5.ciunit = *nout;
#line 168 "UD01ND.f"
	    s_wsfe(&io___5);
#line 168 "UD01ND.f"
	    do_fio(&c__1, text, ltext);
#line 168 "UD01ND.f"
	    i__2 = k - 1;
#line 168 "UD01ND.f"
	    do_fio(&c__1, (char *)&i__2, (ftnlen)sizeof(integer));
#line 168 "UD01ND.f"
	    do_fio(&c__1, (char *)&(*mp), (ftnlen)sizeof(integer));
#line 168 "UD01ND.f"
	    do_fio(&c__1, (char *)&(*np), (ftnlen)sizeof(integer));
#line 168 "UD01ND.f"
	    e_wsfe();
#line 169 "UD01ND.f"
	}
#line 170 "UD01ND.f"
	n1 = (*np - 1) / *l;
#line 171 "UD01ND.f"
	j1 = 1;
#line 172 "UD01ND.f"
	j2 = *l;

#line 174 "UD01ND.f"
	i__2 = n1;
#line 174 "UD01ND.f"
	for (j = 1; j <= i__2; ++j) {
#line 175 "UD01ND.f"
	    io___10.ciunit = *nout;
#line 175 "UD01ND.f"
	    s_wsfe(&io___10);
#line 175 "UD01ND.f"
	    i__3 = j2;
#line 175 "UD01ND.f"
	    for (jj = j1; jj <= i__3; ++jj) {
#line 175 "UD01ND.f"
		do_fio(&c__1, (char *)&jj, (ftnlen)sizeof(integer));
#line 175 "UD01ND.f"
	    }
#line 175 "UD01ND.f"
	    e_wsfe();

#line 177 "UD01ND.f"
	    i__3 = *mp;
#line 177 "UD01ND.f"
	    for (i__ = 1; i__ <= i__3; ++i__) {
#line 178 "UD01ND.f"
		io___13.ciunit = *nout;
#line 178 "UD01ND.f"
		s_wsfe(&io___13);
#line 178 "UD01ND.f"
		do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
#line 178 "UD01ND.f"
		i__4 = j2;
#line 178 "UD01ND.f"
		for (jj = j1; jj <= i__4; ++jj) {
#line 178 "UD01ND.f"
		    do_fio(&c__1, (char *)&p[i__ + (jj + k * p_dim2) * p_dim1]
			    , (ftnlen)sizeof(doublereal));
#line 178 "UD01ND.f"
		}
#line 178 "UD01ND.f"
		e_wsfe();
#line 179 "UD01ND.f"
/* L20: */
#line 179 "UD01ND.f"
	    }

#line 181 "UD01ND.f"
	    j1 += *l;
#line 182 "UD01ND.f"
	    j2 += *l;
#line 183 "UD01ND.f"
/* L30: */
#line 183 "UD01ND.f"
	}

#line 185 "UD01ND.f"
	io___14.ciunit = *nout;
#line 185 "UD01ND.f"
	s_wsfe(&io___14);
#line 185 "UD01ND.f"
	i__2 = *np;
#line 185 "UD01ND.f"
	for (j = j1; j <= i__2; ++j) {
#line 185 "UD01ND.f"
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
#line 185 "UD01ND.f"
	}
#line 185 "UD01ND.f"
	e_wsfe();

#line 187 "UD01ND.f"
	i__2 = *mp;
#line 187 "UD01ND.f"
	for (i__ = 1; i__ <= i__2; ++i__) {
#line 188 "UD01ND.f"
	    io___15.ciunit = *nout;
#line 188 "UD01ND.f"
	    s_wsfe(&io___15);
#line 188 "UD01ND.f"
	    do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
#line 188 "UD01ND.f"
	    i__3 = *np;
#line 188 "UD01ND.f"
	    for (jj = j1; jj <= i__3; ++jj) {
#line 188 "UD01ND.f"
		do_fio(&c__1, (char *)&p[i__ + (jj + k * p_dim2) * p_dim1], (
			ftnlen)sizeof(doublereal));
#line 188 "UD01ND.f"
	    }
#line 188 "UD01ND.f"
	    e_wsfe();
#line 189 "UD01ND.f"
/* L40: */
#line 189 "UD01ND.f"
	}

#line 191 "UD01ND.f"
/* L50: */
#line 191 "UD01ND.f"
    }

#line 193 "UD01ND.f"
    io___16.ciunit = *nout;
#line 193 "UD01ND.f"
    s_wsfe(&io___16);
#line 193 "UD01ND.f"
    e_wsfe();

#line 195 "UD01ND.f"
    return 0;


/* *** Last line of UD01ND *** */
} /* ud01nd_ */

