#line 1 "UE01MD.f"
/* UE01MD.f -- translated by f2c (version 20100827).
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

#line 1 "UE01MD.f"
/* Table of constant values */

static integer c__1 = 1;
static integer c_n1 = -1;
static integer c__2 = 2;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__8 = 8;

integer ue01md_(integer *ispec, char *name__, char *opts, integer *n1, 
	integer *n2, integer *n3, ftnlen name_len, ftnlen opts_len)
{
    /* System generated locals */
    integer ret_val, i__1, i__2;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    static integer i__;
    static char c1[1], c2[2], c3[1];
    static integer ic, nb, iz, nx;
    static logical cname, sname;
    static integer nbmin;
    extern integer ilaenv_(integer *, char *, char *, integer *, integer *, 
	    integer *, integer *, ftnlen, ftnlen);
    static char subnam[6];


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

/*     To provide an extension of the LAPACK routine ILAENV to */
/*     machine-specific parameters for SLICOT routines. */

/*     The default values in this version aim to give good performance on */
/*     a wide range of computers. For optimal performance, however, the */
/*     user is advised to modify this routine. Note that an optimized */
/*     BLAS is a crucial prerequisite for any speed gains. For further */
/*     details, see ILAENV. */

/*     FUNCTION VALUE */

/*     UE01MD  INTEGER */
/*             The function value set according to ISPEC. */

/*     ARGUMENTS */

/*     Input/Output Parameters */

/*     ISPEC   (input) INTEGER */
/*             Specifies the parameter to be returned as the value of */
/*             UE01MD, as follows: */
/*             = 1: the optimal blocksize; if the returned value is 1, an */
/*                  unblocked algorithm will give the best performance; */
/*             = 2: the minimum block size for which the block routine */
/*                  should be used; if the usable block size is less than */
/*                  this value, an unblocked routine should be used; */
/*             = 3: the crossover point (in a block routine, for N less */
/*                  than this value, an unblocked routine should be used) */
/*             = 4: the number of shifts, used in the product eigenvalue */
/*                  routine; */
/*             = 8: the crossover point for the multishift QR method for */
/*                  product eigenvalue problems. */

/*     NAME    (input) CHARACTER*(*) */
/*             The name of the calling subroutine, in either upper case */
/*             or lower case. */

/*     OPTS    (input) CHARACTER*(*) */
/*             The character options to the subroutine NAME, concatenated */
/*             into a single character string. */

/*     N1      (input) INTEGER */
/*     N2      (input) INTEGER */
/*     N3      (input) INTEGER */
/*             Problem dimensions for the subroutine NAME; these may not */
/*             all be required. */

/*     CONTRIBUTORS */

/*     D. Kressner, Technical Univ. Berlin, Germany, and */
/*     P. Benner, Technical Univ. Chemnitz, Germany, December 2003. */

/*     REVISIONS */

/*     V. Sima, June 2008 (SLICOT version of the HAPACK routine ILAHAP). */

/*     ****************************************************************** */

/*     .. Scalar Arguments .. */

/*     .. Local Scalars .. */
/*     .. External Functions .. */
/*     .. Intrinsic Functions .. */

/*     .. Executable Statements .. */

#line 99 "UE01MD.f"
    if (*ispec == 1 || *ispec == 2 || *ispec == 3) {

/*        Convert NAME to upper case if the first character is lower */
/*        case. */

#line 104 "UE01MD.f"
	ret_val = 1;
#line 105 "UE01MD.f"
	s_copy(subnam, name__, (ftnlen)6, name_len);
#line 106 "UE01MD.f"
	ic = *(unsigned char *)subnam;
#line 107 "UE01MD.f"
	iz = 'Z';
#line 108 "UE01MD.f"
	if (iz == 90 || iz == 122) {

/*           ASCII character set. */

#line 112 "UE01MD.f"
	    if (ic >= 97 && ic <= 122) {
#line 113 "UE01MD.f"
		*(unsigned char *)subnam = (char) (ic - 32);
#line 114 "UE01MD.f"
		for (i__ = 2; i__ <= 6; ++i__) {
#line 115 "UE01MD.f"
		    ic = *(unsigned char *)&subnam[i__ - 1];
#line 116 "UE01MD.f"
		    if (ic >= 97 && ic <= 122) {
#line 116 "UE01MD.f"
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 116 "UE01MD.f"
		    }
#line 118 "UE01MD.f"
/* L10: */
#line 118 "UE01MD.f"
		}
#line 119 "UE01MD.f"
	    }

#line 121 "UE01MD.f"
	} else if (iz == 233 || iz == 169) {

/*           EBCDIC character set. */

#line 125 "UE01MD.f"
	    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || ic >= 162 
		    && ic <= 169) {
#line 128 "UE01MD.f"
		*(unsigned char *)subnam = (char) (ic + 64);
#line 129 "UE01MD.f"
		for (i__ = 2; i__ <= 6; ++i__) {
#line 130 "UE01MD.f"
		    ic = *(unsigned char *)&subnam[i__ - 1];
#line 131 "UE01MD.f"
		    if (ic >= 129 && ic <= 137 || ic >= 145 && ic <= 153 || 
			    ic >= 162 && ic <= 169) {
#line 131 "UE01MD.f"
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic + 64);
#line 131 "UE01MD.f"
		    }
#line 135 "UE01MD.f"
/* L20: */
#line 135 "UE01MD.f"
		}
#line 136 "UE01MD.f"
	    }

#line 138 "UE01MD.f"
	} else if (iz == 218 || iz == 250) {

/*           Prime machines:  ASCII+128. */

#line 142 "UE01MD.f"
	    if (ic >= 225 && ic <= 250) {
#line 143 "UE01MD.f"
		*(unsigned char *)subnam = (char) (ic - 32);
#line 144 "UE01MD.f"
		for (i__ = 2; i__ <= 6; ++i__) {
#line 145 "UE01MD.f"
		    ic = *(unsigned char *)&subnam[i__ - 1];
#line 146 "UE01MD.f"
		    if (ic >= 225 && ic <= 250) {
#line 146 "UE01MD.f"
			*(unsigned char *)&subnam[i__ - 1] = (char) (ic - 32);
#line 146 "UE01MD.f"
		    }
#line 148 "UE01MD.f"
/* L30: */
#line 148 "UE01MD.f"
		}
#line 149 "UE01MD.f"
	    }
#line 150 "UE01MD.f"
	}

#line 152 "UE01MD.f"
	*(unsigned char *)c1 = *(unsigned char *)subnam;
#line 153 "UE01MD.f"
	sname = *(unsigned char *)c1 == 'S' || *(unsigned char *)c1 == 'D';
#line 154 "UE01MD.f"
	cname = *(unsigned char *)c1 == 'C' || *(unsigned char *)c1 == 'Z';
#line 155 "UE01MD.f"
	if (! (cname || sname)) {
#line 155 "UE01MD.f"
	    return ret_val;
#line 155 "UE01MD.f"
	}
#line 157 "UE01MD.f"
	s_copy(c2, subnam + 3, (ftnlen)2, (ftnlen)2);
#line 158 "UE01MD.f"
	*(unsigned char *)c3 = *(unsigned char *)&subnam[5];

#line 160 "UE01MD.f"
	if (*ispec == 1) {

/*           Block size. */

#line 164 "UE01MD.f"
	    nb = 1;
#line 165 "UE01MD.f"
	    if (s_cmp(c2, "4S", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "4T", 
		    (ftnlen)2, (ftnlen)2) == 0) {
#line 166 "UE01MD.f"
		if (*(unsigned char *)c3 == 'B') {
#line 167 "UE01MD.f"
		    nb = ilaenv_(&c__1, "DGEQRF", " ", n1, n2, &c_n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
#line 168 "UE01MD.f"
		} else if (*(unsigned char *)c3 == 'T') {
#line 169 "UE01MD.f"
		    nb = ilaenv_(&c__1, "DGEHRD", " ", n1, n2, n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 4;
#line 170 "UE01MD.f"
		}
#line 171 "UE01MD.f"
	    } else if (s_cmp(c2, "4P", (ftnlen)2, (ftnlen)2) == 0) {
#line 172 "UE01MD.f"
		if (*(unsigned char *)c3 == 'B') {
#line 173 "UE01MD.f"
		    nb = ilaenv_(&c__1, "DGEHRD", " ", n1, n2, n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
#line 174 "UE01MD.f"
		}
#line 175 "UE01MD.f"
	    } else if (s_cmp(c2, "4W", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2,
		     "4Q", (ftnlen)2, (ftnlen)2) == 0) {
#line 176 "UE01MD.f"
		if (*(unsigned char *)c3 == 'D') {
#line 177 "UE01MD.f"
		    nb = ilaenv_(&c__1, "DORGQR", " ", n1, n2, n3, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
#line 178 "UE01MD.f"
		} else if (*(unsigned char *)c3 == 'B') {
#line 179 "UE01MD.f"
		    nb = ilaenv_(&c__1, "DORMQR", " ", n1, n2, n3, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
#line 180 "UE01MD.f"
		}
/* *          ELSE IF ( C2.EQ.'SH' ) THEN */
/* *             IF ( C3.EQ.'PVB' ) THEN */
/* *                NB = ILAENV( 1, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 2 */
/* *             END IF */
#line 185 "UE01MD.f"
	    }
#line 186 "UE01MD.f"
	    ret_val = nb;
#line 187 "UE01MD.f"
	} else if (*ispec == 2) {

/*           Minimum block size. */

#line 191 "UE01MD.f"
	    nbmin = 2;
#line 192 "UE01MD.f"
	    if (s_cmp(c2, "4S", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "4T", 
		    (ftnlen)2, (ftnlen)2) == 0) {
#line 193 "UE01MD.f"
		if (*(unsigned char *)c3 == 'B') {
/* Computing MAX */
#line 194 "UE01MD.f"
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DGEQRF", " ", n1, n2, &
			    c_n1, &c_n1, (ftnlen)6, (ftnlen)1) / 2;
#line 194 "UE01MD.f"
		    nbmin = max(i__1,i__2);
#line 196 "UE01MD.f"
		} else if (*(unsigned char *)c3 == 'T') {
/* Computing MAX */
#line 197 "UE01MD.f"
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DGEHRD", " ", n1, n2, n1,
			     &c_n1, (ftnlen)6, (ftnlen)1) / 4;
#line 197 "UE01MD.f"
		    nbmin = max(i__1,i__2);
#line 199 "UE01MD.f"
		}
#line 200 "UE01MD.f"
	    } else if (s_cmp(c2, "4P", (ftnlen)2, (ftnlen)2) == 0) {
#line 201 "UE01MD.f"
		if (*(unsigned char *)c3 == 'B') {
/* Computing MAX */
#line 202 "UE01MD.f"
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DGEHRD", " ", n1, n2, n1,
			     &c_n1, (ftnlen)6, (ftnlen)1) / 4;
#line 202 "UE01MD.f"
		    nbmin = max(i__1,i__2);
#line 204 "UE01MD.f"
		}
#line 205 "UE01MD.f"
	    } else if (s_cmp(c2, "4W", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2,
		     "4Q", (ftnlen)2, (ftnlen)2) == 0) {
#line 206 "UE01MD.f"
		if (*(unsigned char *)c3 == 'D') {
/* Computing MAX */
#line 207 "UE01MD.f"
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQR", " ", n1, n2, n3,
			     &c_n1, (ftnlen)6, (ftnlen)1) / 2;
#line 207 "UE01MD.f"
		    nbmin = max(i__1,i__2);
#line 209 "UE01MD.f"
		} else if (*(unsigned char *)c3 == 'B') {
/* Computing MAX */
#line 210 "UE01MD.f"
		    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMQR", " ", n1, n2, n3,
			     &c_n1, (ftnlen)6, (ftnlen)1) / 2;
#line 210 "UE01MD.f"
		    nbmin = max(i__1,i__2);
#line 212 "UE01MD.f"
		}
/* *          ELSE IF ( C2.EQ.'SH' ) THEN */
/* *             IF ( C3.EQ.'PVB' ) THEN */
/* *                NBMIN = MAX( 2, ILAENV( 2, 'DGEHRD', ' ', N1, N2, N1, */
/* *   $                                    -1 ) / 4 ) */
/* *             END IF */
#line 218 "UE01MD.f"
	    }
#line 219 "UE01MD.f"
	    ret_val = nbmin;
#line 220 "UE01MD.f"
	} else if (*ispec == 3) {

/*           Crossover point. */

#line 224 "UE01MD.f"
	    nx = 0;
#line 225 "UE01MD.f"
	    if (s_cmp(c2, "4S", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2, "4T", 
		    (ftnlen)2, (ftnlen)2) == 0) {
#line 226 "UE01MD.f"
		if (*(unsigned char *)c3 == 'B') {
#line 227 "UE01MD.f"
		    nx = ilaenv_(&c__3, "DGEQRF", " ", n1, n2, &c_n1, &c_n1, (
			    ftnlen)6, (ftnlen)1);
#line 228 "UE01MD.f"
		} else if (*(unsigned char *)c3 == 'T') {
#line 229 "UE01MD.f"
		    nx = ilaenv_(&c__3, "DGEHRD", " ", n1, n2, n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
#line 230 "UE01MD.f"
		}
#line 231 "UE01MD.f"
	    } else if (s_cmp(c2, "4P", (ftnlen)2, (ftnlen)2) == 0) {
#line 232 "UE01MD.f"
		if (*(unsigned char *)c3 == 'B') {
#line 233 "UE01MD.f"
		    nx = ilaenv_(&c__3, "DGEHRD", " ", n1, n2, n1, &c_n1, (
			    ftnlen)6, (ftnlen)1) / 2;
#line 234 "UE01MD.f"
		}
#line 235 "UE01MD.f"
	    } else if (s_cmp(c2, "4W", (ftnlen)2, (ftnlen)2) == 0 || s_cmp(c2,
		     "4Q", (ftnlen)2, (ftnlen)2) == 0) {
#line 236 "UE01MD.f"
		if (*(unsigned char *)c3 == 'D') {
#line 237 "UE01MD.f"
		    nx = ilaenv_(&c__3, "DORGQR", " ", n1, n2, n3, &c_n1, (
			    ftnlen)6, (ftnlen)1);
#line 238 "UE01MD.f"
		} else if (*(unsigned char *)c3 == 'B') {
#line 239 "UE01MD.f"
		    nx = ilaenv_(&c__3, "DORGQR", " ", n1, n2, n3, &c_n1, (
			    ftnlen)6, (ftnlen)1);
#line 240 "UE01MD.f"
		}
/* *          ELSE IF ( C2.EQ.'SH' ) THEN */
/* *             IF ( C3.EQ.'PVB' ) THEN */
/* *                NX = ILAENV( 3, 'DGEHRD', ' ', N1, N2, N1, -1 ) / 2 */
/* *             END IF */
#line 245 "UE01MD.f"
	    }
#line 246 "UE01MD.f"
	    ret_val = nx;
#line 247 "UE01MD.f"
	}
#line 248 "UE01MD.f"
    } else if (*ispec == 4) {

/*        Number of shifts (used by MB03XP). */

#line 252 "UE01MD.f"
	ret_val = ilaenv_(&c__4, "DHSEQR", opts, n1, n2, n3, &c_n1, (ftnlen)6,
		 opts_len);
#line 253 "UE01MD.f"
    } else if (*ispec == 8) {

/*        Crossover point for multishift (used by MB03XP). */

#line 257 "UE01MD.f"
	ret_val = ilaenv_(&c__8, "DHSEQR", opts, n1, n2, n3, &c_n1, (ftnlen)6,
		 opts_len);
#line 258 "UE01MD.f"
    } else {

/*        Invalid value for ISPEC. */

#line 262 "UE01MD.f"
	ret_val = -1;
#line 263 "UE01MD.f"
    }
#line 264 "UE01MD.f"
    return ret_val;
/* *** Last line of UE01MD *** */
} /* ue01md_ */

