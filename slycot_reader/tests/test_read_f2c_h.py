import os
import unittest

import slycot_reader.read_f2c_h as ch


def read_f2c_h_text():
    path_to_f2c_h = os.path.abspath(os.path.join(os.pardir, os.pardir, 'slycot', 'src-f2c', 'f2c.h'))
    assert os.path.exists(path_to_f2c_h)
    assert os.path.isfile(path_to_f2c_h)
    with open(path_to_f2c_h, 'r') as f:
        input_txt = f.read()
    return input_txt


class TestF2cH(unittest.TestCase):
    input_txt = read_f2c_h_text()

    def setUp(self):
        self.reader = ch.ReadF2cHeader()

    def test_remove_c_comment(self):
        result = self.reader.remove_c_comment(self.input_txt)

        self.assertFalse('/*' in result)
        self.assertFalse('*/' in result)

    def test_get_c_comment_pattern(self):
        r = self.reader.get_c_comment_pattern()

        result_list = r.findall(self.input_txt)

        expected_list = ['/* f2c.h  --  Standard Fortran to C header file */',
                         '''/**  barf  [ba:rf]  2.  "He suggested using FORTRAN, and everybody barfed."

	- From The Shogakukan DICTIONARY OF NEW ENGLISH (Second edition) */''',
                         '/* Adjust for integer*8. */',
                         '/* system-dependent */',
                         '/* system-dependent */',
                         '/* Extern is for use with -E */',
                         '/* I/O stuff */',
                         '/* for -i2 */',
                         '/*external read, write*/',
                         '/*internal read, write*/',
                         '/*open*/',
                         '/*close*/',
                         '/*rewind, backspace, endfile*/',
                         '/* inquire */',
                         "/*parameters in standard's order*/",
                         '/* for multiple entry points */',
                         '/* longint j; */',
                         '/*typedef long int Long;*/',
                         '/* No longer used; formerly in Namelist */',
                         '/* for Namelist */',
                         '/* procedure parameter types for -A and -C++ */',
                         '/* Unknown procedure type */',
                         '/* Complex */',
                         '/* Double Complex */',
                         '/* Character */',
                         '/* Subroutine */',
                         '/* Unknown procedure type */',
                         '/* Complex */',
                         '/* Double Complex */',
                         '/* Character */',
                         '/* Subroutine */',
                         '/* E_fp is for real functions when -R is not specified */',
                         '/* complex function */',
                         '/* character function */',
                         '/* double complex function */',
                         '/* real function with -R not specified */',
                         '/* undef any lower-case symbols that your C compiler predefines, e.g.: */',
                         ]
        self.assertSequenceEqual(expected_list, result_list)

    def test_get_c_typedef_struct_parenthesis(self):
        r = self.reader.get_c_typedef_struct_parenthesis()

        result_list = r.findall(self.input_txt)

        expected_list = ['complex', 'doublecomplex', 'cilist', 'icilist', 'olist', 'cllist', 'alist', 'inlist', ]

        self.assertSequenceEqual(expected_list, result_list)

    def test_remove_comments(self):
        input_txt = '''#ifdef INTEGER_STAR_8	/* Adjust for integer*8. */
typedef long long longint;		/* system-dependent */
typedef unsigned long long ulongint;	/* system-dependent */
#define qbit_clear(a,b)	((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)	((a) |  ((ulongint)1 << (b)))
#endif
'''
        result = self.reader.remove_c_comment(input_txt)

        expected = '''#ifdef INTEGER_STAR_8	
typedef long long longint;		
typedef unsigned long long ulongint;	
#define qbit_clear(a,b)	((a) & ~((ulongint)1 << (b)))
#define qbit_set(a,b)	((a) |  ((ulongint)1 << (b)))
#endif
'''

        self.assertEqual(expected, result)

    def test_replace_typedef_struct(self):
        input_txt_list = [
            '''typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef long int logical;
''',
            '''typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;
#endif

/*external read, write*/
typedef struct
{	flag cierr;
	ftnint ciunit;
	flag ciend;
	char *cifmt;
	ftnint cirec;
} cilist;

/*internal read, write*/
typedef struct
{	flag icierr;
	char *iciunit;
	flag iciend;
	char *icifmt;
	ftnint icirlen;
	ftnint icirnum;
} icilist;

/*open*/
typedef struct
{	flag oerr;
	ftnint ounit;
	char *ofnm;
	ftnlen ofnmlen;
	char *osta;
	char *oacc;
	char *ofm;
	ftnint orl;
	char *oblnk;
} olist;

/*close*/
typedef struct
{	flag cerr;
	ftnint cunit;
	char *csta;
} cllist;

/*rewind, backspace, endfile*/
typedef struct
{	flag aerr;
	ftnint aunit;
} alist;

/* inquire */
typedef struct
{	flag inerr;
	ftnint inunit;
	char *infile;
	ftnlen infilen;
	ftnint	*inex;	/*parameters in standard's order*/
	ftnint	*inopen;
	ftnint	*innum;
	ftnint	*innamed;
	char	*inname;
	ftnlen	innamlen;
	char	*inacc;
	ftnlen	inacclen;
	char	*inseq;
	ftnlen	inseqlen;
	char 	*indir;
	ftnlen	indirlen;
	char	*infmt;
	ftnlen	infmtlen;
	char	*inform;
	ftnint	informlen;
	char	*inunf;
	ftnlen	inunflen;
	ftnint	*inrecl;
	ftnint	*innrec;
	char	*inblank;
	ftnlen	inblanklen;
} inlist;

#define VOID void

union Multitype {	/* for multiple entry points */
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};
''',
        ]

        expected_list = [
            '''typedef double doublereal;
ctypedef struct complex:
    pass
ctypedef struct doublecomplex:
    pass
typedef long int logical;
''',
            '''typedef long int flag;
typedef long int ftnlen;
typedef long int ftnint;
#endif

/*external read, write*/
ctypedef struct cilist:
    pass

/*internal read, write*/
ctypedef struct icilist:
    pass

/*open*/
ctypedef struct olist:
    pass

/*close*/
ctypedef struct cllist:
    pass

/*rewind, backspace, endfile*/
ctypedef struct alist:
    pass

/* inquire */
ctypedef struct inlist:
    pass

#define VOID void

union Multitype {	/* for multiple entry points */
	integer1 g;
	shortint h;
	integer i;
	/* longint j; */
	real r;
	doublereal d;
	complex c;
	doublecomplex z;
	};
''',
        ]

        for input_txt, expected in zip(input_txt_list, expected_list):
            result = self.reader.replace_typedef_struct(input_txt)

            self.assertEqual(expected, result)
