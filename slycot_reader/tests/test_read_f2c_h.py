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
