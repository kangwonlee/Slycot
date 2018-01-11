import re


class ReadF2cHeader(object):
    def __init__(self):
        self.re_c_comment = self.get_c_comment_pattern()
        self.re_typedef_struct = self.get_c_typedef_struct_parenthesis()

    def __del__(self):
        del self.re_c_comment
        del self.re_typedef_struct

    @staticmethod
    def get_c_comment_pattern():
        # https://stackoverflow.com/questions/3534507/python-regex-matching-pattern-over-multiple-lines-why-isnt-this-working
        return re.compile(r'/\*.*?\*/', flags=re.S)

    @staticmethod
    def get_c_typedef_struct_parenthesis():
        return re.compile(r'typedef\s+struct\s+{.*?}\s*(?P<name>\w+)\s*;', flags=re.S)

    def remove_c_comment(self, txt):
        """
        From C source code text, remove /* ... */ comments

        :param str txt: C source code text
        :return:
        """
        # remove comments
        # https://stackoverflow.com/questions/17620301/how-does-the-python-re-sub-function-work
        txt_without_comments = self.re_c_comment.sub('', txt,)

        return txt_without_comments

    def replace_typedef_struct_parenthesis(self, txt):
        """
        From C source code text, replace 'typedef struct {} <type name> ;' with 'ctypedef struct <name>:\n pass'

        :param str txt: C source code text
        :return:
        """
        # first, get names
        names_list = self.re_typedef_struct.findall(txt)
        for name in names_list:
            new_text = 'ctypedef struct %s:\n    pass' % name
            # https://stackoverflow.com/questions/17620301/how-does-the-python-re-sub-function-work
            txt = self.re_typedef_struct.sub(new_text, txt, 1)

        return txt

    @staticmethod
    def get_header():
        return 'cdef extern from "f2c.h":'

    def read_f2c_h(self, f2c_h_path):
        """
        Objective : read f2c.h file, collect information, and produce a pxd file

        :param str f2c_h_path: path to the f2c.h file
        :return:
        """

        with open(f2c_h_path) as f:
            txt = f.read()

        txt_without_comments = self.remove_c_comment(txt)

        f2c_h_lines = txt_without_comments.splitlines()
