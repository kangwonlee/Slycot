import re


class ReadF2cHeader(object):
    def __init__(self):
        self.re_c_comment = self.get_c_comment_pattern()

    @staticmethod
    def get_c_comment_pattern():
        # https://stackoverflow.com/questions/3534507/python-regex-matching-pattern-over-multiple-lines-why-isnt-this-working
        return re.compile(r'/\*.*?\*/', flags=re.S)

    @staticmethod
    def get_c_typedef_struct_parenthesis():
        return re.compile(r'typedef\s+struct\s+\{.*?\}\s*(?P<name>\w+)\s*;', flags=re.S)

    def remove_c_comment(self, txt):
        """
        From C source code text, remove /* ... */ comments

        :param str txt: C source code text
        :return:
        """
        # remove comments
        txt_without_comments = self.re_c_comment.sub(txt, '')

        return txt_without_comments

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
