import re


def get_c_comment_pattern():
    # https://stackoverflow.com/questions/3534507/python-regex-matching-pattern-over-multiple-lines-why-isnt-this-working
    return re.compile(r'/\*.*?\*/', flags=re.S)


def read_f2c_h(f2c_h_path):
    """
    Objective : read f2c.h file, collect information, and produce a pxd file

    :param str f2c_h_path: path to the f2c.h file
    :return:
    """

    with open(f2c_h_path) as f:
        txt = f.read()

    # remove comments
