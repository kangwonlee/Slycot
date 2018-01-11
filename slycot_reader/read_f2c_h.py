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

    txt_without_comments = remove_c_comment(txt)


def remove_c_comment(txt):
    """
    From C source code text, remove /* ... */ comments

    :param str txt: C source code text
    :return:
    """
    re_comments = get_c_comment_pattern()
    # remove comments
    txt_without_comments = re_comments.sub(txt, '')

    return txt_without_comments
