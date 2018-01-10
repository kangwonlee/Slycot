import os
import re

default_path_dict = {
    'slycot': {
        'src': os.path.join(os.pardir, 'slycot', 'src'),
        'f2c': os.path.join(os.pardir, 'slycot', 'src-f2c'),
    },
    'lapack': {
        'src': os.path.join(os.pardir, 'slycot', 'src-f2c', 'lapack', 'SRC'),
        'f2c': os.path.join(os.pardir, 'slycot', 'src-f2c', 'lapack', 'src-f2c'),
    },
    'blas': {
        'src': os.path.join(os.pardir, 'slycot', 'src-f2c', 'lapack', 'BLAS', 'SRC'),
        'f2c': os.path.join(os.pardir, 'slycot', 'src-f2c', 'lapack', 'BLAS', 'src-f2c'),
    },
}


def get_function_name_pattern():
    return re.compile(r'.*\s+(.*?)\s*\(')


def parse_f2c_p(f2c_p_file_path):
    with open(f2c_p_file_path) as f:
        lines = f.readlines()
    # first line : c definitions
    # second line and after : list of other functions called

    r = get_function_name_pattern()

    result = r.search(lines[0])
    print(result.groups())


if __name__ == '__main__':
    parse_f2c_p(os.path.join(default_path_dict['slycot']['f2c'], 'AB09AD.P'))
