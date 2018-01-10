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


class F2cpReader(object):
    def __init__(self):
        self.re_function_name = self.get_function_name_pattern()

    def __del__(self):
        del self.re_function_name

    @staticmethod
    def get_function_name_pattern():
        return re.compile(r'.*\s+(.*?)\s*\(')

    @staticmethod
    def get_calling_function_name_pattern():
        # for lines 2nd ~
        return re.compile(r'/\*:ref:\s+(.*?)\s')

    def parse_f2c_p(self, f2c_p_file_path):
        with open(f2c_p_file_path) as f:
            lines = f.readlines()
        # first line : c definitions
        # second line and after : list of other functions called

        result = self.find_c_function_name(lines[0])
        print(result)
        print(lines[0][result['start']:result['end']])

    def find_c_function_name(self, f2c_p_first_line):
        """
        From the first line of the f2c P file, find function name using regex

        :param str f2c_p_first_line: first line of f2c P file
        :return: function name string
        """
        match = self.re_function_name.search(f2c_p_first_line)
        result = {
            'name': match.groups()[0],
            'start': match.regs[1][0],
            'end': match.regs[1][1],
        }

        return result


if __name__ == '__main__':
    reader = F2cpReader()
    reader.parse_f2c_p(os.path.join(default_path_dict['slycot']['f2c'], 'AB09AD.P'))
