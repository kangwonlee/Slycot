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
        self.re_latter_line = self.get_latter_lines_pattern()
        self.re_first_line = self.get_first_line_pattern()
        self.re_arg_type_name_split = self.get_arg_type_name_split()

    def __del__(self):
        del self.re_function_name
        del self.re_latter_line
        del self.re_first_line
        del self.re_arg_type_name_split

    @staticmethod
    def get_function_name_pattern():
        return re.compile(r'\w\s+\w+\s+(\w+?)\s*\(')

    @staticmethod
    def get_first_line_pattern():
        return re.compile(r'\w\s+(?P<return_type>\w+)\s+(?P<name>\w+?)\s*\((?P<arg_list>.*)\);')

    @staticmethod
    def get_calling_function_name_pattern():
        # for lines 2nd ~
        return re.compile(r'/\*:ref:\s+(.*?)\s')

    @staticmethod
    def get_latter_lines_pattern():
        # for lines 2nd ~
        return re.compile(r'/\*:ref:\s+(?P<name>.*?)\s(?P<return_type>\d+)\s(?P<no_args>\d+)\s(?P<arg_types>.+)\s+\*/')

    @staticmethod
    def get_arg_type_name_split():
        return re.compile(r'(?P<type>(\(.+\s+\*\))|(.+\s+\*)|(\w+))\s?(?P<name>.+)')

    def parse_f2c_p(self, f2c_p_file_path):
        with open(f2c_p_file_path) as f:
            lines = f.readlines()
        # first line : c definitions
        # second line and after : list of other functions called

        result = self.find_function_info(lines[0])
        print(result)

        # functions used inside
        for latter_line in lines[1:]:
            info = self.find_calling_function_info(latter_line)
            print(info)

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

    def find_function_info(self, f2c_p_first_line):
        """
        Collect information about the function from the first line of the P file

        :param str f2c_p_first_line:
        :return: {'name': str, 'return type': str, '# arg': int, 'arg list': [str]}
        """
        match = self.re_first_line.search(f2c_p_first_line)
        arg_list = [s.strip() for s in match.group('arg_list').split(',')]

        # identify argument type and name
        arg_type_name_list = []
        for arg_type_name_str in arg_list:
            split = self.re_arg_type_name_split.search(arg_type_name_str)
            arg_type_name_list.append((split.group('type'), split.group('name')))

        result = {
            'name': match.group('name'),
            'return type': match.group('return_type'),
            '# arg': len(arg_list),
            'arg list': arg_type_name_list,
        }

        return result

    def find_calling_function_info(self, f2c_p_latter_line):
        match = self.re_latter_line.search(f2c_p_latter_line)
        result = {
            'name': match.group('name'),
            'return type': int(match.group('return_type')),
            '# arg': int(match.group('no_args')),
            'arg types': [int(s) for s in match.group('arg_types').split()],
        }

        return result


if __name__ == '__main__':
    reader = F2cpReader()
    reader.parse_f2c_p(os.path.join(default_path_dict['slycot']['f2c'], 'AB09AD.P'))
