import os
import pprint
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

f2c_path_dict = {
    'src': {},
    'f2c': {},
}

for lib, path_dict in default_path_dict.items():
    f2c_path_dict['src'][lib] = path_dict['src']
    f2c_path_dict['f2c'][lib] = path_dict['f2c']


class F2cpReader(object):
    def __init__(self):
        self.re_function_name = self.get_function_name_pattern()
        self.re_latter_line = self.get_latter_lines_pattern()
        self.re_first_line = self.get_first_line_pattern()
        self.re_arg_type_name_split = self.get_arg_type_name_split()
        self.big_table = {}
        self.arg_type_lookup = {}

    def __del__(self):
        del self.re_function_name
        del self.re_latter_line
        del self.re_first_line
        del self.re_arg_type_name_split
        del self.big_table
        del self.arg_type_lookup

    @staticmethod
    def get_function_name_pattern():
        return re.compile(r'\w\s+\w+\s+(\w+?)\s*\(')

    @staticmethod
    def get_first_line_pattern():
        return re.compile(r'\w\s+(?P<return_type>\w+)\s+(?P<name>\w+?)\s*\((?P<arg_list>.*)\);')

    @staticmethod
    def get_calling_function_name_pattern():
        # for lines 2nd ~
        return re.compile(r'/\*:ref:\s+(?P<name>.*?)\s')

    @staticmethod
    def get_latter_lines_pattern():
        # for lines 2nd ~
        return re.compile(r'/\*:ref:\s+(?P<name>.*?)\s(?P<return_type>\d+)\s(?P<no_args>\d+)\s(?P<arg_types>.+)\s+\*/')

    @staticmethod
    def get_arg_type_name_split():
        # '(char *)a' -> ('(char *)', 'a')
        # 'char *a' -> ('char *', 'a')
        # 'char a' -> ('char', 'a')
        return re.compile(r'(?P<type>(\(.+\s+\*\))|(.+\s+\*)|(\w+))\s?(?P<name>.+)')

    def parse_f2c_p(self, f2c_p_file_path):
        with open(f2c_p_file_path) as f:
            lines = f.readlines()
        # first line : c definitions
        # second line and after : list of other functions called

        for line in lines:
            line = line.strip()

            if not line.startswith('/*'):
                # functions defined
                info = self.find_function_info(line)
                print(info)
            else:
                # functions used inside
                info = self.find_calling_function_info(line)
            self.update_big_table(info)

    def update_big_table(self, info_dict):
        # begin update big table using the info_dict
        big_table_entry = self.big_table.get(info_dict['name'], {})
        big_table_entry.update(info_dict)
        del big_table_entry['name']
        self.big_table[info_dict['name']] = big_table_entry
        # end update big table using the first line info

        self.update_arg_type_lookup(big_table_entry)

    def update_arg_type_lookup(self, big_table_entry):
        if ('arg types' in big_table_entry) and ('arg list' in big_table_entry):
            # begin update arg_type_lookup
            for type_id, arg_type_name in zip(big_table_entry['arg types'], big_table_entry['arg list']):
                arg_type_lookup_entry = self.arg_type_lookup.get(type_id, {})
                # count number of each case
                arg_type_lookup_entry[arg_type_name[0]] = arg_type_lookup_entry.get(arg_type_name[0], 0) + 1
                self.arg_type_lookup[type_id] = arg_type_lookup_entry
            # end update arg_type_lookup

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
        if match is not None:
            result = {
                'name': match.group('name'),
                'return type': int(match.group('return_type')),
                '# arg': int(match.group('no_args')),
                'arg types': [int(s) for s in match.group('arg_types').split()],
            }
        else:
            match = self.get_calling_function_name_pattern().search(f2c_p_latter_line)
            result = {'name': match.group('name')}

        return result

    def find_any_missing_function(self):
        definition_missing = set(self.big_table.keys())
        never_called = set(self.big_table.keys())

        # function loop
        for function_name, function_info in self.big_table.items():
            if 'arg list' in function_info:
                definition_missing.remove(function_name)
            if 'arg types' in function_info:
                never_called.remove(function_name)

        return definition_missing, never_called


def scan_f2c():
    reader = F2cpReader()
    for lib, lib_path in f2c_path_dict['f2c'].items():
        for dir_name, dir_list, file_list in os.walk(lib_path):
            for file_name in file_list:
                if '.P' == os.path.splitext(file_name)[-1]:
                    reader.parse_f2c_p(os.path.join(dir_name, file_name))
    return reader


def main():
    # scan through f2c folders
    reader = scan_f2c()
    # argument_type_id vs argument_type lookup table
    pprint.pprint(reader.arg_type_lookup)
    # size of the big table
    print('total functions: %d' % len(reader.big_table))
    # find functions not defined or not used
    definition_missing_set, never_called_set = reader.find_any_missing_function()
    print('never used %d' % len(never_called_set))
    print(never_called_set)
    print('not defined %d' % len(definition_missing_set))
    print(definition_missing_set)


if __name__ == '__main__':
    main()
