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
        self.p_file_name = ''
        self.lib_name = ''
        self.p_file_path = ''

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

    def parse_f2c_p(self, f2c_p_file_path, b_verbose=False):

        self.get_lib_name_from_p_file_path(f2c_p_file_path)

        with open(f2c_p_file_path) as f:
            lines = f.readlines()
        # first line : c definitions
        # second line and after : list of other functions called

        for line in lines:
            line = line.strip()

            if not line.startswith('/*'):
                # functions defined
                info = self.find_function_info(line)
                if b_verbose:
                    print(info)
            else:
                # functions used inside
                info = self.find_calling_function_info(line)
            self.update_big_table(info)

    def get_lib_name_from_p_file_path(self, f2c_p_file_path):
        """

        :param str f2c_p_file_path: path to the f2c P file
        :return:
        """

        path_lower = f2c_p_file_path.lower()

        # because lapack is a subfolder of slycot/f2c
        # and blas is a a subfolder of lapack
        if 'blas' in path_lower:
            self.lib_name = 'blas'
        elif 'lapack' in path_lower:
            self.lib_name = 'lapack'
        elif 'slycot' in path_lower:
            self.lib_name = 'slycot'
        else:
            raise ValueError('library unknown')

    def update_big_table(self, info_dict):
        # begin update big table using the info_dict
        big_table_entry = self.big_table.get(info_dict['name'], {})

        # if already know return type in a string, do not update
        if ('return type' in big_table_entry) and (isinstance(big_table_entry['return type'], str)):
            info_dict.pop('return type', None)
        # end if already know return type in a string, do not update

        big_table_entry.update(info_dict)
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
            'lib': self.lib_name,
            'path': self.p_file_path,
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
                'lib': self.lib_name,
                'path': self.p_file_path,
            }
        else:
            match = self.get_calling_function_name_pattern().search(f2c_p_latter_line)
            result = {'name': match.group('name')}

        return result

    def find_any_missing_function(self):
        definition_missing_set = set(self.big_table.keys())
        never_called_set = set(self.big_table.keys())

        # function loop
        for function_name, function_info in self.big_table.items():
            if 'arg list' in function_info:
                # found definition
                definition_missing_set.remove(function_name)
            if 'arg types' in function_info:
                # found usage
                never_called_set.remove(function_name)

        definition_missing_dict = {}
        for function_name in definition_missing_set:
            definition_missing_dict[function_name] = self.big_table[function_name]
            definition_missing_dict[function_name].pop('arg types', None)

        never_called_dict = {}
        for function_name in never_called_set:
            never_called_dict[function_name] = self.big_table[function_name]
            never_called_dict[function_name].pop('arg types', None)

        return definition_missing_dict, never_called_dict


class Dict2MDTable(object):
    """
    >>> table = Dict2MDTable(
        {   # table data
            'a': {'b': 1, 'c': 2, 'd': 3},
            'e': {'b': 4, 'c': 5, 'd': 6},
        },
        [   # column order
            {'name':'b'}, {'name':'c', 'align': 'right'}, {'name':'d', 'align': 'left'}
        ],
    )
    >>> print(table)
    |    | b | c | d |
    |:-----:|:-----:|------:|:------|
    | a | 1 | 2 | 3 |
    | e | 4 | 5 | 6 |
    """
    align = {
        'right': '------:',
        'center': ':-----:',
        'left': ':------',
    }

    def __init__(self, input_dict, column_order_list=None, row_selection_list=None):
        """

        :param input_dict:
        :param column_order_list:
        :param list | tuple | set row_selection_list: if not given, all rows
        """
        self.input_dict = input_dict

        # to select rows of interest
        if row_selection_list is None:
            self.row_selection_list = list(self.input_dict.keys())
        elif isinstance(row_selection_list, (list, tuple, set)):
            self.row_selection_list = row_selection_list
        else:
            raise ValueError('expect row_selection_list to be one of list, tuple, and set')

        # column selection
        if column_order_list is None:
            one_sample = self.input_dict[set(self.input_dict.keys()).pop()]
            self.column_order_list = {'name': key for key in one_sample}
        elif isinstance(column_order_list, (list, tuple)):
            self.column_order_list = column_order_list
        else:
            raise ValueError('expect column_order_list to be a list or tuple')

    def first_row(self):
        space = '    '
        row_list = ['', space]
        for column in self.column_order_list:
            row_list.append(' %s ' % column.get('name', space))

        row_list.append('')

        result = '|'.join(row_list)

        return result

    def second_row(self):
        row_list = ['', self.align['center']]
        for column in self.column_order_list:
            row_list.append(self.align[column.get('align', 'center')])
        row_list.append('')
        result = '|'.join(row_list)

        return result

    def third_and_latter_row(self):
        row_list = []

        # row loop
        for function_name in self.row_selection_list:
            function_info_dict = self.input_dict[function_name]

            column_list = self.get_column_list_third_and_latter_row(function_info_dict, function_name)

            row_text = ' '.join(column_list)
            # this completes one row

            row_list.append(row_text)
        # this completes all rows

        result = '\n'.join(row_list)

        return result

    def get_column_list_third_and_latter_row(self, function_info_dict, function_name):
        column_list = ['|', str(function_name), '|']
        # first column
        # loop for the following columns
        for column in self.column_order_list:
            column_list.append(str(function_info_dict.get(column['name'], '')))
            column_list.append('|')
        return column_list

    def __str__(self):
        table_list = [
            self.first_row(),
            self.second_row(),
            self.third_and_latter_row(),
        ]

        result = '\n'.join(table_list)

        return result


def scan_f2c():
    reader = F2cpReader()
    for lib, lib_path in f2c_path_dict['f2c'].items():
        reader.lib_name = lib
        for dir_name, dir_list, file_list in os.walk(lib_path):
            for file_name in file_list:
                if '.P' == os.path.splitext(file_name)[-1]:
                    reader.p_file_path = dir_name
                    reader.p_file_name = file_name
                    reader.parse_f2c_p(os.path.join(dir_name, file_name))
    return reader


def main():
    # scan through f2c folders
    reader = scan_f2c()
    # argument_type_id vs argument_type lookup table
    pprint.pprint(reader.arg_type_lookup)
    # size of the big table
    print('total functions: %d\n' % len(reader.big_table))
    big_table_printer = Dict2MDTable(
        reader.big_table,
        [{'name': 'return type'}, {'name': '# arg'}, {'name': 'arg list'}, {'name': 'lib'},
         {'name': 'path', 'align': 'left'},
         ],
        ['sb03md_', 'sb04md_', 'sg03ad_', 'sb04qd_', 'sb02md_', 'sb02mt_', 'sg02ad_', 'ab09md_', 'ab09md_', 'ab09nd_',
         'sb10hd_', 'sb10hd_', 'sb10hd_', 'sb03od_', 'tb01pd_', 'td04ad_', 'td04ad_', 'sb02od_', ],
    )
    print(big_table_printer)

    # find functions not defined or not used
    definition_missing, never_called = reader.find_any_missing_function()

    # never called table
    print('never used %d\n' % len(never_called))
    never_called_table_converter = Dict2MDTable(
        never_called,
        [{'name': 'lib'}, {'name': '# arg'}, {'name': 'return type'}, {'name': 'path'}, ]
    )
    print(never_called_table_converter)

    # not defined table
    print('not defined %d\n' % len(definition_missing))
    not_defined_table_converter = Dict2MDTable(
        definition_missing,
        [{'name': 'lib'}, {'name': '# arg'}, {'name': 'return type'}, {'name': 'path'}, ]
    )
    print(not_defined_table_converter)


if __name__ == '__main__':
    main()
