import glob
import os
import subprocess
import sys


class CFlow(object):
    def __init__(self, cflow_path=os.path.join('/', 'usr', 'bin', 'cflow')):
        self.cflow_path = cflow_path
        self.calls_dict = {}
        self.called_dict = {}
        self.result_str = ''

    def run(self, cmd):
        self.run_cmd_list([self.cflow_path, cmd])

        result_lines = self.result_str.splitlines()

        value = None
        key = 'del_this'

        for line in result_lines:
            if not line.startswith(' '):
               # store result list
               self.update_calls_dict(key, value)
               # reset result list
               key = line.split()[0]
               callee_list = []
               value = {'info':line, 'list':callee_list}
            else:
               callee_name = line.strip().rstrip('_()')
               callee_list.append(callee_name)
               caller_set = self.called_dict.get(callee_name, set())
               caller_set.add(key)
               self.called_dict[callee_name] = caller_set

        self.update_calls_dict(key, value)
        del self.calls_dict['del_this']

    def update_calls_dict(self, key, value):
        self.calls_dict.update({key: value})

    def run_files(self, files):
        for file in files:
            self.run(file)

    def run_files_altogether(self, files):
        command_list = [self.cflow_path, ] + list(files)
        self.run_cmd_list(command_list)

        return self.result_str

    def run_cmd_list(self, command_list):
        cp = subprocess.run(command_list, stdout=subprocess.PIPE)
        self.result_str = cp.stdout.decode()

    def run_files_altogether_main(self, main_function_name, files):
        command_list = [self.cflow_path, '--main', main_function_name] + list(files)
        self.run_cmd_list(command_list)

        return self.result_str


def tree_path(files):
    cflow = CFlow()
    cflow.run_files_altogether(files)

    return cflow


def main():
    if 3 <= len(sys.argv):
        lapack_path = sys.argv[1]
        blas_path = sys.argv[2]

    slycot_function_list_filename = 'slycot_functions_list.txt'
    not_fully_expanded_function_list_filename = 'still_not_expanded.txt'

    with open(slycot_function_list_filename, 'rt') as python_control_slycot_file:
        slycot_set = set([line.strip().upper() +'.c' for line in python_control_slycot_file.readlines()])

    # build top level call tree
    cflow = tree_path(slycot_set)

    print('calls dict (%4d)'.ljust(60, '#') % len(cflow.calls_dict))
    # pprint.pprint(cflow.calls_dict)
    print('called dict (%4d)'.ljust(60, '=') % len(cflow.called_dict))
    # pprint.pprint(cflow.called_dict)

    slycot_files_set = glob.glob('*.c')
    print("# slycot files =", len(slycot_files_set))
    print(list(slycot_files_set)[:10])

    lapack_files_set = glob.glob(os.path.join(lapack_path, '*.c'))
    print("# lapack files =", len(lapack_files_set))
    print(list(lapack_files_set)[:10])

    blas_files_set = glob.glob(os.path.join(blas_path, '*.c'))
    print("# blas files =", len(blas_files_set))
    print(list(blas_files_set)[:10])

    slycot_dict = {}
    lapack_dict = {}
    blas_dict = {}

    slycot_function_path_set = set()
    lapack_function_path_set = set()
    blas_function_path_set = set()
    unknown_set = set()

    function_name_set = set(cflow.called_dict.keys())

    # initialize previously_unknown_set
    previously_unknown_set = init_prev_unknown_set(not_fully_expanded_function_list_filename)

    function_name_set.update(previously_unknown_set)

    for function_name in function_name_set:
        slycot_function_path = get_slycot_function_path_uppercase(function_name)
        lapack_function_path = get_function_path(lapack_path, function_name)
        blas_function_path = get_function_path(blas_path, function_name)

        if slycot_function_path in slycot_files_set:
            slycot_function_path_set.add(slycot_function_path)
        elif lapack_function_path in lapack_files_set:
            lapack_function_path_set.add(lapack_function_path)
        elif blas_function_path in blas_files_set:
            blas_function_path_set.add(blas_function_path)
        else:
            unknown_set.add(function_name)

    print("# functions in slycot =", len(slycot_dict))
    print("# functions in lapack =", len(lapack_dict))
    print("# functions in blas =", len(blas_dict))
    print("# functions unknown =", len(cflow.called_dict))
    print(sorted(list(unknown_set)))

    # build extended function set
    result_replaced = build_extended_call_tree(slycot_set, slycot_function_path_set, lapack_function_path_set,
                                               blas_function_path_set, blas_path, lapack_path, cflow)

    print(result_replaced)

    not_fully_expanded_function_set = find_not_fully_expanded_functions(result_replaced, slycot_set, unknown_set)

    print('may need to run for these (%d)' % len(not_fully_expanded_function_set),
          sorted(list(not_fully_expanded_function_set)))

    # write not fully expanded functions to a separate file
    with open(not_fully_expanded_function_list_filename, 'w') as output_file:
        for function_name in sorted(list(not_fully_expanded_function_set.union(previously_unknown_set))):
            output_file.write('%s\n' % function_name)


def init_prev_unknown_set(not_fully_expanded_function_list_filename):
    if os.path.exists(not_fully_expanded_function_list_filename):
        with open(not_fully_expanded_function_list_filename, 'r') as expand_these_file:
            previously_unknown_set = set([line.strip() for line in expand_these_file.readlines()])
    else:
        previously_unknown_set = set()
    return previously_unknown_set


def find_not_fully_expanded_functions(result_replaced, slycot_set, unknown_set):
    may_need_to_add_set = set()
    for line in result_replaced.splitlines():
        line_split = line.split()
        if 1 < len(line_split):
            continue
        else:
            not_fully_expanded_function_name = line.strip().rstrip('()').rstrip('_')
            if (not_fully_expanded_function_name not in unknown_set) and (
                        not_fully_expanded_function_name not in slycot_set):
                may_need_to_add_set.add(not_fully_expanded_function_name)
    return may_need_to_add_set


def build_extended_call_tree(slycot_set, slycot_function_path_set, lapack_function_path_set, blas_function_path_set,
                             blas_path, lapack_path, cflow):
    big_set = slycot_function_path_set.union(lapack_function_path_set.union(blas_function_path_set.union(slycot_set)))
    result = cflow.run_files_altogether(big_set)
    result_replaced = result.replace(blas_path, '[BLAS]').replace(lapack_path, '[LAPACK]')
    with open('python-control-slycot_call_tree.txt', 'wt') as output_file:
        output_file.write(result_replaced)
    return result_replaced


def get_function_path(function_folder, c_function_name):
    function_c_file_name = c_function_name + '.c'
    return os.path.join(function_folder, function_c_file_name)


def get_slycot_function_path_uppercase(c_function_name):
    return c_function_name.upper() + '.c'


if '__main__' == __name__:
    if 2 > len(sys.argv):
        with open('call_tree.ini', 'rt') as argv:
            sys.argv.extend([line.strip() for line in argv.readlines()])
    main()
