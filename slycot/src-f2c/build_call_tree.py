import subprocess
import os
import pprint
import sys
import glob


class CFlow(object):
    def __init__(self, cflow_path=os.path.join('/', 'usr', 'bin', 'cflow')):
        self.cflow_path = cflow_path
        self.calls_dict = {}
        self.called_dict = {}

    def run(self, cmd):
        cp = subprocess.run([self.cflow_path, cmd], stdout=subprocess.PIPE)

        result_lines = cp.stdout.decode().splitlines()

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


def tree_path(path):
    cflow = CFlow()
    files = set(glob.glob(os.path.join(path,'*.c')))
    cflow.run_files(files)

    return cflow

def main():
    if 3 <= len(sys.argv):
        lapack_path = sys.argv[1]
        blas_path = sys.argv[2]

    cflow = tree_path(os.curdir)

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

    function_name_set = set(cflow.called_dict.keys())

    for function_name in function_name_set:
        function_c_file_name = function_name + '.c'
        if function_name.upper() + '.c' in slycot_files_set:
            slycot_dict[function_name] = cflow.called_dict[function_name]
            del cflow.called_dict[function_name]
        elif os.path.join(lapack_path,function_c_file_name) in lapack_files_set:
            lapack_dict[function_name] = cflow.called_dict[function_name]
            del cflow.called_dict[function_name]
        elif os.path.join(blas_path, function_c_file_name) in blas_files_set:
            blas_dict[function_name] = cflow.called_dict[function_name]
            del cflow.called_dict[function_name]

    print("# functions in slycot =", len(slycot_dict))
    print("# functions in lapack =", len(lapack_dict))
    print("# functions in blas =", len(blas_dict))
    print("# functions unknown =", len(cflow.called_dict))
    print(sorted(list(cflow.called_dict.keys())))


if '__main__' == __name__:
    if 2 > len(sys.argv):
        with open('call_tree.ini', 'rt') as argv:
            sys.argv.extend([line.strip() for line in argv.readlines()])
    main()
