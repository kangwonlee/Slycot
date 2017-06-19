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
    if 2 <= len(sys.argv):
        lapack_path = sys.argv[1]

    cflow = tree_path(os.curdir)

    print('calls dict '.ljust(60, '#'))
    pprint.pprint(cflow.calls_dict)
    print('called dict '.ljust(60, '='))
    pprint.pprint(cflow.called_dict)


    lapack_files_set = glob.glob(os.path.join(lapack_path, '*.c'))
    print("# lapack files =", len(lapack_files_set))


if '__main__' == __name__:
    main()
