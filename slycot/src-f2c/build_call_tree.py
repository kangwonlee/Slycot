import subprocess
import os
import pprint
import sys


class CFlow(object):
    def __init__(self, cflow_path=os.path.join('/', 'usr', 'bin', 'cflow')):
        self.cflow_path = cflow_path
        self.calls_dict = {}
        self.called_dict = {}
        self.run_result = ''

    def run(self, cmd):
        cp = subprocess.run([self.cflow_path, cmd], stdout=subprocess.PIPE)
        self.run_result = cp.stdout

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
               callee_name = line.strip()
               callee_list.append(callee_name)
               caller_set = self.called_dict.get(callee_name, set())
               caller_set.add(key)
               self.called_dict[callee_name] = caller_set

        self.update_calls_dict(key, value)
        del self.calls_dict['del_this']

    def update_calls_dict(self, key, value):
        self.calls_dict.update({key: value})


def main():
    cflow = CFlow()
    cflow.run('SB02MT.c')
    print('calls dict '.ljust(60, '#'))
    pprint.pprint(cflow.calls_dict)
    print('called dict '.ljust(60, '='))
    pprint.pprint(cflow.called_dict)


if '__main__' == __name__:
    main()
