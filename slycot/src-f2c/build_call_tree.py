import subprocess
import os
import pprint


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

        callee_list = None
        key = 'del_this'

        for line in result_lines:
            if not line.startswith(' '):
               # store result list
               self.calls_dict[key] = str(callee_list)
               # reset result list
               key = line.split()[0]
               callee_list = [line]
            else:
               callee_list.append(line.strip())

        self.calls_dict.update({key:callee_list})
        del self.calls_dict['del_this']


def main():
    cflow = CFlow()
    cflow.run('SB02MT.c')
    pprint.pprint(cflow.calls_dict)


if '__main__' == __name__:
    main()
