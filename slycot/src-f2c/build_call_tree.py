import subprocess
import os



class CFlow(object):
    def __init__(self, cflow_path=os.path.join('/', 'usr', 'bin', 'cflow')):
        self.cflow_path = cflow_path

    def run(self, cmd):
        cp = subprocess.run([self.cflow_path, cmd])
        print(cp.stdout)


def main():
    cflow = CFlow()
    cflow.run('SB02MT.c')


if '__main__' == __name__:
    main()
