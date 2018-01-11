# to build cython code
# http://cython.readthedocs.io/en/latest/src/quickstart/build.html#building-a-cython-module-using-distutils
from distutils.core import setup, Extension

import numpy
from Cython.Distutils import build_ext

# http://cython.readthedocs.io/en/latest/src/tutorial/clibraries.html#compiling-and-linking

sourcefiles = ['analysis.pyx', 'SB03MD.c']

setup(cmdclass={'build_ext': build_ext},
      ext_modules=[Extension('analysis',
                             sources=sourcefiles,
                             include_dirs=[numpy.get_include()])],
      )
