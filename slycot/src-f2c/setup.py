# to build cython code
# http://cython.readthedocs.io/en/latest/src/quickstart/build.html#building-a-cython-module-using-distutils
from distutils.core import setup

from Cython.Build import cythonize

setup(
    name='Hello world app',
    ext_modules=cythonize('first.pyx')
)
