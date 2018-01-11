# Valentin Haenel, 2.8.2.1. Example, 2.8.2. Python-C-Api, Scipy Lectures, Oct 18 2016, [Online]
#   Available: http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#example
# Valentin Haenel, 2.8.4.1. Example, 2.8.4. SWIG, Scipy Lectures, Oct 18 2016, [Online]
#   Available: , http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id8
# Valentin Haenel, 2.8.5.1. Example, 2.8.5. Cython, Scipy Lectures, Oct 18 2016, [Online]
#   Available: , http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id10
# Valentin Haenel, 2.8.5.2. Numpy Support, 2.8.5. Cython, Scipy Lectures, Oct 18 2016, [Online]
#   Available: , http://www.scipy-lectures.org/advanced/interfacing_with_c/interfacing_with_c.html#id13

from distutils.core import setup, Extension

import numpy
from Cython.Distutils import build_ext

print('for NumPy Support of Cython '.ljust(60, '#'))
setup(cmdclass={'build_ext': build_ext},
      ext_modules=[Extension("cos_cython_numpy",
                             sources=['_cos_cython_numpy.pyx', "cos_cython_numpy.c"],
                             include_dirs=[numpy.get_include()])],
      )
