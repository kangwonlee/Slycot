# SWIG
# http://intermediate-and-advanced-software-carpetnry.readthedocs.io/en/latest/c++-wrapping.html
all:
    # run WIG to generate the C code extension
    swig -python -c++ -o SB02MT.c sb02mt.i
    # actually build it
    python setup.py build_ext --inplace
