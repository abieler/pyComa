'''
run the script with:

python setupOSX.py build_ext --inplace

in the same directory as the createRay.pyx file is. This compiles the cython pyx
module into createRay.so which can directly imported from python
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy as np

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension('createRay',
                             ['createRay.pyx'],
                             include_dirs=[np.get_include()])])
