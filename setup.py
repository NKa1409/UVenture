from distutils.core import setup
from Cython.Build import cythonize
import numpy


setup(
    ext_modules=cythonize(
        ["c_functions/MS_functions.pyx", 
         "c_functions/peakdetection_funcs.pyx",
         "c_functions/UVenture.pyx"], compiler_directives={"language_level": "3"}
    ),
    include_dirs=[numpy.get_include()],
)