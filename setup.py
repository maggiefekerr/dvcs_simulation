from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    name='dvcs_sim',
    ext_modules=cythonize("km15gen.pyx", compiler_directives={'language_level': "3"}),
    include_dirs=[numpy.get_include()],
    scripts=['main.py'],
)