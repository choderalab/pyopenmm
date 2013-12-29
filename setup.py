"""
pyopenmm

A pure Python implementation of the OpenMM molecular dynamics API and engine.

"""

from __future__ import print_function
DOCLINES = __doc__.split("\n")

import os
import sys
import shutil
import tempfile
import subprocess
from distutils.ccompiler import new_compiler
from setuptools import setup, Extension

import numpy
try:
    from Cython.Distutils import build_ext
    setup_kwargs = {'cmdclass': {'build_ext': build_ext}}
    cython_extension = 'pyx'
except ImportError:
    setup_kwargs = {}
    cython_extension = 'c'



##########################
VERSION = "0.1"
ISRELEASED = False
__version__ = VERSION
##########################


CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Lesser General Public License (LGPL), Version 3
Programming Language :: C
Programming Language :: OpenCL
Programming Language :: CUDA
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""

extensions = []

setup(name='pyopenmm',
      author='John D. Chodera',
      author_email='jchodera@gmail.com',
      description=DOCLINES[0],
      long_description="\n".join(DOCLINES[2:]),
      version=__version__,
      license='LGPLv3+',
      url='http://github.com/choderalab/pyopenmm',
      platforms=['Linux', 'Mac OS-X', 'Unix'],
      classifiers=CLASSIFIERS.splitlines(),
      packages=["pyopenmm"],
      install_requires=['numpy', 'pyopencl', 'nose', 'nose-exclude', 'mako'],
      zip_safe=False,
      scripts=[],
      ext_modules=extensions,
      package_data={'pyopenmm': ['data/*/*']},  # Install all data directories of the form testsystems/data/X/
      **setup_kwargs
      )
