#!/usr/bin/env python

import numpy as np
from setuptools import setup
from Cython.Build import cythonize
import sys


kwargs = {
    'include_dirs': [np.get_include()]
}

if sys.platform != 'win32':
    # Cython is used to add resonance reconstruction and fast float_endf, but
    # it can only be compiled on Linux/MacOS systems currently.
    kwargs['ext_modules'] = cythonize('openmc/data/*.pyx')

setup(**kwargs)
