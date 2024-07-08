#!/usr/bin/env python

import glob
import sys
import numpy as np

from setuptools import setup, find_packages
from Cython.Build import cythonize


kwargs = {
    # Cython is used to add resonance reconstruction and fast float_endf
    'ext_modules': cythonize('openmc/data/*.pyx'),
    'include_dirs': [np.get_include()]
}

setup(**kwargs)
