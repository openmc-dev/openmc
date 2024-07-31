#!/usr/bin/env python

import numpy as np
from setuptools import setup
from Cython.Build import cythonize


kwargs = {
    # Cython is used to add resonance reconstruction and fast float_endf
    'ext_modules': cythonize('openmc/data/*.pyx'),
    'include_dirs': [np.get_include()]
}

setup(**kwargs)
