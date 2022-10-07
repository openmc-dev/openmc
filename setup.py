#!/usr/bin/env python

import glob
import sys
import numpy as np

from setuptools import setup, find_packages
try:
    from Cython.Build import cythonize
    have_cython = True
except ImportError:
    have_cython = False


# Determine shared library suffix
if sys.platform == 'darwin':
    suffix = 'dylib'
else:
    suffix = 'so'

# Get version information from __init__.py. This is ugly, but more reliable than
# using an import.
with open('src/openmc/__init__.py', 'r') as f:
    version = f.readlines()[-1].split()[-1].strip("'")

kwargs = {
    'version': version,
    'packages': find_packages(exclude=['tests*']),
    'scripts': glob.glob('scripts/openmc-*'),

    # Data files and libraries
    'package_data': {
        'openmc.lib': ['libopenmc.{suffix}'],
        'openmc.data': ['*.{suffix}', 'mass16.txt', 'BREMX.DAT', 'half_life.json', '*.h5'],
        'openmc.data.effective_dose': ['*.txt']
    }
}

# If Cython is present, add resonance reconstruction and fast float_endf
if have_cython:
    kwargs.update({
        'ext_modules': cythonize('src/openmc/data/*.pyx'),
        'include_dirs': [np.get_include()]
    })

setup(**kwargs)
