#!/usr/bin/env python

import glob
import sys
import numpy as np

from setuptools import find_packages

from skbuild import setup
from Cython.Build import cythonize


# Determine shared library suffix
if sys.platform == 'darwin':
    suffix = 'dylib'
else:
    suffix = 'so'

# Get version information from __init__.py. This is ugly, but more reliable
# than using an import.
with open('openmc/__init__.py', 'r') as f:
    version = f.readlines()[-1].split()[-1].strip("'")

kwargs = {
    'name': 'openmc',
    'version': version,
    'packages': find_packages(exclude=['tests*']),
    'scripts': glob.glob('scripts/openmc-*'),

    # Data files and libraries
    'package_data': {
        'openmc.lib': ['libopenmc.{}'.format(suffix)],
        'openmc.data': ['*.pyx', '*.c', 'mass16.txt', 'BREMX.DAT', 'half_life.json', '*.h5'],
        'openmc.data.effective_dose': ['*.txt']
    },
    'python_requires': '>=3.7',
    'include_package_data': True,
    'install_requires': [
        'numpy>=1.9', 'h5py', 'scipy', 'ipython', 'matplotlib',
        'pandas', 'lxml', 'uncertainties', 'cython'
    ],
    'ext_modules': cythonize('openmc/data/*.pyx'),
    'include_dirs': [np.get_include()]
}

setup(**kwargs)
