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
with open('openmc/__init__.py', 'r') as f:
    version = f.readlines()[-1].split()[-1].strip("'")

kwargs = {
    'name': 'openmc',
    'version': version,
    'packages': find_packages(exclude=['tests*']),
    'scripts': glob.glob('scripts/openmc-*'),

    # Data files and librarries
    'package_data': {
        'openmc.capi': ['libopenmc.{}'.format(suffix)],
        'openmc.data': ['mass.mas12', 'fission_Q_data_endfb71.h5']
    },

    # Metadata
    'author': 'The OpenMC Development Team',
    'author_email': 'openmc-dev@googlegroups.com',
    'description': 'OpenMC',
    'url': 'https://github.com/mit-crpg/openmc',
    'classifiers': [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: End Users/Desktop',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering'
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],

    # Required dependencies
    'install_requires': [
        'numpy>=1.9', 'h5py', 'scipy', 'ipython', 'matplotlib',
        'pandas', 'lxml', 'uncertainties'
    ],

    # Optional dependencies
    'extras_require': {
        'test': ['pytest', 'pytest-cov'],
        'vtk': ['vtk', 'silomesh'],
    },
}

# If Cython is present, add resonance reconstruction capability
if have_cython:
    kwargs.update({
        'ext_modules': cythonize('openmc/data/reconstruct.pyx'),
        'include_dirs': [np.get_include()]
    })

setup(**kwargs)
