#!/usr/bin/env python

import glob
import numpy as np
try:
    from setuptools import setup
    have_setuptools = True
except ImportError:
    from distutils.core import setup
    have_setuptools = False

try:
    from Cython.Build import cythonize
    have_cython = True
except ImportError:
    have_cython = False

kwargs = {'name': 'openmc',
          'version': '0.8.0',
          'packages': ['openmc', 'openmc.data', 'openmc.mgxs', 'openmc.model',
                       'openmc.stats'],
          'scripts': glob.glob('scripts/openmc-*'),

          # Metadata
          'author': 'Will Boyd',
          'author_email': 'wbinventor@gmail.com',
          'description': 'OpenMC Python API',
          'url': 'https://github.com/mit-crpg/openmc',
          'classifiers': [
              'Intended Audience :: Developers',
              'Intended Audience :: End Users/Desktop',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: MIT License',
              'Natural Language :: English',
              'Programming Language :: Python',
              'Topic :: Scientific/Engineering'
          ]}

if have_setuptools:
    kwargs.update({
        # Required dependencies
        'install_requires': ['six', 'numpy>=1.9', 'h5py', 'matplotlib'],

        # Optional dependencies
        'extras_require': {
            'decay': ['uncertainties'],
            'pandas': ['pandas>=0.17.0'],
            'sparse' : ['scipy'],
            'vtk': ['vtk', 'silomesh'],
            'validate': ['lxml']
        },

        # Data files
        'package_data': {
            'openmc.data': ['mass.mas12', 'fission_Q_data_endfb71.h5']
        },
    })

# If Cython is present, add resonance reconstruction capability
if have_cython:
    kwargs.update({
        'ext_modules': cythonize('openmc/data/reconstruct.pyx'),
        'include_dirs': [np.get_include()]
    })

setup(**kwargs)
