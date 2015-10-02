#!/usr/bin/env python

import glob
import os
try:
    from setuptools import setup
    have_setuptools = True
except ImportError:
    from distutils.core import setup
    have_setuptools = False

kwargs = {'name': 'openmc',
          'version': '0.7.0',
          'packages': ['openmc'],
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
        'install_requires': ['numpy', 'h5py', 'matplotlib'],

        # Optional dependencies
        'extras_require': {
            'pandas': ['pandas'],
            'vtk': ['vtk', 'silomesh'],
            'validate': ['lxml']
        }})

setup(**kwargs)
