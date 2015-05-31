#!/usr/bin/env python

import glob
import os
from setuptools import setup, find_packages

setup(name='openmc',
      version='0.6.2',
      packages=find_packages(),
      scripts=glob.glob('scripts/openmc-*'),

      # Required dependencies
      install_requires=['numpy', 'scipy', 'h5py', 'matplotlib'],

      # Optional dependencies
      extras_require={
          'pandas': ['pandas'],
          'vtk': ['vtk', 'silomesh'],
          'validate': ['lxml']
      },

      # Metadata
      author='Will Boyd',
      author_email='wbinventor@gmail.com',
      description='OpenMC Python API',
      url='https://github.com/mit-crpg/openmc',
      classifiers=[
          'Intended Audience :: Developers',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Natural Language :: English',
          'Programming Language :: Python',
          'Topic :: Scientific/Engineering'
      ]
)
