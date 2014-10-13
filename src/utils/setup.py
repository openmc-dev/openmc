#!/usr/bin/env python

from distutils.core import setup

setup(name='statepoint',
      version='0.6.0',
      description='OpenMC StatePoint',
      author='Paul Romano',
      author_email='paul.k.romano@gmail.com',
      url='https://github.com/mit-crpg/openmc',
      py_modules=['statepoint'])

setup(name='openmc',
      version='0.6.0',
      description='OpenMC Python API',
      author='Will Boyd',
      author_email='wbinventor@gmail.com',
      url='https://github.com/mit-crpg/openmc',
      packages=['openmc'])
