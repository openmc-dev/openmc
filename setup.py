#!/usr/bin/env python

import os
from distutils.core import setup

setup(name='openmc',
      version='0.6.2',
      description='OpenMC Python API',
      author='Will Boyd',
      author_email='wbinventor@gmail.com',
      url='https://github.com/mit-crpg/openmc',
      packages=['openmc'],
      scripts=[os.path.join('scripts', f) for f in os.listdir('scripts')])
