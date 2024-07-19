#!/usr/bin/env python

import os
import numpy as np
from setuptools import setup, Extension
from Cython.Build import cythonize


class OpenMCExtension(Extension):
    def __init__(self, name, cmake_lists_dir=".", sources=[], **kwa):
        Extension.__init__(self, name, sources=sources, **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)


kwargs = {
    # Cython is used to add resonance reconstruction and fast float_endf
    'ext_modules': cythonize('openmc/data/*.pyx') + [OpenMCExtension('libopenmc')],
    'include_dirs': [np.get_include()]
}

setup(**kwargs)
