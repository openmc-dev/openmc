#!/usr/bin/env python

import os
import numpy as np
from setuptools import setup, Extension


class OpenMCExtension(Extension):
    def __init__(self, name, cmake_lists_dir=".", sources=[], **kwa):
        Extension.__init__(self, name, sources=sources, **kwa)
        self.cmake_lists_dir = os.path.abspath(cmake_lists_dir)


kwargs = {
    'ext_modules': [OpenMCExtension('libopenmc')],
    'include_dirs': [np.get_include()]
}

setup(**kwargs)
