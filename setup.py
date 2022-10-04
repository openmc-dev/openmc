#!/usr/bin/env python

import glob
import sys

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

setup(
    name="openmc",
    version=version,
    packages=find_packages(exclude=['tests*']),
    include_package_data=True,
    python_requires=">=3.7",
    install_requires=[
        "matplotlib",
        "numpy",
        "scipy",
        "ipython",
        "matplotlib",
        "uncertainties",
        "lxml",
        "pandas",
        "h5py"],
    scripts=glob.glob('scripts/openmc-*'),
    package_data={
        'openmc.lib': ['libopenmc.{}'.format(suffix)],
        'openmc.data': ['mass16.txt', 'BREMX.DAT', 'half_life.json', '*.h5'],
        'openmc.data.effective_dose': ['*.txt']
    },
    ext_modules= cythonize('openmc/data/*.pyx')
)
