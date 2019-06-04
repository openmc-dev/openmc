#!/bin/bash
set -ex

# Install NJOY 2016
./tools/ci/travis-install-njoy.sh

# Install DAGMC if needed
if [[ $DAGMC = 'y' ]]; then
    ./tools/ci/travis-install-dagmc.sh
fi

# Upgrade pip, pytest, numpy before doing anything else
pip install --upgrade pip
pip install --upgrade pytest
pip install --upgrade numpy

# Install mpi4py for MPI configurations
if [[ $MPI == 'y' ]]; then
    pip install --no-binary=mpi4py mpi4py
fi

# Build and install OpenMC executable
python tools/ci/travis-install.py

# Install Python API in editable mode
pip install -e .[test,vtk]

# For coverage testing of the C++ source files
pip install cpp-coveralls

# For coverage testing of the Python source files
pip install coveralls
