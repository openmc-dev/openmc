#!/bin/bash
set -ex

# Install NJOY 2016
./tools/ci/travis-install-njoy.sh

# Upgrade pip before doing anything else
pip install --upgrade pip

# pytest installed by default -- make sure we get latest
pip install --upgrade pytest

# Install mpi4py for MPI configurations
if [[ $MPI == 'y' ]]; then
    pip install --no-binary=mpi4py mpi4py
fi

# Build and install OpenMC executable
python tools/ci/travis-install.py

# Install Python API in editable mode
pip install -e .[test]

# For uploading to coveralls
pip install python-coveralls
