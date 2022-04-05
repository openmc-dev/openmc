#!/bin/bash
set -ex

# Upgrade pip, pytest, numpy before doing anything else.
# TODO: numpy 1.22 results in several failing tests, so we force a lower version
# for now (similar change made in pyproject.toml). When this is removed, those
# tests will need to be updated.
pip install --upgrade pip
pip install --upgrade pytest
pip install --upgrade "numpy<1.22"

# Install NJOY 2016
./tools/ci/gha-install-njoy.sh

# Install DAGMC if needed
if [[ $DAGMC = 'y' ]]; then
    ./tools/ci/gha-install-dagmc.sh
fi

# Install vectfit for WMP generation if needed
if [[ $VECTFIT = 'y' ]]; then
    ./tools/ci/gha-install-vectfit.sh
fi

# Install libMesh if needed
if [[ $LIBMESH = 'y' ]]; then
    ./tools/ci/gha-install-libmesh.sh
fi

# For MPI configurations, make sure mpi4py and h5py are built against the
# correct version of MPI
if [[ $MPI == 'y' ]]; then
    pip install --no-binary=mpi4py mpi4py

    export CC=mpicc
    export HDF5_MPI=ON
    export HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/mpich
    pip install --no-binary=h5py h5py
fi

# Build and install OpenMC executable
python tools/ci/gha-install.py

# Install Python API in editable mode
pip install -e .[test,vtk]

# For coverage testing of the C++ source files
pip install cpp-coveralls

# For coverage testing of the Python source files
pip install coveralls
