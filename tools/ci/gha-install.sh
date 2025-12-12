#!/bin/bash
set -ex

# Upgrade pip, pytest, numpy before doing anything else.
pip install --upgrade pip
pip install --upgrade pytest
pip install --upgrade numpy

# Install NJOY 2016
./tools/ci/gha-install-njoy.sh

# Install DAGMC if needed
if [[ $DAGMC = 'y' ]]; then
    ./tools/ci/gha-install-dagmc.sh
fi

# Install NCrystal and verify installation
pip install 'ncrystal>=4.1.0'
nctool --test

# Install vectfit for WMP generation if needed
if [[ $VECTFIT = 'y' ]]; then
    ./tools/ci/gha-install-vectfit.sh
fi

# Install libMesh if needed
if [[ $LIBMESH = 'y' ]]; then
    ./tools/ci/gha-install-libmesh.sh
fi

# Install MCPL
pip install mcpl

# For MPI configurations, make sure mpi4py and h5py are built against the
# correct version of MPI
if [[ $MPI == 'y' ]]; then
    pip install --no-binary=mpi4py mpi4py

    export CC=mpicc
    export HDF5_MPI=ON
    export HDF5_DIR=/usr/lib/x86_64-linux-gnu/hdf5/mpich
    pip install --no-binary=h5py h5py
fi

# Build and install OpenMC
python tools/ci/gha-install.py

