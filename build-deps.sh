#!/bin/bash
set -eu

# pip install numpy==1.25.1 cython==0.29.36 # required for MOAB install (silly)

# It's important to use LLVM clang rather than Apple clang
# CC=/opt/homebrew/Cellar/llvm/16.0.6/bin/clang
# CXX=/opt/homebrew/Cellar/llvm/16.0.6/bin/clang++

NUM_CORES=2

# 5.5.0 breaks with pyne std::isnan bug
MOAB_TAG='5.3.0'
MOAB_REPO='https://bitbucket.org/fathomteam/moab/'

DAGMC_TAG='v3.2.2'
DAGMC_REPO='https://github.com/svalinn/DAGMC'

# MOAB
# WARNING: CMake often gets wrong python path. Check logs in case of error.
git clone  --single-branch -b ${MOAB_TAG} --depth 1 ${MOAB_REPO}
pushd moab
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DENABLE_HDF5=ON \
    -DENABLE_NETCDF=OFF \
    -DBUILD_SHARED_LIBS=ON \
    -DENABLE_FORTRAN=OFF \
    -DENABLE_BLASLAPACK=OFF \
    -DENABLE_PYMOAB=OFF
cmake --build build -j $NUM_CORES
sudo cmake --build build -t install
popd

# DAGMC
git clone --shallow-submodules --recurse-submodules --single-branch -b ${DAGMC_TAG} --depth 1 ${DAGMC_REPO}
pushd DAGMC
# MANUAL: Fix isnan problem in pyne.h and pyne.cpp, changing isnan to std::isnan
cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_TALLY=ON \
    -DMOAB_DIR=/usr/local \
    -DBUILD_STATIC_LIBS=OFF
sudo cmake --build build -j $NUM_CORES -t install
popd

# OpenMC
# build dir is claimed by bdist_wheel
cmake -S . -B bld \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_BUILD_TESTS=OFF \
    -DOPENMC_USE_DAGMC=ON
cmake --build bld -j 8