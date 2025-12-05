#!/bin/bash
set -ex

echo "NJOY_HASH=$(git ls-remote https://github.com/njoy/NJOY2016.git HEAD |\
    head -c 15)" >> $GITHUB_ENV

if [[ $DAGMC = 'y' ]]; then
    echo "MOAB_HASH=$(git ls-remote https://bitbucket.org/fathomteam/moab.git \
    --branch Version5.1.0 HEAD | head -c 15)" >> $GITHUB_ENV

    echo "DAGMC_HASH=$(git ls-remote https://github.com/svalinn/dagmc.git HEAD |\
    head -c 15)" >> $GITHUB_ENV
fi

if [[ $VECTFIT = 'y' ]]; then
    echo "PYBIND_HASH=$(git ls-remote https://github.com/pybind/pybind11 \
    --branch master HEAD | head -c 15)" >> $GITHUB_ENV

    echo "VF_HASH=$(git ls-remote https://github.com/liangjg/vectfit.git HEAD |\
    head -c 15)" >> $GITHUB_ENV
fi

if [[ $LIBMESH = 'y' ]]; then
    echo "LB_HASH=$(git ls-remote https://github.com/libmesh/libmesh \
    --branch v1.7.1 HEAD | head -c 15)" >> $GITHUB_ENV
fi
