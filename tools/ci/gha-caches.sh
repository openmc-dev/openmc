#!/bin/bash
set -ex

echo "NJOY_HASH=$(git ls-remote https://github.com/njoy/NJOY2016.git | cut -f 1 |\
    head -c 15)" >> $GITHUB_ENV

if [[ $DAGMC = 'y' ]]; then
    echo "MOAB_HASH=$(git ls-remote https://bitbucket.org/fathomteam/moab.git \
    --branch Version5.1.0 | cut -f 1 | head -c 15)" >> $GITHUB_ENV

    echo "DAGMC_HASH=$(git ls-remote https://github.com/svalinn/dagmc.git |\
    cut -f 1 | head -c 15)" >> $GITHUB_ENV
fi

if [[ $NCRYSTAL = 'y' ]]; then
    echo "NC_HASH=$(git ls-remote https://github.com/mctools/ncrystal \
    --branch develop | cut -f 1 | head -c 15)" >> $GITHUB_ENV
fi

if [[ $VECTFIT = 'y' ]]; then
    echo "PYBIND_HASH=$(git ls-remote https://github.com/pybind/pybind11 \
    --branch master | cut -f 1 | head -c 15)" >> $GITHUB_ENV

    echo "VF_HASH=$(git ls-remote https://github.com/liangjg/vectfit.git |\
    cut -f 1 | head -c 15)" >> $GITHUB_ENV
fi

if [[ $LIBMESH = 'y' ]]; then
    echo "NC_HASH=$(git ls-remote https://github.com/libmesh/libmesh \
    --branch v1.7.1 | cut -f 1 | head -c 15)" >> $GITHUB_ENV
fi

echo "MCPL_HASH=$(git ls-remote https://github.com/mctools/mcpl.git |\
    cut -f 1 | head -c 15)" >> $GITHUB_ENV
