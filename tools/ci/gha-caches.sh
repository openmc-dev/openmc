#!/bin/bash
set -ex

echo "NJOY_HASH=$(git ls-remote https://github.com/njoy/NJOY2016.git | cut -f 1 |\
    head -c 15)" >> $GITHUB_ENV

if [[ $DAGMC = 'y' ]]; then
    echo "MOAB_HASH=$(git ls-remote https://bitbucket.org/fathomteam/moab.git |\
    cut -f 1 | head -c 15)" >> $GITHUB_ENV

    echo "DAGMC_HASH=$(git ls-remote https://github.com/svalinn/dagmc.git |\
     cut -f 1 | head -c 15)" >> $GITHUB_ENV
fi

# Check for NCrystal cache if needed
if [[ $NCRYSTAL = 'y' ]]; then
    echo "NC_HASH=$(git ls-remote https://github.com/mctools/ncrystal \
    --branch develop | cut -f 1 | head -c 15)" >> $GITHUB_ENV
fi

# Install vectfit for WMP generation if needed
#if [[ $VECTFIT = 'y' ]]; then
#    ./tools/ci/gha-install-vectfit.sh
#fi

# Check for libmesh cache if needed
#if [[ $LIBMESH = 'y' ]]; then
#    ./tools/ci/gha-install-libmesh.sh
#fi

echo "MCPL_HASH=$(git ls-remote https://github.com/mctools/mcpl.git | cut -f 1 |\
    head -c 15)" >> $GITHUB_ENV

# may be able to use this to control the directories to be cached.