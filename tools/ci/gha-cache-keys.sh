#!/bin/bash
set -ex

echo "NJOY_HASH=$(git ls-remote https://github.com/njoy/NJOY2016.git HEAD \
    --tags 2016.78 | head -c 15)" >> $GITHUB_ENV

if [[ "$DAGMC" == 'y' ]]; then
    echo "MOAB_HASH=$(git ls-remote https://bitbucket.org/fathomteam/moab.git \
    --tags Version5.1.0 | head -c 15)" >> $GITHUB_ENV

    echo "DAGMC_HASH=$(git ls-remote https://github.com/svalinn/dagmc.git HEAD |\
    head -c 15)" >> $GITHUB_ENV
fi

if [[ "$LIBMESH" == 'y' ]]; then
    echo "LB_HASH=$(git ls-remote https://github.com/libmesh/libmesh \
    --tags v1.7.1 | head -c 15)" >> $GITHUB_ENV
fi

if [[ "$MPI" == 'y' ]]; then
    dpkg-query -W libmpich-dev libhdf5-mpich-dev > $GITHUB_WORKSPACE/tools/ci/mpi-deps.txt
    cat $GITHUB_WORKSPACE/tools/ci/mpi-deps.txt
    sha256sum $GITHUB_WORKSPACE/tools/ci/requirements-mpi-y.txt
    sha256sum $GITHUB_WORKSPACE/tools/ci/mpi-deps.txt
else
    touch $GITHUB_WORKSPACE/tools/ci/mpi-deps.txt
    sha256sum $GITHUB_WORKSPACE/tools/ci/requirements-mpi-n.txt
    sha256sum $GITHUB_WORKSPACE/tools/ci/mpi-deps.txt
fi