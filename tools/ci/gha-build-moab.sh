#!/bin/bash
set -ex

MOAB_BRANCH='5.5.1'
MOAB_REPO='https://bitbucket.org/fathomteam/moab/'
MOAB_INSTALL_DIR=$HOME/MOAB/

pushd $HOME

if [[ -d $MOAB_INSTALL_DIR/lib ]]; then
    popd
    exit 0
fi

mkdir -p MOAB
cd MOAB
git clone -b $MOAB_BRANCH $MOAB_REPO moab
mkdir build
cd build
cmake ../moab \
    -DENABLE_HDF5=ON \
    -DENABLE_NETCDF=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DCMAKE_INSTALL_PREFIX=$MOAB_INSTALL_DIR
make -j4
make -j4 install
rm -rf $HOME/MOAB/moab $HOME/MOAB/build

popd
