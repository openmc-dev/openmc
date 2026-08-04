#!/bin/bash
set -ex

# DAGMC Variables
DAGMC_BRANCH='develop'
DAGMC_REPO='https://github.com/svalinn/dagmc'
DAGMC_INSTALL_DIR=$HOME/DAGMC/

CURRENT_DIR=$(pwd)

# DAGMC Install
cd $HOME
mkdir DAGMC && cd DAGMC
git clone -b $DAGMC_BRANCH $DAGMC_REPO
mkdir build && cd build
cmake ../dagmc -DBUILD_TALLY=ON -DCMAKE_INSTALL_PREFIX=$DAGMC_INSTALL_DIR -DBUILD_STATIC_LIBS=OFF -DMOAB_DIR=$MOAB_INSTALL_DIR
make -j install
rm -rf $HOME/DAGMC/dagmc $HOME/DAGMC/build

cd $CURRENT_DIR
