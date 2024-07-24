
#!/bin/bash
set -ex

# MOAB Variables
MOAB_BRANCH='Version5.1.0'
MOAB_REPO='https://bitbucket.org/fathomteam/moab/'
MOAB_INSTALL_DIR=$HOME/MOAB/

# DAGMC Variables
DAGMC_BRANCH='develop'
DAGMC_REPO='https://github.com/svalinn/dagmc'
DAGMC_INSTALL_DIR=$HOME/DAGMC/

CURRENT_DIR=$(pwd)

# MOAB Install
cd $HOME
mkdir MOAB && cd MOAB
git clone -b $MOAB_BRANCH $MOAB_REPO
mkdir build && cd build
cmake ../moab -DENABLE_HDF5=ON -DENABLE_NETCDF=ON -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$MOAB_INSTALL_DIR -DENABLE_BLASLAPACK=OFF
make -j && make -j install
#rm -rf $HOME/MOAB/moab $HOME/MOAB/build

# DAGMC Install
cd $HOME
mkdir DAGMC && cd DAGMC
git clone -b $DAGMC_BRANCH $DAGMC_REPO
mkdir build && cd build
cmake ../dagmc -DBUILD_TALLY=ON -DCMAKE_INSTALL_PREFIX=$DAGMC_INSTALL_DIR -DBUILD_STATIC_LIBS=OFF -DMOAB_DIR=$MOAB_INSTALL_DIR
make -j install
#rm -rf $HOME/DAGMC/dagmc $HOME/DAGMC/build

cd $CURRENT_DIR
