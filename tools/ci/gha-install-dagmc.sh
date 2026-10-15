
#!/bin/bash
set -ex

# DAGMC Variables
DAGMC_BRANCH='develop'
DAGMC_REPO='https://github.com/svalinn/dagmc'
DAGMC_INSTALL_DIR=$HOME/DAGMC/

# DAGMC Install
pushd $HOME
mkdir -p DAGMC
cd DAGMC
git clone -b $DAGMC_BRANCH $DAGMC_REPO dagmc
mkdir build
cd build
cmake ../dagmc \
    -DBUILD_TALLY=ON \
    -DCMAKE_INSTALL_PREFIX=$DAGMC_INSTALL_DIR \
    -DBUILD_STATIC_LIBS=OFF \
    -DMOAB_DIR=$HOME/MOAB
make -j4 install
rm -rf $HOME/DAGMC/dagmc $HOME/DAGMC/build
popd
