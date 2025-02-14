#!/bin/bash
set -ex
cd $HOME

# NCrystal Variables
NCRYSTAL_BRANCH='v4.0.0'
NCRYSTAL_REPO='https://github.com/mctools/ncrystal'
NCRYSTAL_INSTALL_DIR=$HOME/NCRYSTAL/

CURRENT_DIR=$(pwd)

# NCrystal Install
cd $HOME
git clone -b $NCRYSTAL_BRANCH $NCRYSTAL_REPO ncrystal
cd ncrystal
mkdir build && cd build
cmake .. \
    -DBUILD_SHARED_LIBS=ON \
    -DNCRYSTAL_NOTOUCH_CMAKE_BUILD_TYPE=ON \
    -DNCRYSTAL_MODIFY_RPATH=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DNCRYSTAL_ENABLE_EXAMPLES=OFF \
    -DNCRYSTAL_ENABLE_SETUPSH=OFF \
    -DNCRYSTAL_ENABLE_DATA=EMBED \
    -DCMAKE_INSTALL_PREFIX="$NCRYSTAL_INSTALL_DIR"
make -j4 && make install
rm -rf $HOME/ncrystal
cd $CURRENT_DIR
