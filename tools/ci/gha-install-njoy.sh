#!/bin/bash
set -ex

# NJOY variables
NJOY_BRANCH='main'
NJOY_REPO='https://github.com/njoy/NJOY2016'
NJOY_INSTALL_DIR=$HOME/NJOY2016/

cd $HOME
git clone -b $NJOY_BRANCH $NJOY_REPO
cd $NJOY_INSTALL_DIR
mkdir build && cd build
cmake -Dstatic=on -DCMAKE_INSTALL_PREFIX=$NJOY_INSTALL_DIR .. 
make 2>/dev/null && sudo make install