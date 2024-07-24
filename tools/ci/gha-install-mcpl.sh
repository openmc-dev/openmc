#!/bin/bash
set -ex

MCPL_BRANCH='main'
MCPL_REPO='https://github.com/mctools/mcpl'
MCPL_INSTALL_DIR=$HOME/mcpl/

cd $HOME
git clone -b $MCPL_BRANCH $MCPL_REPO
cd $MCPL_INSTALL_DIR
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$MCPL_INSTALL_DIR
make 2>/dev/null && sudo make install