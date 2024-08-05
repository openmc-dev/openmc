#!/bin/bash
set -ex

# MCPL variables
MCPL_REPO='https://github.com/mctools/mcpl'

cd $HOME
git clone $MCPL_REPO
cd mcpl
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$MCPL_INSTALL_DIR &&
make 2>/dev/null && sudo make install