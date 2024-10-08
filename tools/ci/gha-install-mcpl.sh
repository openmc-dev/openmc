#!/bin/bash
set -ex

# MCPL variables
MCPL_REPO='https://github.com/mctools/mcpl'

cd $HOME
git clone $MCPL_REPO
cd mcpl
mkdir build && cd build
cmake .. && make 2>/dev/null && sudo make install