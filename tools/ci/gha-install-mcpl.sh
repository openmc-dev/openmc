#!/bin/bash
set -ex
cd $HOME
git clone https://github.com/mctools/mcpl
cd mcpl
mkdir build && cd build
cmake .. && make 2>/dev/null && sudo make install
