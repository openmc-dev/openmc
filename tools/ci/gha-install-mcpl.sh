#!/bin/bash
set -ex
cd $HOME
git clone https://github.com/mctools/mcpl
cd mcpl
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Debug .. && make 2>/dev/null && sudo make install
