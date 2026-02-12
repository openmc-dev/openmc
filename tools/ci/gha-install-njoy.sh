#!/bin/bash
set -ex
cd $HOME
git clone --depth 1 https://github.com/njoy/NJOY2016
cd NJOY2016
mkdir build && cd build
cmake -Dstatic=on .. && make -j `nproc` 2>/dev/null && sudo make -j `nproc` install
