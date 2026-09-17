#!/bin/bash
set -ex
cd $HOME
git clone -b 2016.78 https://github.com/njoy/NJOY2016
cd NJOY2016
mkdir build && cd build
cmake -Dstatic=on .. && make 2>/dev/null && sudo make install
