#!/bin/bash
set -ex
cd $HOME
git clone https://github.com/njoy/NJOY2016
cd NJOY2016
mkdir build && cd build
cmake -Dstatic=on .. && make 2>/dev/null && sudo make install
