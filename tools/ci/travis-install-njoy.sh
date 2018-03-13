#!/bin/bash
set -ex
cd $HOME
git clone https://github.com/njoy/NJOY2016
cd NJOY2016
sed -i -e 's/5\.1/4.8/' CMakeLists.txt
mkdir build && cd build
cmake -Dstatic=on .. && make 2>/dev/null && sudo make install
