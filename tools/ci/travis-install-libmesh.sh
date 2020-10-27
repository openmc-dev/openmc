
#!/bin/bash
set -ex

CURRENT_DIR=$(pwd)

# LibMESH install
cd $HOME
mkdir LIBMESH && cd LIBMESH
git clone https://github.com/libmesh/libmesh -b v1.6.0 --recurse-submodules
mkdir build && cd build
../libmesh/configure --prefix=$HOME/LIBMESH --enable-exodus --disable-mpi
make -j4 install
rm -rf $HOME/LIBMESH/build

cd $CURRENT_DIR
