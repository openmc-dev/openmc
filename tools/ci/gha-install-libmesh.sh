
#!/bin/bash
set -ex

# libMESH install
pushd $HOME
mkdir LIBMESH && cd LIBMESH
git clone https://github.com/libmesh/libmesh -b v1.7.1 --recurse-submodules
mkdir build && cd build
export METHODS="opt"

if [[ $MPI == 'y' ]]; then
    ../libmesh/configure --prefix=$HOME/LIBMESH CXX=mpicxx CC=mpicc FC=mpifort F77=mpif77 \
        --enable-exodus --disable-netcdf-4 --disable-eigen --disable-lapack
else
    ../libmesh/configure --prefix=$HOME/LIBMESH --enable-exodus --disable-netcdf-4 \
    --disable-eigen --disable-lapack --disable-mpi
fi
make -j4 install
export LIBMESH_PC=$HOME/LIBMESH/lib/pkgconfig/
#rm -rf $HOME/LIBMESH/build

popd
