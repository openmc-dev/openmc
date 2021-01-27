
#!/bin/bash
set -ex

CURRENT_DIR=$(pwd)

# libMESH install
cd $HOME
mkdir LIBMESH && cd LIBMESH
git clone https://github.com/libmesh/libmesh -b v1.6.0 --recurse-submodules
mkdir build && cd build
export METHODS="opt"


if [[ $MPI == 'y' ]]; then
../libmesh/configure --prefix=$HOME/LIBMESH CXX=mpicxx.mpich CC=mpicc.mpich FC=mpifort.mpich F77=mpif77.mpich \
--enable-exodus --disable-netcdf-4 --disable-eigen --disable-lapack
else
../libmesh/configure --prefix=$HOME/LIBMESH CXX=mpicxx.mpich CC=mpicc.mpich FC=mpifort.mpich F77=mpif77.mpich \
--enable-exodus --disable-netcdf-4 --disable-eigen --disable-lapack
fi
make -j4 install
export LIBMESH_PC=$HOME/LIBMESH/lib/pkgconfig/
rm -rf $HOME/LIBMESH/build

cd $CURRENT_DIR
