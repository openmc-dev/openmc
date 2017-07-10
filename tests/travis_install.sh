#!/bin/bash

set -ev

# Build MPICH
if [[ ! -e $HOME/mpich_install/bin/mpiexec ]]; then
    wget -q http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
    tar -xzvf mpich-3.2.tar.gz >/dev/null 2>&1
    cd mpich-3.2
    ./configure --prefix=$HOME/mpich_install -q
    make -j 2 &> /dev/null
    make install &> /dev/null
    cd ..
fi

# Build PHDF5
if [[ ! -e $HOME/phdf5_install/bin/h5pfc ]]; then
    wget -q https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar.bz2
    tar -xjf hdf5-1.10.1.tar.bz2
    mv hdf5-1.10.1 phdf5-1.10.1
    cd phdf5-1.10.1
    CC=$HOME/mpich_install/bin/mpicc FC=$HOME/mpich_install/bin/mpif90 \
      ./configure --prefix=$HOME/phdf5_install -q \
      --enable-fortran --enable-parallel
    make -j 2 &> /dev/null
    make install &> /dev/null
    cd ..
fi

# Build HDF5
if [[ ! -e $HOME/hdf5_install/bin/h5fc ]]; then
    tar -xjf hdf5-1.10.1.tar.bz2
    cd hdf5-1.10.1
    CC=gcc FC=gfortran ./configure --prefix=$HOME/hdf5_install -q \
                                   --enable-fortran
    make -j 2 &> /dev/null
    make install &> /dev/null
    cd ..
fi
