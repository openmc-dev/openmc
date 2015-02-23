#!/bin/sh

set -ev

# Build HDF5 and PETSc for rest of debug tests
if [ "$TRAVIS_PULL_REQUEST" != "false" ]; then

  # Build MPICH
  wget -q http://www.mpich.org/static/downloads/3.1.3/mpich-3.1.3.tar.gz
  tar -xzvf mpich-3.1.3.tar.gz >/dev/null 2>&1
  cd mpich-3.1.3
  ./configure --prefix=$PWD/../mpich_install -q
  make -j >/dev/null 2>&1
  make install >/dev/null 2>&1
  cd ..

  # Build PHDF5
  wget -q http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.14.tar.gz
  tar -xzvf hdf5-1.8.14.tar.gz >/dev/null 2>&1
  mv hdf5-1.8.14 phdf5-1.8.14; cd phdf5-1.8.14
  CC=$PWD/../mpich_install/bin/mpicc FC=$PWD/../mpich_install/bin/mpif90 \
     ./configure \
     --prefix=$PWD/../phdf5_install -q --enable-fortran \
     --enable-fortran2003 --enable-parallel
  make -j >/dev/null 2>&1
  make install >/dev/null 2>&1
  cd ..

  # Build HDF5
  tar -xzvf hdf5-1.8.14.tar.gz >/dev/null 2>&1
  cd  hdf5-1.8.14
  CC=gcc FC=gfortran ./configure --prefix=$PWD/../hdf5_install -q \
                                 --enable-fortran --enable-fortran2003
  make -j >/dev/null 2>&1
  make install >/dev/null 2>&1
  cd ..

  # Build PETSc
  wget -q http://ftp.mcs.anl.gov/pub/petsc/release-snapshots/petsc-lite-3.5.3.tar.gz
  tar -xzvf petsc-lite-3.5.3.tar.gz >/dev/null 2>&1
  cd petsc-3.5.3
  ./configure --prefix=$PWD/../petsc_install -q --download-fblaslapack \
              --with-mpi-dir=$PWD/../mpich_install --with-share-libraries \
              --with-fortran-datatypes
  make PETSC_DIR=$PWD \
       PETSC_ARCH=arch-linux2-c-debug all >/dev/null 2>&1
  make PETSC_DIR=$PWD \
       PETSC_ARCH=arch-linux2-c-debug install >/dev/null 2>&1
  make install >/dev/null 2>&1
  cd ..

fi
