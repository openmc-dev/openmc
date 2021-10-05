# To build with OpenMC and by default this Dockerfile builds the master branch of OpenMC.
# docker build -t openmc .

# To build with OpenMC develop branch
# docker build -t openmc_develop --build-arg OPENMC_BRANCH=develop .

# To build with OpenMC and DAGMC enabled
# docker build -t openmc_dagmc --build-arg include_dagmc=true --build-arg compile_cores=4 .

# To build with OpenMC and Libmesh support
# docker build -t openmc_libmesh --build-arg include_libmesh=true --build-arg compile_cores=4 .

# sudo docker run image_name:tag_name or ID with no tag sudo docker run ID number

FROM ubuntu:latest

# By default this Dockerfile builds OpenMC without DAGMC and LIBMESH support
ARG include_dagmc=false
ARG include_libmesh=false

# By default one core is used to compile
ARG compile_cores=2

# OpenMC variables
ARG OPENMC_BRANCH=master
ENV OPENMC_REPO='https://github.com/openmc-dev/openmc'

# Embree variables
ENV EMBREE_TAG='v3.12.2'
ENV EMBREE_REPO='https://github.com/embree/embree'
ENV EMBREE_INSTALL_DIR=$HOME/EMBREE/

# MOAB variables
ENV MOAB_BRANCH='master'
ENV MOAB_REPO='https://bitbucket.org/fathomteam/moab/'

# Double-Down variables
ENV DD_BRANCH='main'
ENV DD_REPO=' https://github.com/pshriwise/double-down'
ENV DD_INSTALL_DIR=$HOME/Double_down

# DAGMC variables
ENV DAGMC_TAG='3.2.0'
ENV DAGMC_REPO='https://github.com/svalinn/DAGMC'
ENV DAGMC_INSTALL_DIR=$HOME/DAGMC/

# LIBMESH variables
ENV LIBMESH_TAG='v1.6.0'
ENV LIBMESH_REPO='https://github.com/libMesh/libmesh'
ENV LIBMESH_INSTALL_DIR=$HOME/LIBMESH

# Setup environment variables for Docker image
ENV CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx \
    LD_LIBRARY_PATH=${DAGMC_INSTALL_DIR}/lib:$LD_LIBRARY_PATH \
    OPENMC_CROSS_SECTIONS=/root/nndc_hdf5/cross_sections.xml \
    OPENMC_ENDF_DATA=/root/endf-b-vii.1 \
    DEBIAN_FRONTEND=noninteractive

# Install and update dependencies from Debian package manager
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y \
        python3-pip python-is-python3 wget git build-essential cmake \
        mpich libmpich-dev libhdf5-serial-dev libhdf5-mpich-dev \
        imagemagick && \
    apt-get autoremove -y

# Update system-provided pip
RUN pip install --upgrade pip

# Clone and install NJOY2016
RUN cd $HOME && git clone --depth 1 https://github.com/njoy/NJOY2016.git && \
    cd NJOY2016 && mkdir build && cd build && \
    cmake -Dstatic=on .. && make 2>/dev/null -j${compile_cores} install ; \
    rm -rf $HOME/NJOY2016

RUN if [ "$include_dagmc" = "true" ]; then \
        # Install addition packages required for DAGMC
        apt-get -y install libeigen3-dev libnetcdf-dev libtbb-dev libglfw3-dev ; \
        pip install --upgrade numpy cython ; \
        # Clone and install EMBREE
        mkdir -p $HOME/EMBREE && cd $HOME/EMBREE; \
        git clone --single-branch -b ${EMBREE_TAG} --depth 1 ${EMBREE_REPO} ; \
        mkdir build && cd build; \
        cmake ../embree \
                    -DCMAKE_INSTALL_PREFIX=${EMBREE_INSTALL_DIR} \
                    -DEMBREE_ISPC_SUPPORT=OFF ; \
        make 2>/dev/null -j${compile_cores} install ; \
        rm -rf ${EMBREE_INSTALL_DIR}/build ${EMBREE_INSTALL_DIR}/embree ; \
        # Clone and install MOAB
        mkdir -p $HOME/MOAB && cd $HOME/MOAB ; \
        git clone  --single-branch -b ${MOAB_BRANCH} --depth 1 ${MOAB_REPO} ; \
        mkdir build && cd build ; \
        cmake ../moab \
                    -DENABLE_HDF5=ON \
                    -DENABLE_NETCDF=ON \
                    -DBUILD_SHARED_LIBS=OFF \
                    -DENABLE_FORTRAN=OFF \
                    -DENABLE_BLASLAPACK=OFF ; \
        make 2>/dev/null -j${compile_cores} install ; \
        cmake ../moab \
                    -DENABLE_PYMOAB=ON \
                    -DBUILD_SHARED_LIBS=ON \
        make 2>/dev/null -j${compile_cores} install ; \
        cd pymoab && bash install.sh; \
        python setup.py install; \
        rm -rf $HOME/MOAB ; \
        # Clone and install Double-Down
        mkdir -p $HOME/Double_down && cd $HOME/Double_down; \
        git clone --single-branch -b ${DD_BRANCH} --depth 1 ${DD_REPO} ; \
        mkdir build && cd build; \
        cmake ../double-down \
                    -DCMAKE_INSTALL_PREFIX=${DD_INSTALL_DIR} \
                    -DMOAB_DIR=/usr/local \
                    -DEMBREE_DIR=${EMBREE_INSTALL_DIR} ; \
        make 2>/dev/null -j${compile_cores} install ; \
        rm -rf ${DD_INSTALL_DIR}/build && ${DD_INSTALL_DIR}/double-down; \
        # Clone and install DAGMC
        mkdir -p $HOME/DAGMC && cd $HOME/DAGMC; \
        git clone --single-branch -b ${DAGMC_TAG} --depth 1 ${DAGMC_REPO} ; \
        mkdir build && cd build ; \
        cmake ../DAGMC \
                    -DBUILD_TALLY=ON \
                    -DCMAKE_INSTALL_PREFIX=${DAGMC_INSTALL_DIR} \
                    -DMOAB_DIR=/usr/local \
                    -DDOUBLE_DOWN=ON \
                    -DDOUBLE_DOWN_DIR=${DD_INSTALL_DIR} \
                    -DCMAKE_PREFIX_PATH=${DD_INSTALL_DIR}/lib \
                    -DBUILD_STATIC_LIBS=OFF \
                    -DBUILD_STATIC_EXE=OFF ; \
        make 2>/dev/null -j${compile_cores} install ; \
        rm -rf ${DAGMC_INSTALL_DIR}/DAGMC ${DAGMC_INSTALL_DIR}/build ; \
        # Clone and install OpenMC with DAGMC support
        mkdir -p $HOME/OpenMC && cd $HOME/OpenMC; \
        git clone --shallow-submodules --recurse-submodules -b ${OPENMC_BRANCH} --depth=1 ${OPENMC_REPO} ; \
        mkdir build && cd build ; \
        cmake ../openmc \
                    -Doptimize=on \
                    -Ddagmc=ON \
                    -Ddebug=on \
                    -DDAGMC_DIR=${DAGMC_INSTALL_DIR} \
                    -DHDF5_PREFER_PARALLEL=on ; \
        make 2>/dev/null -j${compile_cores} install ; \
        cd ../openmc && pip install -e .[test] ; \
    elif [ "$include_libmesh" = "true" ]; then \
        # Install addition packages required for LIBMESH
        apt-get -y install m4 libnetcdf-dev libpnetcdf-dev ; \
        # Install LIBMESH
        mkdir -p $HOME/LIBMESH && cd $HOME/LIBMESH; \
        git clone --shallow-submodules --recurse-submodules --single-branch -b ${LIBMESH_TAG} --depth 1 ${LIBMESH_REPO} ; \
        mkdir build && cd build; \
        ../libmesh/configure \
                    --prefix=${LIBMESH_INSTALL_DIR} CXX=mpicxx CC=mpicc FC=mpifort F77=mpif77 \
                    --enable-exodus \
                    --disable-netcdf-4 \
                    --disable-eigen \
                    --disable-fortran \
                    --disable-lapack \
                    --enable-hdf5 \
                    --with-hdf5=/usr/lib/x86_64-linux-gnu/hdf5/serial/ \
                    --enable-mpi; \
        make 2>/dev/null -j${compile_cores} install ; \
        rm -rf ${LIBMESH_INSTALL_DIR}/build ${LIBMESH_INSTALL_DIR}/libmesh ; \
        # Clone and install OpenMC with LIBMESH support
        mkdir -p $HOME/OpenMC && cd $HOME/OpenMC; \
        git clone --shallow-submodules --recurse-submodules -b ${OPENMC_BRANCH} --depth=1 ${OPENMC_REPO} ; \
        mkdir build && cd build ; \
        cmake ../openmc \
                    -Doptimize=on \
                    -Dlibmesh=on \
                    -Ddebug=on \
                    -DCMAKE_PREFIX_PATH=${LIBMESH_INSTALL_DIR} \
                    -DHDF5_PREFER_PARALLEL=on ; \
        make 2>/dev/null -j${compile_cores} install ; \
        cd ../openmc && pip install -e .[test] ; \
    else \
        mkdir -p $HOME/OpenMC && cd $HOME/OpenMC; \
        git clone --shallow-submodules --recurse-submodules -b ${OPENMC_BRANCH} --depth=1 ${OPENMC_REPO} ; \
        mkdir build && cd build; \
        cmake -Doptimize=on -DHDF5_PREFER_PARALLEL=on ../openmc ; \
        make 2>/dev/null -j${compile_cores} install ; \
        cd ../openmc && pip install -e .[test] ; \
    fi ;

# Download cross sections (NNDC and WMP) and ENDF data needed by test suite
RUN $HOME/OpenMC/openmc/tools/ci/download-xs.sh
