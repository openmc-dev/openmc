# To build with OpenMC and by default this Dockerfile builds the master branch of OpenMC.
# docker build -t openmc .

# To build with OpenMC develop branch
# docker build -t openmc_develop --build-arg openmc_branch=develop .

# To build with OpenMC and DAGMC enabled
# docker build -t openmc_dagmc --build-arg build_dagmc=on --build-arg compile_cores=4 .

# To build with OpenMC and Libmesh enabled
# docker build -t openmc_libmesh --build-arg build_libmesh=on --build-arg compile_cores=4 .

# To build with both DAGMC and Libmesh enabled
# docker build -t openmc_dagmc_libmesh --build-arg build_dagmc=on --build-arg build_libmesh=on --build-arg compile_cores=4 .

# sudo docker run image_name:tag_name or ID with no tag sudo docker run ID number


# global ARG as these ARGS are used in multiple stages
# By default one core is used to compile
ARG compile_cores=1

# By default this Dockerfile builds OpenMC without DAGMC and LIBMESH support
ARG build_dagmc=off
ARG build_libmesh=off

FROM debian:bookworm-slim AS dependencies

ARG compile_cores
ARG build_dagmc
ARG build_libmesh

# Set default value of HOME to /root
ENV HOME=/root

# Embree variables
ENV EMBREE_TAG='v4.3.1'
ENV EMBREE_REPO='https://github.com/embree/embree'
ENV EMBREE_INSTALL_DIR=$HOME/EMBREE/

# MOAB variables
ENV MOAB_TAG='5.5.1'
ENV MOAB_REPO='https://bitbucket.org/fathomteam/moab/'

# Double-Down variables
ENV DD_TAG='v1.1.0'
ENV DD_REPO='https://github.com/pshriwise/double-down'
ENV DD_INSTALL_DIR=$HOME/Double_down

# DAGMC variables
ENV DAGMC_BRANCH='v3.2.3'
ENV DAGMC_REPO='https://github.com/svalinn/DAGMC'
ENV DAGMC_INSTALL_DIR=$HOME/DAGMC/

# LIBMESH variables
ENV LIBMESH_TAG='v1.7.1'
ENV LIBMESH_REPO='https://github.com/libMesh/libmesh'
ENV LIBMESH_INSTALL_DIR=$HOME/LIBMESH

# NJOY variables
ENV NJOY_REPO='https://github.com/njoy/NJOY2016'

# Setup environment variables for Docker image
ENV LD_LIBRARY_PATH=${DAGMC_INSTALL_DIR}/lib:$LD_LIBRARY_PATH \
    OPENMC_ENDF_DATA=/root/endf-b-vii.1 \
    DEBIAN_FRONTEND=noninteractive

# Install and update dependencies from Debian package manager
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y \
        python3-pip python-is-python3 wget git build-essential cmake \
        mpich libmpich-dev libhdf5-serial-dev libhdf5-mpich-dev \
        libpng-dev python3-venv && \
    apt-get autoremove

# create virtual enviroment to avoid externally managed environment error
RUN python3 -m venv openmc_venv
ENV PATH=/openmc_venv/bin:$PATH

# Update system-provided pip
RUN pip install --upgrade pip

# Clone and install NJOY2016
RUN cd $HOME \
    && git clone --single-branch --depth 1 ${NJOY_REPO} \
    && cd NJOY2016 \
    && mkdir build \
    && cd build \
    && cmake -Dstatic=on .. \
    && make 2>/dev/null -j${compile_cores} install \
    && rm -rf $HOME/NJOY2016


RUN if [ "$build_dagmc" = "on" ]; then \
        # Install addition packages required for DAGMC
        apt-get -y install libeigen3-dev libnetcdf-dev libtbb-dev libglfw3-dev \
        && pip install --upgrade numpy \
        # Clone and install EMBREE
        && mkdir -p $HOME/EMBREE && cd $HOME/EMBREE \
        && git clone --single-branch -b ${EMBREE_TAG} --depth 1 ${EMBREE_REPO} \
        && mkdir build && cd build \
        && cmake ../embree \
                    -DCMAKE_INSTALL_PREFIX=${EMBREE_INSTALL_DIR} \
                    -DEMBREE_MAX_ISA=NONE \
                    -DEMBREE_ISA_SSE42=ON \
                    -DEMBREE_ISPC_SUPPORT=OFF \
        && make 2>/dev/null -j${compile_cores} install \
        && rm -rf ${EMBREE_INSTALL_DIR}/build ${EMBREE_INSTALL_DIR}/embree ; \
        # Clone and install MOAB
        mkdir -p $HOME/MOAB && cd $HOME/MOAB \
        && git clone  --single-branch -b ${MOAB_TAG} --depth 1 ${MOAB_REPO} \
        && mkdir build && cd build \
        && cmake ../moab -DCMAKE_BUILD_TYPE=Release \
                      -DENABLE_HDF5=ON \
                      -DENABLE_NETCDF=ON \
                      -DBUILD_SHARED_LIBS=OFF \
                      -DENABLE_FORTRAN=OFF \
                      -DENABLE_BLASLAPACK=OFF \
        && make 2>/dev/null -j${compile_cores} install \
        && cmake ../moab \
                    -DENABLE_PYMOAB=ON \
                    -DBUILD_SHARED_LIBS=ON \
        && make 2>/dev/null -j${compile_cores} install \
        && cd pymoab && bash install.sh \
        && python setup.py install \
        && python -c "import pymoab" \
        && rm -rf $HOME/MOAB ; \
        # Clone and install Double-Down
        mkdir -p $HOME/Double_down && cd $HOME/Double_down \
        && git clone --single-branch -b ${DD_TAG} --depth 1 ${DD_REPO} \
        && mkdir build && cd build \
        && cmake ../double-down -DCMAKE_INSTALL_PREFIX=${DD_INSTALL_DIR} \
                             -DMOAB_DIR=/usr/local \
                             -DEMBREE_DIR=${EMBREE_INSTALL_DIR} \
        && make 2>/dev/null -j${compile_cores} install \
        && rm -rf ${DD_INSTALL_DIR}/build ${DD_INSTALL_DIR}/double-down ; \
        # Clone and install DAGMC
        mkdir -p $HOME/DAGMC && cd $HOME/DAGMC \
        && git clone --single-branch -b ${DAGMC_BRANCH} --depth 1 ${DAGMC_REPO} \
        && mkdir build && cd build \
        && cmake ../DAGMC -DBUILD_TALLY=ON \
                       -DCMAKE_INSTALL_PREFIX=${DAGMC_INSTALL_DIR} \
                       -DMOAB_DIR=/usr/local \
                       -DDOUBLE_DOWN=ON \
                       -DDOUBLE_DOWN_DIR=${DD_INSTALL_DIR} \
                       -DCMAKE_PREFIX_PATH=${DD_INSTALL_DIR}/lib \
                       -DBUILD_STATIC_LIBS=OFF \
        && make 2>/dev/null -j${compile_cores} install \
        && rm -rf ${DAGMC_INSTALL_DIR}/DAGMC ${DAGMC_INSTALL_DIR}/build ; \
    fi


RUN if [ "$build_libmesh" = "on" ]; then \
        # Install addition packages required for LIBMESH
        apt-get -y install m4 libnetcdf-dev libpnetcdf-dev \
        # Install LIBMESH
        && mkdir -p $HOME/LIBMESH && cd $HOME/LIBMESH \
        && git clone --shallow-submodules --recurse-submodules --single-branch -b ${LIBMESH_TAG} --depth 1 ${LIBMESH_REPO} \
        && mkdir build && cd build \
        && ../libmesh/configure \
                    --prefix=${LIBMESH_INSTALL_DIR} CXX=mpicxx CC=mpicc FC=mpifort F77=mpif77 \
                    --enable-exodus \
                    --enable-mpi \
                    --enable-silent-rules \
                    --enable-unique-id \
                    --disable-eigen \
                    --disable-fortran \
                    --disable-lapack \
                    --disable-examples \
                    --disable-warnings \
                    --disable-maintainer-mode \
                    --disable-metaphysicl \
                    --with-methods="opt" \
                    --without-gdb-command \
                    --with-cxx-std-min=2014 \
        && make 2>/dev/null -j${compile_cores} install \
        && rm -rf ${LIBMESH_INSTALL_DIR}/build ${LIBMESH_INSTALL_DIR}/libmesh ; \
    fi

FROM dependencies AS build

ENV HOME=/root

ARG openmc_branch=master
ENV OPENMC_REPO='https://github.com/openmc-dev/openmc'

ARG compile_cores
ARG build_dagmc
ARG build_libmesh

ENV DAGMC_INSTALL_DIR=$HOME/DAGMC/
ENV LIBMESH_INSTALL_DIR=$HOME/LIBMESH

# clone and install openmc
RUN mkdir -p ${HOME}/OpenMC && cd ${HOME}/OpenMC \
    && git clone --shallow-submodules --recurse-submodules --single-branch -b ${openmc_branch} --depth=1 ${OPENMC_REPO} \
    && mkdir build && cd build ; \
    if [ ${build_dagmc} = "on" ] && [ ${build_libmesh} = "on" ]; then \
        cmake ../openmc \
            -DCMAKE_CXX_COMPILER=mpicxx \
            -DOPENMC_USE_MPI=on \
            -DHDF5_PREFER_PARALLEL=on \
            -DOPENMC_USE_DAGMC=on \
            -DOPENMC_USE_LIBMESH=on \
            -DCMAKE_PREFIX_PATH="${DAGMC_INSTALL_DIR};${LIBMESH_INSTALL_DIR}" ; \
    fi ; \
    if [ ${build_dagmc} = "on" ] && [ ${build_libmesh} = "off" ]; then \
        cmake ../openmc \
            -DCMAKE_CXX_COMPILER=mpicxx \
            -DOPENMC_USE_MPI=on \
            -DHDF5_PREFER_PARALLEL=on \
            -DOPENMC_USE_DAGMC=ON \
            -DCMAKE_PREFIX_PATH=${DAGMC_INSTALL_DIR} ; \
    fi ; \
    if [ ${build_dagmc} = "off" ] && [ ${build_libmesh} = "on" ]; then \
        cmake ../openmc \
            -DCMAKE_CXX_COMPILER=mpicxx \
            -DOPENMC_USE_MPI=on \
            -DHDF5_PREFER_PARALLEL=on \
            -DOPENMC_USE_LIBMESH=on \
            -DCMAKE_PREFIX_PATH=${LIBMESH_INSTALL_DIR} ; \
    fi ; \
    if [ ${build_dagmc} = "off" ] && [ ${build_libmesh} = "off" ]; then \
        cmake ../openmc \
            -DCMAKE_CXX_COMPILER=mpicxx \
            -DOPENMC_USE_MPI=on \
            -DHDF5_PREFER_PARALLEL=on ; \
    fi ; \
    make 2>/dev/null -j${compile_cores} install \
    && cd ../openmc && pip install .[test,depletion-mpi] \
    && python -c "import openmc"

FROM build AS release

ENV HOME=/root
ENV OPENMC_CROSS_SECTIONS=/root/nndc_hdf5/cross_sections.xml

# Download cross sections (NNDC and WMP) and ENDF data needed by test suite
RUN ${HOME}/OpenMC/openmc/tools/ci/download-xs.sh
