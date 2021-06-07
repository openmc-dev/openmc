# To build with OpenMC
# docker build -t openmc .

# To build with OpenMC and DAGMC enabled
# docker build -t openmc_dagmc --build-arg include_dagmc=true .

# To make use of multiple cores during the compile stages of the docker build
# docker build -t openmc_dagmc --build-arg compile_cores=8 .

FROM ubuntu:latest

# By default this Dockerfile builds OpenMC without dagmc
ARG include_dagmc=false

# By default one core is used to compile
ARG compile_cores=1

# Setup environment variables for Docker image
ENV CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx \
    PATH=/opt/openmc/bin:/opt/NJOY2016/build:$PATH \
    LD_LIBRARY_PATH=/opt/openmc/lib:$LD_LIBRARY_PATH \
    OPENMC_CROSS_SECTIONS=/root/nndc_hdf5/cross_sections.xml \
    OPENMC_ENDF_DATA=/root/endf-b-vii.1 \
    DEBIAN_FRONTEND=noninteractive

# Install dependencies from Debian package manager
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y \
        python3-pip python-is-python3 wget git gfortran g++ cmake \
        mpich libmpich-dev libhdf5-serial-dev libhdf5-mpich-dev \
        imagemagick && \
    apt-get autoremove

# Update system-provided pip
RUN pip3 install --upgrade pip

# Clone and install NJOY2016
RUN git clone https://github.com/njoy/NJOY2016.git /opt/NJOY2016 && \
    cd /opt/NJOY2016 && \
    mkdir build && cd build && \
    cmake -Dstatic=on .. && make 2>/dev/null && make install

# Clone and install OpenMC without DAGMC
RUN if [ "$include_dagmc" = "false" ] ; \
    then git clone --recurse-submodules https://github.com/openmc-dev/openmc.git /opt/openmc ; \
    cd /opt/openmc ; \
    mkdir -p build ; \
    cd build ; \
    cmake -Doptimize=on \
          -DHDF5_PREFER_PARALLEL=on .. ; \
    make -j"$compile_cores" ; \
    make -j"$compile_cores" install ; \
    cd ..  ; \
    pip install -e .[test] ; \
    fi

# install addition packages required for DAGMC
RUN if [ "$include_dagmc" = "true" ] ; \
    then apt-get --yes install libeigen3-dev ; \
    apt-get --yes install libnetcdf-dev ; \
    apt-get --yes install libtbb-dev ; \
    apt-get --yes install libglfw3-dev ; \
    fi

# Clone and install Embree
RUN if [ "$include_dagmc" = "true" ] ; \
    then git clone --single-branch --branch v3.12.2 https://github.com/embree/embree.git ; \
    cd embree ; \
    mkdir build ; \
    cd build ; \
    cmake .. -DCMAKE_INSTALL_PREFIX=.. \
             -DEMBREE_ISPC_SUPPORT=OFF ; \
    make -j"$compile_cores" ; \
    make -j"$compile_cores" install ; \
    fi

# Clone and install MOAB
RUN if [ "$include_dagmc" = "true" ] ; \
    then pip install --upgrade numpy cython ; \
    mkdir MOAB ; \
    cd MOAB ; \
    mkdir build ; \
    git clone  --single-branch --branch 5.2.1 https://bitbucket.org/fathomteam/moab.git ; \
    cd build ; \
    cmake ../moab -DENABLE_HDF5=ON \
                  -DENABLE_NETCDF=ON \
                  -DBUILD_SHARED_LIBS=OFF \
                  -DENABLE_FORTRAN=OFF \
                  -DENABLE_BLASLAPACK=OFF ; \
    make -j"$compile_cores" ; \
    make -j"$compile_cores" install ; \
    rm -rf * ; \
    cmake ../moab -DBUILD_SHARED_LIBS=ON \
                  -DENABLE_HDF5=ON \
                  -DENABLE_PYMOAB=ON \
                  -DENABLE_FORTRAN=OFF \
                  -DENABLE_BLASLAPACK=OFF ; \
    make -j"$compile_cores" ; \
    make -j"$compile_cores" install ; \
    fi

# Clone and install Double-Down
RUN if [ "$include_dagmc" = "true" ] ; \
    then git clone --single-branch --branch main https://github.com/pshriwise/double-down.git ; \
    cd double-down ; \
    mkdir build ; \
    cd build ; \
    cmake .. -DCMAKE_INSTALL_PREFIX=.. \
             -DMOAB_DIR=/usr/local \
             -DEMBREE_DIR=/embree ; \
    make -j"$compile_cores" ; \
    make -j"$compile_cores" install ; \
    fi

# Clone and install DAGMC
RUN if [ "$include_dagmc" = "true" ] ; \
    then mkdir DAGMC ; \
    cd DAGMC ; \
    git clone --single-branch --branch 3.2.0 https://github.com/svalinn/DAGMC.git ; \
    mkdir build ; \
    cd build ; \
    cmake ../DAGMC -DBUILD_TALLY=ON \
                   -DCMAKE_INSTALL_PREFIX=/DAGMC/ \
                   -DMOAB_DIR=/usr/local \
                   -DBUILD_STATIC_LIBS=OFF \
                   -DBUILD_STATIC_EXE=OFF ; \
    make -j"$compile_cores" install ; \
    rm -rf /DAGMC/DAGMC /DAGMC/build ; \
    fi

# Clone and install OpenMC with DAGMC
RUN if [ "$include_dagmc" = "true" ] ; \
    then git clone --recurse-submodules https://github.com/openmc-dev/openmc.git /opt/openmc ; \
    cd /opt/openmc ; \
    mkdir build ; \
    cd build ; \
    cmake -Doptimize=on \
          -Ddagmc=ON \
          -DDAGMC_DIR=/DAGMC/ \
          -DHDF5_PREFER_PARALLEL=on ..  ; \
    make -j"$compile_cores" ; \
    make -j"$compile_cores" install ; \
    cd ..  ; \
    pip install -e .[test] ; \
    fi

# Download cross sections (NNDC and WMP) and ENDF data needed by test suite
RUN /opt/openmc/tools/ci/download-xs.sh
