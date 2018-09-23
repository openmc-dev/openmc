FROM ubuntu:latest
MAINTAINER Will Boyd <boyd.william.r@gmail.com>

# Setup environment variables for Docker image
ENV FC=/usr/bin/mpif90 CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx \
    PATH=/opt/openmc/bin:/opt/conda/bin:/opt/NJOY2016/build:$PATH \
    LD_LIBRARY_PATH=/opt/openmc/lib:$LD_LIBRARY_PATH \
    OPENMC_CROSS_SECTIONS=/opt/openmc/data/nndc_hdf5/cross_sections.xml \
    OPENMC_MULTIPOLE_LIBRARY=/opt/openmc/data/WMP_Library \
    OPENMC_ENDF_DATA=/opt/openmc/data/endf-b-vii.1

# Install dependencies from Debian package manager
RUN apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y wget git emacs && \
    apt-get install -y gfortran g++ cmake && \
    apt-get install -y mpich libmpich-dev && \
    apt-get install -y libhdf5-serial-dev libhdf5-mpich-dev && \
    apt-get autoremove

# Download Miniconda3 and install Python dependencies
RUN wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm -rf ~/miniconda.sh
RUN pip install --upgrade pip && \
    pip install six numpy scipy pandas h5py matplotlib uncertainties lxml cython jupyter pytest && \
    pip install --no-binary=mpi4py mpi4py

# Clone and install NJOY2016
RUN git clone https://github.com/njoy/NJOY2016 /opt/NJOY2016 && \
    cd /opt/NJOY2016 && \
    mkdir build && cd build && \
    cmake -Dstatic=on .. && make 2>/dev/null && make install

# Clone and install OpenMC
RUN git clone https://github.com/openmc-dev/openmc.git /opt/openmc && \
    cd /opt/openmc && mkdir -p build && cd build && \
    cmake  -Doptimize=on -DHDF5_PREFER_PARALLEL=on .. && \
    make && make install && \
    cd .. && pip install -e .

# Download cross sections (NNDC and WMP) and ENDF data needed by test suite
RUN mkdir /opt/openmc/data
RUN wget https://anl.box.com/shared/static/a0eflty17atnpd0pp7460exagr3nuhm7.xz -O - | tar -C /opt/openmc/data -xvJ
RUN wget https://anl.box.com/shared/static/4kd2gxnf4gtk4w1c8eua5fsua22kvgjb.xz -O - | tar -C /opt/openmc/data -xvJ
RUN wget https://github.com/mit-crpg/WMP_Library/releases/download/v1.0/WMP_Library_v1.0.tar.gz && \
   tar -C /opt/openmc/data -xzvf WMP_Library_v1.0.tar.gz && \
   rm WMP_Library_v1.0.tar.gz