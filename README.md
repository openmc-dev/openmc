# TPMS branch specific (INL developers)

<img src="https://media1.tenor.com/m/98dCgWDs4rsAAAAC/epic-sax-guy-sax.gif" width="100px" align="center">

## Installation

##### create a new conda env:
```
conda create -n openmc-TPMS moose-dev=2024.10.01=mpich
conda deactivate
conda activate openmc-TPMS
conda install boost=1.85.0
```
##### clone the project to your local repo:
```
git clone git@github.inl.gov:paul-ferney/openmc.git
cd openmc
git checkout TPMS
```
##### install openmc python API in dev mode:
```
python -m pip install -e .[test]
```
##### check that your compilers are from your conda env:
```
which $CC && which $CXX
>>> /opt/anaconda3/envs/openmc-TPMS/bin/mpicc
>>> /opt/anaconda3/envs/openmc-TPMS/bin/mpicxx
```
##### create a build directory and make your project:
```
mkdir ./build && cd ./build
cmake .. && make -j 12
cd ..
```
##### test openmc:
```
cd tests/TPMS/base_no_tpms
python base_no_tpms.py
../../../build/bin/openmc
```

## Trouble shooting
### Shared libraries conflict
When running cmake, you might get: 
Warning Cannot generate a safe runtime search path for target libopenmc because files in some directories may conflict with libraries in implicit directories:
```
CMake Warning at CMakeLists.txt:463 (add_library):
  Cannot generate a safe runtime search path for target libopenmc because
  files in some directories may conflict with libraries in implicit
  directories:

    runtime library [libpng.dylib] in /opt/anaconda3/envs/openmc-TPMS/lib may be hidden by files in:
      /opt/homebrew/lib
    runtime library [libpugixml.1.dylib] in /opt/anaconda3/envs/openmc-TPMS/lib may be hidden by files in:
      /opt/homebrew/lib

  Some of these libraries may not be found correctly.
```
There is a conflict between your conda env shared libraries and shared libraries installed somewhere else. You want cmake to strictly use your conda libs. Run instead:
```
cmake .. -DCMAKE_PREFIX_PATH=/opt/anaconda3/envs/openmc-TPMS/ -DOPENMC_USE_MPI=ON -DCMAKE_FIND_USE_CMAKE_SYSTEM_PATH=FALSE -DCMAKE_FIND_USE_SYSTEM_ENVIRONMENT_PATH=FALSE 
```

# OpenMC Monte Carlo Particle Transport Code

[![License](https://img.shields.io/badge/license-MIT-green)](https://docs.openmc.org/en/latest/license.html)
[![GitHub Actions build status (Linux)](https://github.com/openmc-dev/openmc/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/openmc-dev/openmc/actions/workflows/ci.yml)
[![Code Coverage](https://coveralls.io/repos/github/openmc-dev/openmc/badge.svg?branch=develop)](https://coveralls.io/github/openmc-dev/openmc?branch=develop)
[![dockerhub-publish-develop-dagmc](https://github.com/openmc-dev/openmc/workflows/dockerhub-publish-develop-dagmc/badge.svg)](https://github.com/openmc-dev/openmc/actions?query=workflow%3Adockerhub-publish-develop-dagmc)
[![dockerhub-publish-develop](https://github.com/openmc-dev/openmc/workflows/dockerhub-publish-develop/badge.svg)](https://github.com/openmc-dev/openmc/actions?query=workflow%3Adockerhub-publish-develop)
[![conda-pacakge](https://anaconda.org/conda-forge/openmc/badges/version.svg)](https://anaconda.org/conda-forge/openmc)

The OpenMC project aims to provide a fully-featured Monte Carlo particle
transport code based on modern methods. It is a constructive solid geometry,
continuous-energy transport code that uses HDF5 format cross sections. The
project started under the Computational Reactor Physics Group at MIT.

Complete documentation on the usage of OpenMC is hosted on Read the Docs (both
for the [latest release](https://docs.openmc.org/en/stable/) and
[developmental](https://docs.openmc.org/en/latest/) version). If you are
interested in the project, or would like to help and contribute, please get in
touch on the OpenMC [discussion forum](https://openmc.discourse.group/).

## Installation

Detailed [installation
instructions](https://docs.openmc.org/en/stable/usersguide/install.html)
can be found in the User's Guide.

## Citing

If you use OpenMC in your research, please consider giving proper attribution by
citing the following publication:

- Paul K. Romano, Nicholas E. Horelik, Bryan R. Herman, Adam G. Nelson, Benoit
  Forget, and Kord Smith, "[OpenMC: A State-of-the-Art Monte Carlo Code for
  Research and Development](https://doi.org/10.1016/j.anucene.2014.07.048),"
  *Ann. Nucl. Energy*, **82**, 90--97 (2015).

## Troubleshooting

If you run into problems compiling, installing, or running OpenMC, first check
the [Troubleshooting
section](https://docs.openmc.org/en/stable/usersguide/troubleshoot.html) in the
User's Guide. If you are not able to find a solution to your problem there,
please post to the [discussion forum](https://openmc.discourse.group/).

## Reporting Bugs

OpenMC is hosted on GitHub and all bugs are reported and tracked through the
[Issues](https://github.com/openmc-dev/openmc/issues) feature on GitHub.
However, GitHub Issues should not be used for common troubleshooting purposes.
If you are having trouble installing the code or getting your model to run
properly, you should first send a message to the [discussion
forum](https://openmc.discourse.group/). If it turns out your issue really is a
bug in the code, an issue will then be created on GitHub. If you want to request
that a feature be added to the code, you may create an Issue on github.

## License

OpenMC is distributed under the MIT/X
[license](https://docs.openmc.org/en/stable/license.html).
