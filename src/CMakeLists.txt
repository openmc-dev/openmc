cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
project(openmc Fortran)

#===============================================================================
# Command line options
#===============================================================================

option(openmp   "Enable shared-memory parallelism with OpenMP"   OFF)
option(profile  "Compile with profiling flags"                   OFF)
option(hdf5     "Enable HDF5 binary output"                      OFF)
option(petsc    "Enable PETSC for use in CMFD acceleration"      OFF)
option(debug    "Compile with debug flags"                       OFF)
option(optimize "Turn on all compiler optimization flags"        OFF)

#===============================================================================
# MPI for distributed-memory parallelism
#===============================================================================

if($ENV{FC} MATCHES "mpi.*")
  find_package(MPI REQUIRED)
  add_definitions(-DMPI)
  include_directories(BEFORE "${MPI_Fortran_INCLUDE_PATH}")
  set(f90flags  ${MPI_Fortran_COMPILE_FLAGS} ${f90flags})
  set(libraries ${MPI_Fortran_LIBRARIES} ${libraries})
endif()

#===============================================================================
# Set compile/link flags based on which compiler is being used
#===============================================================================

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  # GNU Fortran compiler options
  set(f90flags  "-cpp -std=f2008 -fbacktrace")
  if(debug)
    set(f90flags  "-g -Wall -pedantic -fbounds-check -ffpe-trap=invalid,overflow,underflow ${f90flags}")
    set(ldflags   "-g")
  endif()
  if(profile)
    set(f90flags  "-pg ${f90flags}")
    set(ldflags   "-pg ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-O3 ${f90flags}")
  endif()
  if(openmp)
    set(f90flags  "-fopenmp ${f90flags}")
    set(ldflags   "-fopenmp ${ldflags}")
    add_definitions(-DOPENMP)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  # Intel Fortran compiler options
  set(f90flags  "-fpp -warn -assume byterecl -traceback")
  if(debug)
    set(f90flags "-g -ftrapuv -fp-stack-check -check all -fpe0 ${f90flags}")
    set(ldflags  "-g")
  endif()
  if(profile)
    set(f90flags  "-pg ${f90flags}")
    set(ldflags   "-pg ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-O3 ${f90flags}")
  endif()
  if(openmp)
    set(f90flags  "-openmp ${f90flags}")
    set(ldflags   "-openmp ${ldflags}")
    add_definitions(-DOPENMP)
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  # PGI Fortran compiler options
  set(f90flags  "-Mpreprocess -Minform=inform -traceback")
  add_definitions(-DNO_F2008)
  if(debug)
    set(f90flags  "-g -Mbounds -Mchkptr -Mchkstk ${f90flags}")
    set(ldflags   "-g")
  endif()
  if(profile)
    set(f90flags  "-pg ${f90flags}")
    set(ldflags   "-pg ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-fast -Mipa ${f90flags}")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "XL")
  # IBM XL compiler options
  set(f90flags  "-WF,-DNO_F2008 -O2")
  if(debug)
    set(f90flags  "-g -C -qflag=i:i -u")
    set(ldflags   "-g")
  endif()
  if(profile)
    set(f90flags  "-p ${f90flags}")
    set(ldflags   "-p ${ldflags}")
  endif()
  if(optimize)
    set(f90flags  "-O3 ${f90flags}")
  endif()
  if(openmp)
    set(f90flags  "-qsmp=omp -WF,-DOPENMP ${f90flags}")
    set(ldflags   "-qsmp=omp ${ldflags}")
  endif()

elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Cray")
  # Cray Fortran compiler options
  set(f90flags  "-e Z -m 0")
  if(debug)
    set(f90flags  "-g -R abcnsp -O0 ${f90flags}")
    set(ldflags   "-g")
  endif()

endif()

#===============================================================================
# HDF5 for binary output
#===============================================================================

if(hdf5)
  find_package(HDF5 COMPONENTS Fortran Fortran_HL REQUIRED)
  add_definitions(-DHDF5)
  include_directories(BEFORE ${HDF5_INCLUDE_DIRS})
  set(libraries ${HDF5_LIBRARIES} ${HDF5_Fortran_LIBRARIES} ${HDF5_Fortran_HL_LIBRARIES} ${libraries})
endif()

#===============================================================================
# PETSc for CMFD functionality
#===============================================================================

#===============================================================================
# FoX Fortran XML Library
#===============================================================================

file(GLOB_RECURSE source_fox xml/*.F90)
add_library(fox STATIC ${source_fox})

#===============================================================================
# Build OpenMC executable
#===============================================================================

set(program "openmc")
file(GLOB source *.F90)
add_executable(${program} ${source})
target_link_libraries(${program} ${libraries} fox)
set_target_properties(${program} PROPERTIES
  COMPILE_FLAGS "${f90flags}"
  COMPILE_DEFINITIONS "${definitions}")
